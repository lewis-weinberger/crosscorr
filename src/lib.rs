extern crate fftw;
extern crate ron;
extern crate serde;
extern crate byteorder;
extern crate num_complex;

use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::f64::consts::PI;
use fftw::array::AlignedVec;
use fftw::types::*;
use fftw::plan::*;
use byteorder::{ReadBytesExt, NativeEndian};
use ron::de::from_reader;
use serde::Deserialize;
use num_complex::Complex;

/// A struct containing the configuration information to run the program, read
/// at runtime from a RON file.
/// 
/// # Examples
///
/// ```
/// let config = Config {
///     grid1_filename: String::from("/path/to/grid1"),
///     grid2_filename: String::from("/path/to/grid2"),
///     output_filename: String::from("/path/to/output"),
///     ngrid: 2048,
///     boxsize: 160.0,
///     }
/// ```
#[derive(Debug, Deserialize)]
pub struct Config {
    pub grid1_filename: String,
    pub grid2_filename: String,
    pub output_filename: String,
    pub ngrid: u32,
    pub boxsize: f32,
}

impl Config {
    /// Reads the configuration file passed at runtime.
    ///
    /// # Examples
    ///
    /// ```
    /// ```
    pub fn new(mut args: std::env::Args) -> Result<Config, &'static str> {
        args.next();
        
        // Match command-line argument for configuration filename
        let config_filename = match args.next() {
            Some(arg) => arg,
            None => return Err("Incorrect command-line argument."),
        };

        // Open configuration file
        println!("\nReading configuration file: {}", config_filename);
	let f = match File::open(&config_filename) {
            Ok(file) => file,
            Err(_) => return Err("Unable to open configuration file."),
        };
	
        // Decode RON format of configuration file
	let config: Config = match from_reader(f) {
	    Ok(x) => x,
	    Err(_) => return Err("Unable to read configuration from file.")
	};

        // Print configuration
        println!("\ngrid1 path:  {}", config.grid1_filename);
        println!("grid1 path:  {}", config.grid2_filename);
        println!("output path: {}", config.output_filename);
        println!("ngrid:       {} cells on a side", config.ngrid);
        println!("boxsize:     {} cMpc/h", config.boxsize);

        Ok(config)
    }
}

/// A struct containing the final output vectors.
/// 
/// # Examples
///
/// ```
/// ```
#[derive(Debug)]
pub struct Output {
    pub w: Vec<f64>,
    pub pow_spec: Vec<f64>,
    pub deltasqk: Vec<f64>,
    pub iweights: Vec<i64>,
}

impl Output {
    /// Saves the power spectrum to a formatted txt file.
    ///
    /// # Examples
    ///
    /// ```
    /// ```
    pub fn save_result(&self, config: &Config) -> Result<(), &'static str> {       
        println!("\nSaving results to: {}", &config.output_filename);

        // Open output file
        let mut f = match File::create(&config.output_filename) {
            Ok(file) => file,
            Err(_) => return Err("Unable to open output file!"),
        };
        write!(f,"# w pow_spec deltasqk iweights\n").unwrap();

        let nhalf: usize = (config.ngrid / 2) as usize;
        for n in 0..nhalf {
            write!(f, "{} {} {} {}\n", self.w[n], self.pow_spec[n], self.deltasqk[n], self.iweights[n]).unwrap();
        }
        
        Ok(())
    }
}

/// Loads a grid stored at `filename` (in a custom binary format) into an 
/// `fftw::array::AlignedVec` object. This custom format stores the 3D grid as
/// a 1D array of values. The data should be stored as deviations from the mean,
/// i.e. delta = (x - mean(x)) / mean(x).
///
/// # Examples
///
/// ```
/// ```
pub fn load_grid(config: &Config, num: usize) -> Result<AlignedVec<c64>, &'static str> {
    let filename = match num {
        1 => &config.grid1_filename,
        2 => &config.grid2_filename,
        _ => return Err("Need to load either grid 1 or 2!")
    };
    println!("\nOpening grid from file: {}", filename);
    let ngrid: usize = config.ngrid as usize;

    // Allocate AlignedVec array to hold grid
    let ngrid3 = ngrid*ngrid*ngrid;
    let mut grid = AlignedVec::new(ngrid3);

    // Open binary file
    let f = match File::open(filename) {
        Ok(file) => file,
        Err(_) => return Err("Unable to open grid file!"),
    };
    let mut buf_reader = BufReader::new(f); 

    // Read in array from binary file
    for i in 0..ngrid3 {
        let cell = match buf_reader.read_f64::<NativeEndian>() {
            Ok(val) => val,
            Err(_) => return Err("Problem reading values from file!"),
        };
        grid[i] = c64::new(cell, 0.0);
    }
    println!("Successfully read {} cells!", ngrid3);

    Ok(grid)
}

/// Calculates the cross power spectrum of the given 3D grids (note if the same
/// grid is given twice then this is the auto power spectrum).
///
/// # Examples
///
/// ```
/// ```
pub fn correlate(config: &Config, mut grid1: AlignedVec<c64>, mut grid2: AlignedVec<c64>) -> Result<Output, &'static str> {
    let ngrid: usize = config.ngrid as usize;
    let boxsize: f64 = config.boxsize as f64;
    
    // Create FFTW plan
    let shape = [ngrid, ngrid, ngrid];
    let mut plan: C2CPlan64 = match C2CPlan::aligned(&shape[..], Sign::Forward, Flag::Estimate) {
        Ok(p) => p,
        Err(_) => return Err("Unable to create FFTW plan."),
    };

    // Perform FFT on grids
    let ngrid3 = ngrid * ngrid * ngrid;
    
    let mut out1 = AlignedVec::new(ngrid3);
    match plan.c2c(&mut grid1, &mut out1) {
        Ok(_) => (),
        Err(_) => return Err("Failed to FFT grid1.")
    };
    
    let mut out2 = AlignedVec::new(ngrid3);
    match plan.c2c(&mut grid2, &mut out2) {
        Ok(_) => (),
        Err(_) => return Err("Failed to FFT grid2.")
    };

    // Sanity prints
    println!("\nFFTs performed! Sanity check:");
    for n in 0..10 {
        println!("grid1[{}] = {} + i{}, out1[{}] = {} + i{}", n, grid1[n].re, grid1[n].im, n, out1[n].re, out1[n].im);
        println!("grid2[{}] = {} + i{}, out2[{}] = {} + i{}", n, grid2[n].re, grid2[n].im, n, out2[n].re, out2[n].im);
    }

    // Calculate power spectrum
    let kf: f64 = 2.0 * PI / boxsize;
    let coeff: f64 = (boxsize / (2.0 * PI)).powf(2.0);
    let nhalf: usize = ngrid / 2;

    let mut w: Vec<f64> = Vec::with_capacity(ngrid);
    for i in 0..nhalf {
        w[i] = kf * (i as f64);
    }
    for i in nhalf..ngrid {
        w[i] = kf * ((i - ngrid) as f64);
    }

    let mut pow_spec: Vec<f64> = Vec::with_capacity(ngrid);
    let mut iweights: Vec<i64> = Vec::with_capacity(ngrid);
    for n in 0..ngrid {
        pow_spec[n] = 0.0;
        iweights[n] = 0;
    }
    
    for i in 0..ngrid {
        let iper = if i > nhalf { ngrid - i } else { i };
        for j in 0..ngrid {
            let jper = if j > nhalf { ngrid - j } else { j };
            for k in 0..ngrid {
                let kper = if k > nhalf { ngrid - k } else { k };
                let r: f64 = (iper*iper + jper*jper + kper*kper) as f64;
                let m: usize = (0.5 + r.sqrt()) as usize;
                iweights[m] += 1;
        
                let g = w[i]*w[i] + w[j]*w[j] + w[k]*w[k];
                if g != 0.0 {
                    let scale: usize = (0.5 + (g*coeff).sqrt()) as usize;
                    let index: usize = k + ngrid*(j + ngrid*i);
                    let contrib: Complex<f64> = out1[index] + out2[index].conj() + out1[index].conj() + out2[index];
                    pow_spec[scale] += contrib.re / 2.0; 
                }
            }
        }
    }
    println!("Power spectrum calculated. Normalising...");
   
    // Normalise power spectrum
    let pisq: f64 = 2.0*PI*PI;
    let mut deltasqk: Vec<f64> = Vec::with_capacity(nhalf);

    for i in 0..nhalf {
        pow_spec[i] *= boxsize.powf(3.0) / (ngrid as f64).powf(6.0);
        pow_spec[i] /= iweights[i] as f64;
        deltasqk[i] = w[i].powf(3.0) * pow_spec[i] / pisq;
    }

    // Return final output
    Ok(Output { w, pow_spec, deltasqk, iweights })
}
