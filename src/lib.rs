use byteorder::{NativeEndian, ReadBytesExt};
use fftw::array::AlignedVec;
use fftw::plan::*;
use fftw::types::*;
use num_complex::Complex;
use ron::de::from_reader;
use serde::Deserialize;
use std::f64::consts::PI;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;

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
/// }
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
    /// Reads the configuration file passed as a command line option  at runtime.
    ///
    /// # Examples
    ///
    /// ```
    /// let config = Config::new(env::args()).unwrap();
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
            Err(_) => return Err("Unable to read configuration from file."),
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
    /// output.save_result(&config).unwrap();
    /// ```
    pub fn save_result(&self, config: &Config) -> Result<(), &'static str> {
        println!("\nSaving results to: {}", &config.output_filename);

        // Open output file
        let mut f = match File::create(&config.output_filename) {
            Ok(file) => file,
            Err(_) => return Err("Unable to open output file!"),
        };
        match writeln!(f, "# w pow_spec deltasqk iweights") {
            Ok(_) => (),
            Err(err) => {
                eprintln!("{}", err);
                return Err("Unable to save output!")
            },
        }

        let nhalf: usize = (config.ngrid / 2) as usize;
        for n in 0..nhalf {
            match writeln!(
                f,
                "{} {} {} {}",
                self.w[n], self.pow_spec[n], self.deltasqk[n], self.iweights[n]
            ) {
                Ok(_) => (),
                Err(err) => {
                    eprintln!("{}", err);
                    return Err("Unable to save output!")
                },
            }
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
/// let grid1 = load_grid(&config, 1).unwrap();
/// ```
pub fn load_grid(config: &Config, num: usize) -> Result<AlignedVec<c64>, &'static str> {
    let filename = match num {
        1 => &config.grid1_filename,
        2 => &config.grid2_filename,
        _ => return Err("Need to load either grid 1 or 2!"),
    };
    println!("\nOpening grid from file: {}", filename);
    let ngrid: usize = config.ngrid as usize;

    // Allocate AlignedVec array to hold grid
    let ngrid3 = ngrid * ngrid * ngrid;
    let mut grid = AlignedVec::new(ngrid3);

    // Open binary file
    let f = match File::open(filename) {
        Ok(file) => file,
        Err(_) => return Err("Unable to open grid file!"),
    };
    let mut buf_reader = BufReader::new(f);

    // Read in array from binary file
    for elem in grid.iter_mut() {
        let cell = match buf_reader.read_f32::<NativeEndian>() {
            Ok(val) => val,
            Err(_) => return Err("Problem reading values from file!"),
        };
        *elem = c64::new(f64::from(cell), 0.0);
    }
    println!("Successfully read {} cells!", ngrid3);
    println!("Sanity print:");
    grid[0..5].iter()
        .enumerate()
        .for_each(|(i, elem)| {
        println!("grid1[{}] = {:.3e} + {:.3e}i", i, elem.re, elem.im);
    });

    Ok(grid)
}

/// Performs FFT on grids
///
/// # Examples
///
/// ```
/// let output: Output = correlate(&config, grid1, grid2).unwrap();
/// ```
pub fn perform_fft(
    config: &Config,
    grid1: AlignedVec<c64>,
    grid2: AlignedVec<c64>,
) -> Result<(AlignedVec<c64>, AlignedVec<c64>), &'static str> {
    println!("\nPerforming FFTs...");
    let ngrid: usize = config.ngrid as usize;

    // Create FFTW plan
    let shape = [ngrid, ngrid, ngrid];
    let mut plan: C2CPlan64 = match C2CPlan::aligned(&shape[..], Sign::Forward, Flag::Estimate) {
        Ok(p) => p,
        Err(_) => return Err("Unable to create FFTW plan."),
    };
    println!("Plan created!");

    // Perform FFT on grids
    let ngrid3 = ngrid * ngrid * ngrid;

    let out1 = fft_from_plan(ngrid3, grid1, &mut plan)?;
    println!("First grid FFT complete!");

    let out2 = fft_from_plan(ngrid3, grid2, &mut plan)?;
    println!("Second grid FFT complete!");

    // Sanity prints
    println!("FFTs performed... Sanity check:");
    for n in 0..10 {
        println!("out1[{}] = {:.3e} + {:.3e}i", n, out1[n].re, out1[n].im);
        println!("out2[{}] = {:.3e} + {:.3e}i", n, out2[n].re, out2[n].im);
    }
    Ok((out1, out2))
}

/// Use FFTW3 plan to perform FFT
fn fft_from_plan(
    ngrid3: usize,
    mut grid: AlignedVec<c64>,
    plan: &mut C2CPlan64,
) -> Result<AlignedVec<c64>, &'static str> {
    let mut out = AlignedVec::new(ngrid3);
    match plan.c2c(&mut grid, &mut out) {
        Ok(_) => (),
        Err(_) => return Err("Failed to FFT grid."),
    };
    Ok(out)
}

/// Calculates the cross power spectrum of the given 3D grids (note if the same
/// grid is given twice then this is the auto power spectrum).
///
/// # Examples
///
/// ```
/// let output: Output = correlate(&config, grid1, grid2).unwrap();
/// ```
pub fn correlate(
    config: &Config,
    out1: AlignedVec<c64>,
    out2: AlignedVec<c64>,
) -> Result<Output, &'static str> {
    println!("\nCalculating power spectrum...");

    if cfg!(feature = "ngp_correction_single") {
        println!("Correcting for NGP mass assignment of one field!");
    } else if cfg!(feature = "cic_correction_single") {
        println!("Correcting for CIC mass assignment of one field!");
    } else if cfg!(feature = "ngp_correction_both") {
        println!("Correcting for NGP mass assignment of both fields!");
    } else if cfg!(feature = "cic_correction_both") {
        println!("Correcting for CIC mass assignment of both fields!");
    }

    let ngrid: usize = config.ngrid as usize;
    let boxsize: f64 = f64::from(config.boxsize);

    // Calculate power spectrum
    let kf: f64 = 2.0 * PI / boxsize;
    let coeff: f64 = (boxsize / (2.0 * PI)).powf(2.0);
    let nhalf: usize = ngrid / 2;

    #[cfg(any(
        feature = "ngp_correction_single",
        feature = "ngp_correction_both",
        feature = "cic_correction_single",
        feature = "cic_correction_both"
    ))]
    let kny: f64 = PI * config.ngrid as f64 / boxsize;

    let mut w: Vec<f64> = Vec::with_capacity(ngrid);
    for i in 0..=nhalf {
        w.push(kf * (i as f64));
    }
    for i in (nhalf + 1)..ngrid {
        w.push(kf * ((i as isize - ngrid as isize) as f64));
    }

    let mut pow_spec: Vec<f64> = vec![0.0; ngrid];
    let mut iweights: Vec<i64> = vec![0; ngrid];

    for i in 0..ngrid {
        let iper = if i >= nhalf { ngrid - i } else { i };
        for j in 0..ngrid {
            let jper = if j >= nhalf { ngrid - j } else { j };
            for k in 0..ngrid {
                let kper = if k >= nhalf { ngrid - k } else { k };
                let r: f64 = (iper * iper + jper * jper + kper * kper) as f64;
                let m: usize = (0.5 + r.sqrt()) as usize;
                iweights[m] += 1;

                let g = w[i] * w[i] + w[j] * w[j] + w[k] * w[k];
                if g != 0.0 {
                    let scale: usize = (0.5 + (g * coeff).sqrt()) as usize;
                    let index: usize = k + ngrid * (j + ngrid * i);
                    let mut contrib: Complex<f64> =
                        out1[index] * out2[index].conj() + out1[index].conj() * out2[index];

                    #[cfg(feature = "ngp_correction_single")]
                    {
                        // Correct for Nearest-Grid-Point mass assignment
                        let wngp = sinc(PI * w[i] as f64 / (2.0 * kny))
                            * sinc(PI * w[j] as f64 / (2.0 * kny))
                            * sinc(PI * w[k] as f64 / (2.0 * kny));
                        contrib.re /= wngp;
                    }

                    #[cfg(feature = "cic_correction_single")]
                    {
                        // Correct for Cloud-in-Cell mass assignment
                        let wcic = (sinc(PI * w[i] as f64 / (2.0 * kny))
                            * sinc(PI * w[j] as f64 / (2.0 * kny))
                            * sinc(PI * w[k] as f64 / (2.0 * kny)))
                        .powi(2);
                        contrib.re /= wcic;
                    }

                    #[cfg(feature = "ngp_correction_both")]
                    {
                        // Correct for Nearest-Grid-Point mass assignment
                        let wngp = sinc(PI * w[i] as f64 / (2.0 * kny))
                            * sinc(PI * w[j] as f64 / (2.0 * kny))
                            * sinc(PI * w[k] as f64 / (2.0 * kny));
                        contrib.re /= wngp * wngp;
                    }

                    #[cfg(feature = "cic_correction_both")]
                    {
                        // Correct for Cloud-in-Cell mass assignment
                        let wcic = (sinc(PI * w[i] as f64 / (2.0 * kny))
                            * sinc(PI * w[j] as f64 / (2.0 * kny))
                            * sinc(PI * w[k] as f64 / (2.0 * kny)))
                        .powi(2);
                        contrib.re /= wcic * wcic;
                    }

                    pow_spec[scale] += contrib.re / 2.0;
                }
            }
        }
    }
    println!("Power spectrum calculated. Normalising...");

    // Normalise power spectrum
    let pisq: f64 = 2.0 * PI * PI;
    let mut deltasqk: Vec<f64> = Vec::with_capacity(nhalf);

    for i in 0..nhalf {
        pow_spec[i] *= boxsize.powi(3) / (ngrid as f64).powi(6);
        pow_spec[i] /= iweights[i] as f64;
        deltasqk.push(w[i].powf(3.0) * pow_spec[i] / pisq);
    }

    // Return final output
    Ok(Output {
        w,
        pow_spec,
        deltasqk,
        iweights,
    })
}

#[cfg(any(
    feature = "ngp_correction_single",
    feature = "ngp_correction_both",
    feature = "cic_correction_single",
    feature = "cic_correction_both"
))]
fn sinc(theta: f64) -> f64 {
    if theta < 1e-20 {
        1.0
    } else {
        (theta.sin() / theta)
    }
}
