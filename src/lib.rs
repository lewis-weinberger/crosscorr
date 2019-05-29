extern crate fftw;
extern crate ron;
extern crate serde;
extern crate byteorder;

use std::fs::File;
use std::io::BufReader;
use fftw::array::AlignedVec;
use fftw::types::c64;
use byteorder::{ReadBytesExt, NativeEndian};
use ron::de::from_reader;
use serde::Deserialize;

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
///     boxsize: 160,
///     }
/// ```
#[derive(Debug, Deserialize)]
pub struct Config {
    pub grid1_filename: String,
    pub grid2_filename: String,
    pub output_filename: String,
    pub ngrid: u32,
    pub boxsize: u32,
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

/// Loads a grid stored at `filename` (in a custom binary format) into an 
/// `fftw::array::AlignedVec` object. This custom format stores the 3D grid as
/// a 1D array of values. The data should be stored as deviations from the mean,
/// i.e. delta = (x - mean(x)) / mean(x).
///
/// # Examples
///
/// ```
/// ```
pub fn load_grid(filename: &String, ngrid: u32) -> Result<AlignedVec<c64>, &'static str> {
    println!("\nOpening grid from file: {}", filename);

    // Allocate AlignedVec array to hold grid
    let ngrid3 = (ngrid as usize)*(ngrid as usize)*(ngrid as usize);
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
pub fn correlate() {}

/// Saves the power spectrum to a formatted txt file.
///
/// # Examples
///
/// ```
/// ```
pub fn save_result() {}
