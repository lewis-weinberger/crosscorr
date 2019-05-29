use std::env;
use std::process;

use crosscorr;
use crosscorr::{Config, load_grid, correlate, save_result};

fn main() {

    // Load configuration
    let config = Config::new(env::args()).unwrap_or_else(|err| {
        eprintln!("Error: {}", err);
        process::exit(1);
    });

    // Load first grid
    let mut grid1 = load_grid(&config.grid1_filename, config.ngrid).unwrap_or_else(|err| {
        eprintln!("Error: {}", err);
        process::exit(1);
    });

    // Load second grid
    let mut grid2 = load_grid(&config.grid2_filename, config.ngrid).unwrap_or_else(|err| {
        eprintln!("Error: {}", err);
        process::exit(1);
    });

    // Calculate power spectrum
    correlate();
    
    // Save power spectrum to output file
    save_result();
}
