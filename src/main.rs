use std::env;
use std::process;

use crosscorr;
use crosscorr::{correlate, load_grid, perform_fft, Config, Output};

fn main() {
    // Load configuration
    let config = Config::new(env::args()).unwrap_or_else(|err| {
        eprintln!("Error: {}", err);
        process::exit(1);
    });

    // Load first grid
    let grid1 = load_grid(&config, 1).unwrap_or_else(|err| {
        eprintln!("Error: {}", err);
        process::exit(1);
    });

    // Load second grid
    let grid2 = load_grid(&config, 2).unwrap_or_else(|err| {
        eprintln!("Error: {}", err);
        process::exit(1);
    });

    // Perform FFT
    let (out1, out2) = perform_fft(&config, grid1, grid2).unwrap_or_else(|err| {
        eprintln!("Error: {}", err);
        process::exit(1);
    });

    // Calculate power spectrum
    let output: Output = correlate(&config, out1, out2).unwrap_or_else(|err| {
        eprintln!("Error: {}", err);
        process::exit(1);
    });

    // Save power spectrum to output file
    output.save_result(&config).unwrap_or_else(|err| {
        eprintln!("Error: {}", err);
        process::exit(1);
    });
}
