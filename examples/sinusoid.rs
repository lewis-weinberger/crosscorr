use byteorder::{NativeEndian, WriteBytesExt};
use std::env;
use std::f64::consts::PI;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufWriter;

/// Generate a mock 3D data sinusoidal field and the corresponding config
/// file, for use in testing.
fn main() -> std::io::Result<()> {
    let args: Vec<String> = env::args().collect();
    let base = &args[1];
    println!("Working in base directory: {}", base);

    let ngrid = 512;
    let boxsize = 50;
    let wavelength = 0.1 * boxsize;

    // First generate the config file in RON format
    let config_filename = format!("{}/config.ron", &base);
    let mut f = File::create(&config_filename).unwrap();
    write!(f, "(\n")?;
    write!(f, "    grid1_filename: \"{}/test_grid.dat\",\n", base)?;
    write!(f, "    grid2_filename: \"{}/test_grid.dat\",\n", base)?;
    write!(f, "    output_filename: \"{}/output.dat\",\n", base)?;
    write!(f, "    ngrid: {},\n", ngrid)?;
    write!(f, "    boxsize: {},\n", boxsize)?;
    write!(f, ")\n")?;
    println!("Created config file: {}", &config_filename);

    // Create test grid
    let mut data: Vec<f64> = vec![0.0; ngrid * ngrid * ngrid];
    for i in 0..ngrid {
        for j in 0..ngrid {
            for k in 0..ngrid {
                let theta = (i + j + k) as f64 / (ngrid as f64);
                theta *= (boxsize / wavelength);
                data[k + ngrid * (j + ngrid * i)] = (theta * PI / 2.0).sin();
            }
        }
    }
    println!("Sanity print:");
    for i in 0..5 {
        println!("data[{}] = {}", i, data[i]);
    }

    // Write test data to file
    let grid_filename = format!("{}/test_grid.dat", &base);
    let fgrid = File::create(&grid_filename)?;
    let mut buf = BufWriter::new(fgrid);

    for val in data.iter() {
        buf.write_f64::<NativeEndian>(*val)?;
    }
    println!("Created mock data file: {}", &grid_filename);

    Ok(())
}
