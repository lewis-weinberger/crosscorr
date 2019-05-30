use std::env;
use std::f64::consts::PI;
use std::io::prelude::*;
use std::io::BufWriter;
use std::fs::File;
use byteorder::{NativeEndian, WriteBytesExt};

/// Generate a mock 3D data sinusoidal field and the corresponding config
/// file, for use in testing.
fn main() -> std::io::Result<()> {
    let args: Vec<String> = env::args().collect();
    let base = &args[1];
    println!("Working in base directory: {}", base);
    
    let ngrid = 256;
    let boxsize = 160;

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
    let mut data: Vec<f64> = vec![0.0; ngrid*ngrid*ngrid];
    for i in 0..ngrid {
        for j in 0..ngrid {
            for k in 0..ngrid {
                let theta = ((i + j + k) / ngrid) as f64;
                data[k + ngrid*(j + ngrid*i)] = (theta * 4.0 * PI).sin();
            }
        }
    }
    
    // Write test data to file
    let grid_filename = format!("{}/test_grid.dat", &base);
    let fgrid = File::create(&grid_filename)?;
    let mut buf = BufWriter::new(fgrid);

    let ngrid3 = ngrid*ngrid*ngrid;
    for i in 0..ngrid3 {
        buf.write_f64::<NativeEndian>(data[i])?;
    }
    println!("Created mock data file: {}", &grid_filename);

    Ok(())
}
