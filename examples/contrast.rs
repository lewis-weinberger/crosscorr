use byteorder::{NativeEndian, ReadBytesExt, WriteBytesExt};
use std::env;
use std::fs::File;
use std::io::{BufReader, BufWriter};


/// Calculate the contrast from mean, delta = (x - mean(x))/mean(x), for a
/// gridded field
fn main() -> std::io::Result<()> {
    let args: Vec<String> = env::args().collect();
    let arg1 = &args[1];
    let ngrid: usize = arg1.parse().unwrap();

    let grid_in = &args[2];
    let grid_out = &args[3];
    
    // Read in original grid
    let mut data: Vec<f64> = vec![0.0; ngrid*ngrid*ngrid];
    read_binary(&grid_in, &mut data)?;
   
    // Sanity prints
    for i in 0..5 {
        println!("data[{}] = {}", i, data[i]);
    }

    // Update grid
    update_data(&mut data);
    
    // Sanity prints
    println!("After updating...");
    for i in 0..5 {
        println!("data[{}] = {}", i, data[i]);
    }
    let _check: f64 = data_mean(&data);

    // Save grid to new file
    write_binary(&grid_out, data)?;

    Ok(())
}

fn read_binary(filename: &String, grid: &mut Vec<f64>) -> std::io::Result<()> {
    _read_binary_f32(filename, grid)
    // _read_binary_f64(filename, grid)
}
    
fn _read_binary_f64(filename: &String, grid: &mut Vec<f64>) -> std::io::Result<()> {
    let file = File::open(filename)?;
    println!("opened: {}", filename);
    let mut buf = BufReader::new(file);

    for val in grid.iter_mut() {
        *val = buf.read_f64::<NativeEndian>()?;
    }
    println!("finished reading from file!");

    Ok(())
}

fn _read_binary_f32(filename: &String, grid: &mut Vec<f64>) -> std::io::Result<()> {
    let file = File::open(filename)?;
    println!("opened: {}", filename);
    let mut buf = BufReader::new(file);

    for val in grid.iter_mut() {
        *val = buf.read_f32::<NativeEndian>()? as f64;
    }
    println!("finished reading from file!");

    Ok(())
}
    
fn write_binary(filename: &String, grid: Vec<f64>) -> std::io::Result<()> {
    let file = File::create(filename)?;
    println!("opened: {}", filename);
    let mut buf = BufWriter::new(file);

    for val in &grid {
        buf.write_f64::<NativeEndian>(*val)?;
    }
    println!("finished writing to file!");
    
    Ok(())
}

fn data_mean(data: &Vec<f64>) -> f64 {
    let mut mean_val: f64 = 0.0;

    for val in data {
        mean_val += *val;
    }
    mean_val /= data.len() as f64;

    println!("Mean value found: {}", mean_val);
    mean_val
}

fn update_data(data: &mut Vec<f64>) {
    let mean_val: f64 = data_mean(data);

    for val in data.iter_mut() {
        *val -= mean_val;
        *val /= mean_val;
    }
}
