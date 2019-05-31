### crosscorr
[![Build Status](https://travis-ci.com/lewis-weinberger/crosscorr.svg?token=8y51CrY5osH9tmS8LkNZ&branch=master)](https://travis-ci.com/lewis-weinberger/crosscorr)

A power spectrum calculation implemented in Rust, using the [fftw](https://github.com/rust-math/fftw) crate as a 
wrapper to the FFTW3 C library. Takes two 3D fields (stored as row-major 1D
arrays) and uses fourier transforms to calculate the cross-power spectrum of the
fields (or, if the fields are the same, the auto-power spectrum). Credit to 
[Girish Kulkarni](https://github.com/gkulkarni) for the underlying algorithm.

##### Installation
Requires an installation of [Rust](https://www.rust-lang.org/tools/install), and 
the fftw crate requires a C compiler and `make` to compile the underlying source
for the library. With those prerequisites, installation should be as simple as
cloning the repository and using `cargo build --release`. Note compiling the 
`release` version will take longer but provides considerable optimisation which
speeds up the code dramatically (in particular the file IO).

##### Usage
The executable reads in the configuration from a 
[RON](https://github.com/ron-rs/ron) file laid out as follows:

```
(
    grid1_filename: "/path/to/grid1",
    grid2_filename: "/path/to/grid2",
    output_filename: "/path/to/output",
    ngrid: 2048,
    boxsize: 160,
)
```

where the `filename` fields are strings containing the paths to the data, and 
the desired output location; `ngrid` is the number of cells on a side for the
3D data cube; `boxsize` is the physical length of the data cube.

To run the program from the root cargo directory use:

```
cargo run config.ron
```

where `config.ron` is the configuration file explained above. If you have not
compiled the executable then this command will compile before running. The
program will print to stdout/sterr as it progresses.

##### Example
The `examples` directory contains code for preparing the data appropriately, 
`contrast.rs`, as well as for generating mock data, `sinusoid.rs`. There is also
a simple Python plotting script, `plot.py`.
