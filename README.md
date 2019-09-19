# crosscorr
[![Build Status](https://travis-ci.com/lewis-weinberger/crosscorr.svg?token=8y51CrY5osH9tmS8LkNZ&branch=master)](https://travis-ci.com/lewis-weinberger/crosscorr)

> A cosmological power spectrum calculation implemented in Rust.

Calculates the cross-power spectrum of two cosmological fields (or, if the fields are the same, the auto-power spectrum).

Credit to [Girish Kulkarni](https://github.com/gkulkarni) for the underlying algorithm. Uses the [fftw](https://github.com/rust-math/fftw) crate as a wrapper to the FFTW3 C library.

## Installation
Requires an installation of [Rust](https://www.rust-lang.org/tools/install), and the fftw crate requires a C compiler and `make` to compile the underlying source for the [library](http://www.fftw.org/index.html). With those prerequisites, installation should be as simple as cloning the repository and using `cargo build`:

    $ curl https://sh.rustup.rs -sSf | sh                       # install Rust
    $ git clone https://github.com/lewis-weinberger/crosscorr   # clone repository
    $ cd crosscorr                                              # change into source directory
    $ cargo build --release                                     # compile release version

Note compiling the `release` version will take longer but provides considerable optimisation -- this speeds up the code dramatically (in particular the file IO).

## Usage
`crosscorr` is designed to take as input two three-dimensional uniformly gridded fields (stored as row-major 1D arrays). The executable reads in the configuration from a [RON](https://github.com/ron-rs/ron) file laid out as follows:

```
// config.ron
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
3D data cubes; `boxsize` is the physical length of the data cube. Note the data fields should be the same size (i.e. same number of grid cells on a side).

To run the program from the root cargo directory use:

```
$ cargo run config.ron
```

where `config.ron` is the configuration file explained above. If you have not compiled the executable then this command will compile it before running. The program will print to stdout/sterr as it progresses.

#### Window functions
To deconvolve either a nearest-grid-point (NGP) or cloud-in-cell (CIC) mass-assignment for one of the fields, you can pass a feature flag:

    cargo build --release --features=ngp_correction_single

or

    cargo build --release --features=cic_correction_single

(substitute `both` for `single` if you want to correct both fields).

## License

[MIT](./LICENSE)
