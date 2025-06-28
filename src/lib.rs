#![doc = r#"
# R-Math Algorithms

A Rust crate for rare, high-performance mathematical algorithms not commonly found in mainstream libraries aimed for practical use.

## Features
- Berlekamp's Algorithm: Factorization of polynomials over finite fields (GF(p))
- Tonelli-Shanks: Modular square roots (r^2 ≡ n mod p) over prime moduli (constant-time version)
- Cuthill-McKee & Reverse Cuthill-McKee: Bandwidth reduction for sparse symmetric matrices (adjacency list, CSR, and CSC formats)
    - Supports conversion from `sprs` sparse matrix types
- Freivalds' Algorithm: Fast probabilistic verification of matrix multiplication (scalar, modular, and SIMD-accelerated variants)
- Chakravala: Solves Pell's equation x² - d*y² = 1 for a given non-square integer d.

- Well-documented, tested, and benchmarked implementations
- SIMD acceleration for Freivalds' algorithm (nightly Rust required, `simd` feature)
- `no_std` compatible (except for benchmarks and RNG)

## Usage

Add to your `Cargo.toml`:
```toml
[dependencies]
rma = "0.3"
```

Enable SIMD features (nightly Rust required):
```toml
[dependencies]
rma = { version = "0.3", features = ["simd"] }
```

## SIMD Feature
- The `simd` feature enables SIMD-accelerated Freivalds' algorithm. Requires nightly Rust and the `portable_simd` feature.
- On stable Rust, only scalar and modular versions are available.

## Examples

See the documentation for each algorithm.

## License
Licensed under either of
- Apache License, Version 2.0
- MIT license

"#]
#![cfg_attr(feature = "simd", feature(portable_simd))]

pub mod freivalds;

pub use freivalds::freivalds_verify_scalar;
pub use freivalds::freivalds_verify_scalar_mod;
#[cfg(feature = "simd")]
pub use freivalds::freivalds_verify_simd;

pub mod cuthill_mckee;
pub use cuthill_mckee::{
    CscMatrix, CsrMatrix, cuthill_mckee, cuthill_mckee_csc, cuthill_mckee_csr,
    reverse_cuthill_mckee, reverse_cuthill_mckee_csc, reverse_cuthill_mckee_csr,
};

pub mod tonelli_shanks;
pub use tonelli_shanks::{TonelliShanksError, tonelli_shanks_ct};

pub mod berlekamp;

pub mod chakravala;
pub use chakravala::chakravala;
