# R-Math Algorithms

[![crates.io](https://img.shields.io/crates/v/rma.svg)](https://crates.io/crates/rma)
[![docs.rs](https://docs.rs/rma/badge.svg)](https://docs.rs/rma)

A Rust crate for rare, high-performance mathematical algorithms not commonly found in mainstream libraries.

---

## Current Features
- **Chakravala method**: Solves Pell's equation x² - Dy² = 1 for a given non-square integer D.
- **Berlekamp's Algorithm**: Factorization of polynomials over finite fields (GF(p))
    - This implementation supports polynomials over GF(p) for prime p, using u128 for all arithmetic.
    - Uses Karatsuba multiplication for all degree polynomials.
    - On M4Pro12, solves 256 degree polynomial over prime p, the largest 64 bit prime 18_446_744_073_709_551_557 in 
      ```
      berlekamp_roots_dense_deg256_p64b time:   [14.323 ms 14.345 ms 14.366 ms]
      berlekamp_roots_sparse_deg256_p64b time:  [15.582 ms 15.609 ms 15.635 ms]
      ```
    - 128 bit primes are not the part of the benchmark, but values close to the above 64 bit prime will give you similar performance on 256 degree polynomial.

- **Tonelli-Shanks**: Modular square roots (r^2 ≡ n mod p) over prime moduli (constant-time version)
- **Cuthill–McKee & Reverse Cuthill–McKee**: Bandwidth reduction for sparse symmetric matrices (adjacency list, CSR, and CSC formats)
    - **Supports conversion from [`sprs`](https://crates.io/crates/sprs) sparse matrix types**
- **Freivalds' Algorithm**: Fast probabilistic verification of matrix multiplication (scalar, modular, and SIMD-accelerated variants)

- **Well-documented, tested, and benchmarked implementations**
- **SIMD acceleration** for Freivalds' algorithm (nightly Rust required)

### More about the algorithms here

[R-Math Algorithms Overview](R-Math-Algorithms-Overview.md)

---

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
rma = "0.2"
```

Enable SIMD features (nightly Rust required):

```toml
[dependencies]
rma = { version = "0.2", features = ["simd"] }
```
---

## Benchmarks

Benchmarks are provided using [criterion](https://crates.io/crates/criterion) in the `benches/` directory.

**On stable Rust:**
- Run all benchmarks:
  ```sh
  cargo bench
  ```
- Run a specific benchmark (e.g., only Freivalds):
  ```sh
  cargo bench --bench freivalds_bench
  ```

**On nightly Rust (to enable SIMD):**
- Run all benchmarks with SIMD (freivalds_bench is the only one that uses SIMD):
  ```sh
  cargo bench --features simd
  ```
- Run a specific benchmark with SIMD:
  ```sh
  cargo bench --bench freivalds_bench --features simd
  ```

*becnhes naming convention is filename with the `_bench` suffix; you check all in Cargo.toml*

---

## Testing

**On stable Rust:**
- Run all tests:
  ```sh
  cargo test
  ```

**On nightly Rust (to enable SIMD tests):**
- Run all tests with SIMD:
  ```sh
  cargo test --features simd
  ```

---

## Documentation

- [API Documentation (docs.rs)](https://docs.rs/rma)

---

## Minimum Supported Rust Version

- Nightly Rust required for SIMD features (`#![feature(portable_simd)]`)
- Stable Rust for scalar and modular algorithms

---

## Roadmap & Contributions

This crate **will be updated with new algorithms** over time, focusing on rare, advanced, or otherwise unported mathematical algorithms.

---

## License

Licensed under either of
- Apache License, Version 2.0
- MIT license

---

**r-math** aims to be a home for rare, advanced, and high-performance mathematical algorithms in Rust.
