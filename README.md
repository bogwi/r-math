# R-Math Algorithms

[![crates.io](https://img.shields.io/crates/v/rma.svg)](https://crates.io/crates/rma)
[![docs.rs](https://docs.rs/rma/badge.svg)](https://docs.rs/rma)

A Rust crate for rare, high-performance mathematical algorithms not commonly found in mainstream libraries.

---

## Current Features

- **Tonelli-Shanks**: Modular square roots (r^2 ≡ n mod p) over prime moduli (constant-time version)
- **Cuthill–McKee & Reverse Cuthill–McKee**: Bandwidth reduction for sparse symmetric matrices (adjacency list, CSR, and CSC formats)
    - **Supports conversion from [`sprs`](https://crates.io/crates/sprs) sparse matrix types**
- **Freivalds' Algorithm**: Fast probabilistic verification of matrix multiplication (scalar, modular, and SIMD-accelerated variants)

- **Well-documented, tested, and benchmarked implementations**
- **SIMD acceleration** for Freivalds' algorithm (nightly Rust required)

### More about the algorithms here

[R-Math Algorithms Overview](R-Math-Algorithms-Overwiew.md)

---

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
rma = "0.1"
```

Enable SIMD features (nightly Rust required):

```toml
[dependencies]
rma = { version = "0.1", features = ["simd"] }
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
- Run all benchmarks with SIMD:
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
