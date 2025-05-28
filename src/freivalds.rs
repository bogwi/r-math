//! Freivalds' algorithm for fast probabilistic matrix multiplication verification.
//!
//! This implementation provides three variants:
//! - `freivalds_verify_scalar`: Basic scalar implementation
//! - `freivalds_verify_scalar_mod`: Modular arithmetic version
//! - `freivalds_verify_simd`: SIMD-accelerated version (requires nightly Rust and `simd` feature)
//!
//! # Algorithm
//! Given matrices A, B, and C, verifies if A × B = C by:
//! 1. Generating a random vector r with entries 0 or 1
//! 2. Computing A × (B × r) and C × r
//! 3. Comparing the results
//!
//! # Error Probability
//! - If A × B = C, always returns true (no false negatives)
//! - If A × B ≠ C, probability of false positive is at most 1/2 per run
//!
//! # Time Complexity
//! O(n²) vs O(n³) for full matrix multiplication
//!
//! # Example
//! ```
//! use rma::freivalds::freivalds_verify_scalar;
//! let a = vec![vec![1, 2], vec![3, 4]];
//! let b = vec![vec![5, 6], vec![7, 8]];
//! let c = vec![vec![19, 22], vec![43, 50]];
//! assert!(freivalds_verify_scalar(&a, &b, &c, 10));
//! ```
//!
//! # References
//! - [Wikipedia](https://en.wikipedia.org/wiki/Freivalds'_algorithm)
//! - [OpenGenus IQ](https://iq.opengenus.org/freivalds-algorithm/)

#[cfg(feature = "simd")]
use core::simd::Simd;
use num_traits::PrimInt;
use rand::Rng;

// Generic scalar dot product for any integer type.
pub fn scalar_dot<T: PrimInt + Copy>(a: &[T], b: &[T]) -> T {
    a.iter()
        .zip(b.iter())
        .fold(T::zero(), |acc, (x, y)| acc + (*x * *y))
}

/// Scalar dot product with modular arithmetic for any integer type.
fn scalar_dot_mod<T: PrimInt + Copy + RemEuclid>(a: &[T], b: &[T], modulus: T) -> T {
    a.iter().zip(b.iter()).fold(T::zero(), |acc, (x, y)| {
        (acc + (*x * *y)).rem_euclid(modulus)
    })
}

#[cfg(feature = "simd")]
#[inline]
pub fn simd_dot<T, const LANES: usize>(a: &[T], b: &[T]) -> T
where
    T: core::simd::SimdElement
        + core::ops::Mul<Output = T>
        + core::ops::Add<Output = T>
        + Copy
        + Default,
    core::simd::LaneCount<LANES>: core::simd::SupportedLaneCount,
    core::simd::Simd<T, LANES>: core::ops::Mul<core::simd::Simd<T, LANES>, Output = core::simd::Simd<T, LANES>>
        + core::ops::Add<core::simd::Simd<T, LANES>, Output = core::simd::Simd<T, LANES>>,
{
    let mut sum = Simd::<T, LANES>::splat(T::default());
    let chunks = a.len() / LANES;

    (0..chunks).for_each(|i| {
        let ai = Simd::<T, LANES>::from_slice(&a[i * LANES..i * LANES + LANES]);
        let bi = Simd::<T, LANES>::from_slice(&b[i * LANES..i * LANES + LANES]);
        sum = sum + (ai * bi);
    });

    let mut total = sum.as_array().iter().fold(T::default(), |acc, &x| acc + x);
    total = total
        + a[chunks * LANES..]
            .iter()
            .zip(b[chunks * LANES..].iter())
            .fold(T::default(), |acc, (x, y)| acc + (*x * *y));
    total
}

// Helper trait for rem_euclid for both signed and unsigned
pub trait RemEuclid: PrimInt {
    fn rem_euclid(self, modulus: Self) -> Self;
}

impl RemEuclid for i64 {
    fn rem_euclid(self, modulus: Self) -> Self {
        ((self % modulus) + modulus) % modulus
    }
}

impl RemEuclid for u64 {
    fn rem_euclid(self, modulus: Self) -> Self {
        self % modulus
    }
}

impl RemEuclid for i32 {
    fn rem_euclid(self, modulus: Self) -> Self {
        ((self % modulus) + modulus) % modulus
    }
}

impl RemEuclid for u32 {
    fn rem_euclid(self, modulus: Self) -> Self {
        self % modulus
    }
}

/// Verifies if AB = C for n x n matrices using Freivalds' algorithm (SIMD version).
///
/// # Arguments
/// * `LANES` - Number of SIMD lanes to use (must be a power of 2 and supported by the platform)
/// * `a` - Left matrix (n x n)
/// * `b` - Right matrix (n x n)
/// * `c` - Product matrix (n x n)
/// * `k` - Number of iterations (higher = lower error probability)
///
/// # Returns
/// * `true` if the verification passes, `false` otherwise.
///
/// # Example
/// ```
/// use rma::freivalds_verify_simd;
/// let a = vec![vec![1, 2], vec![3, 4]];
/// let b = vec![vec![2, 0], vec![1, 2]];
/// let c = vec![vec![4, 4], vec![10, 8]]; // c = a * b
/// assert!(freivalds_verify_simd::<u64, 4>(&a, &b, &c, 10));
/// ```
///
/// # Note
/// This version uses explicit SIMD for the dot product, but currently only supports unsigned integer types (`u64`).
/// For signed or modular arithmetic, use `freivalds_verify_scalar_mod`.
/// The `LANES` parameter must be a power of 2 and supported by the platform (e.g., 4, 8, 16, 32, 64).
/// Using simd might be slower than the scalar version for any matrix size.
#[cfg(feature = "simd")]
#[inline]
pub fn freivalds_verify_simd<T, const LANES: usize>(
    a: &[Vec<T>],
    b: &[Vec<T>],
    c: &[Vec<T>],
    k: usize,
) -> bool
where
    T: PrimInt
        + core::simd::SimdElement
        + Copy
        + rand::distributions::uniform::SampleUniform
        + Default,
    core::simd::LaneCount<LANES>: core::simd::SupportedLaneCount,
    core::simd::Simd<T, LANES>: core::ops::Mul<core::simd::Simd<T, LANES>, Output = core::simd::Simd<T, LANES>>
        + core::ops::Add<core::simd::Simd<T, LANES>, Output = core::simd::Simd<T, LANES>>,
{
    let n = a.len();
    if n == 0
        || b.len() != n
        || c.len() != n
        || a[0].len() != n
        || b[0].len() != n
        || c[0].len() != n
    {
        return false;
    }
    let mut rng = rand::thread_rng();

    for _ in 0..k {
        let r: Vec<T> = (0..n)
            .map(|_| T::from(rng.gen_range(0u8..=1u8)).unwrap())
            .collect();

        let br: Vec<T> = b.iter().map(|row| simd_dot::<T, LANES>(row, &r)).collect();
        let abr: Vec<T> = a.iter().map(|row| simd_dot::<T, LANES>(row, &br)).collect();
        let cr: Vec<T> = c.iter().map(|row| simd_dot::<T, LANES>(row, &r)).collect();

        if abr != cr {
            return false;
        }
    }
    true
}

/// Generalized scalar Freivalds' algorithm supporting unsigned integer types and arbitrary modulus.
///
/// Note: Only unsigned integer types are supported.
///
/// # Arguments
/// * `a` - Left matrix (n x n)
/// * `b` - Right matrix (n x n)
/// * `c` - Product matrix (n x n)
/// * `k` - Number of iterations (higher = lower error probability)
/// * `modulus` - Modulus for all arithmetic (must be > 0)
///
/// # Returns
/// * `true` if the verification passes, `false` otherwise.
///
/// # Example
/// ```
/// use rma::freivalds_verify_scalar_mod;
/// let a = vec![vec![1i64, 2], vec![3, 4]];
/// let b = vec![vec![2, 0], vec![1, 2]];
/// let c = vec![vec![4, 4], vec![10, 8]]; // c = a * b
/// assert!(freivalds_verify_scalar_mod(&a, &b, &c, 10, 1_000_000_007));
/// ```
///
/// # Note
/// This version supports both signed and unsigned integer types, and arbitrary modulus.
/// The SIMD version (`freivalds_verify_simd`) currently supports only unsigned integers.
pub fn freivalds_verify_scalar_mod<T>(
    a: &[Vec<T>],
    b: &[Vec<T>],
    c: &[Vec<T>],
    k: usize,
    modulus: T,
) -> bool
where
    T: PrimInt + Copy + rand::distributions::uniform::SampleUniform + RemEuclid,
{
    let n = a.len();
    if n == 0
        || b.len() != n
        || c.len() != n
        || a[0].len() != n
        || b[0].len() != n
        || c[0].len() != n
        || modulus <= T::zero()
    {
        return false;
    }
    let mut rng = rand::thread_rng();
    for _ in 0..k {
        // For unsigned: {0,1}, for signed: {0,1} (could be extended to {-1,0,1} if desired)
        let r: Vec<T> = (0..n)
            .map(|_| T::from(rng.gen_range(0u8..=1u8)).unwrap())
            .collect();
        let br: Vec<T> = b
            .iter()
            .map(|row| scalar_dot_mod(row, &r, modulus))
            .collect();
        let abr: Vec<T> = a
            .iter()
            .map(|row| scalar_dot_mod(row, &br, modulus))
            .collect();
        let cr: Vec<T> = c
            .iter()
            .map(|row| scalar_dot_mod(row, &r, modulus))
            .collect();
        if abr != cr {
            return false;
        }
    }
    true
}

/// Verifies if AB = C for n x n matrices using Freivalds' algorithm (scalar, non-modular, signed/unsigned).
///
/// # Arguments
/// * `a` - Left matrix (n x n)
/// * `b` - Right matrix (n x n)
/// * `c` - Product matrix (n x n)
/// * `k` - Number of iterations (higher = lower error probability)
///
/// # Returns
/// * `true` if the verification passes, `false` otherwise.
///
/// # Example
/// ```
/// use rma::freivalds_verify_scalar;
/// let a = vec![vec![1, 2], vec![3, 4]];
/// let b = vec![vec![2, 0], vec![1, 2]];
/// let c = vec![vec![4, 4], vec![10, 8]]; // c = a * b
/// assert!(freivalds_verify_scalar(&a, &b, &c, 10));
/// ```
pub fn freivalds_verify_scalar<T>(a: &[Vec<T>], b: &[Vec<T>], c: &[Vec<T>], k: usize) -> bool
where
    T: PrimInt + Copy + rand::distributions::uniform::SampleUniform,
{
    let n = a.len();
    if n == 0
        || b.len() != n
        || c.len() != n
        || a[0].len() != n
        || b[0].len() != n
        || c[0].len() != n
    {
        return false;
    }
    let mut rng = rand::thread_rng();
    for _ in 0..k {
        // For unsigned: {0,1}, for signed: {0,1}
        let r: Vec<T> = (0..n)
            .map(|_| T::from(rng.gen_range(0u8..=1u8)).unwrap())
            .collect();
        let br: Vec<T> = b.iter().map(|row| scalar_dot(row, &r)).collect();
        let abr: Vec<T> = a.iter().map(|row| scalar_dot(row, &br)).collect();
        let cr: Vec<T> = c.iter().map(|row| scalar_dot(row, &r)).collect();
        if abr != cr {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_freivalds_correct() {
        let a = vec![vec![1, 2], vec![3, 4]];
        let b = vec![vec![2, 0], vec![1, 2]];
        let c = vec![vec![4, 4], vec![10, 8]]; // c = a * b
        #[cfg(feature = "simd")]
        assert!(freivalds_verify_simd::<u64, 4>(&a, &b, &c, 10));
        assert!(freivalds_verify_scalar_mod(&a, &b, &c, 10, 1_000_000_007));
    }

    #[test]
    fn test_freivalds_incorrect() {
        let a = vec![vec![1, 2], vec![3, 4]];
        let b = vec![vec![2, 0], vec![1, 2]];
        let c = vec![vec![4, 4], vec![10, 9]]; // c is not a * b
        #[cfg(feature = "simd")]
        assert!(!freivalds_verify_simd::<u64, 4>(&a, &b, &c, 10));
        assert!(!freivalds_verify_scalar_mod(&a, &b, &c, 10, 1_000_000_007));
    }

    #[test]
    fn test_freivalds_zero_matrix() {
        let a: Vec<Vec<u64>> = vec![];
        let b: Vec<Vec<u64>> = vec![];
        let c: Vec<Vec<u64>> = vec![];
        #[cfg(feature = "simd")]
        assert!(!freivalds_verify_simd::<u64, 4>(&a, &b, &c, 10));
        assert!(!freivalds_verify_scalar_mod(&a, &b, &c, 10, 1_000_000_007));
    }

    #[test]
    fn test_freivalds_mismatched_sizes() {
        let a = vec![vec![1, 2], vec![3, 4]];
        let b = vec![vec![2, 0]]; // wrong size
        let c = vec![vec![4, 4], vec![10, 8]];
        #[cfg(feature = "simd")]
        assert!(!freivalds_verify_simd::<u64, 4>(&a, &b, &c, 10));
        assert!(!freivalds_verify_scalar_mod(&a, &b, &c, 10, 1_000_000_007));
    }

    #[test]
    fn test_freivalds_larger_matrix() {
        let n = 8;
        let mut a = vec![vec![0u64; n]; n];
        let mut b = vec![vec![0u64; n]; n];
        for i in 0..n {
            for j in 0..n {
                a[i][j] = (i * j + 1) as u64;
                b[i][j] = ((i + j) % 5 + 1) as u64;
            }
        }
        // Compute c = a * b
        let mut c = vec![vec![0u64; n]; n];
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    c[i][j] = c[i][j].wrapping_add(a[i][k].wrapping_mul(b[k][j]));
                }
            }
        }
        #[cfg(feature = "simd")]
        assert!(freivalds_verify_simd::<u64, 4>(&a, &b, &c, 10));
        assert!(freivalds_verify_scalar_mod(&a, &b, &c, 10, 1_000_000_007));
    }

    #[test]
    fn test_freivalds_all_zeros() {
        let a = vec![vec![0, 0], vec![0, 0]];
        let b = vec![vec![0, 0], vec![0, 0]];
        let c = vec![vec![0, 0], vec![0, 0]];
        #[cfg(feature = "simd")]
        assert!(freivalds_verify_simd::<u64, 4>(&a, &b, &c, 10));
        assert!(freivalds_verify_scalar_mod(&a, &b, &c, 10, 1_000_000_007));
    }

    #[test]
    fn test_freivalds_all_ones() {
        let a = vec![vec![1, 1], vec![1, 1]];
        let b = vec![vec![1, 1], vec![1, 1]];
        let c = vec![vec![2, 2], vec![2, 2]]; // c = a * b
        #[cfg(feature = "simd")]
        assert!(freivalds_verify_simd::<u64, 4>(&a, &b, &c, 10));
        assert!(freivalds_verify_scalar_mod(&a, &b, &c, 10, 1_000_000_007));
    }

    #[test]
    fn test_freivalds_single_element() {
        let a = vec![vec![7]];
        let b = vec![vec![3]];
        let c = vec![vec![21]];
        #[cfg(feature = "simd")]
        assert!(freivalds_verify_simd::<u64, 4>(&a, &b, &c, 10));
        assert!(freivalds_verify_scalar_mod(&a, &b, &c, 10, 1_000_000_007));
    }

    #[test]
    fn test_freivalds_large_values() {
        let max = u64::MAX;
        let a = vec![vec![max, max], vec![max, max]];
        let b = vec![vec![1, 0], vec![0, 1]]; // identity
        let c = vec![vec![max, max], vec![max, max]];
        #[cfg(feature = "simd")]
        assert!(freivalds_verify_simd::<u64, 4>(&a, &b, &c, 10));
        assert!(freivalds_verify_scalar_mod(&a, &b, &c, 10, 1_000_000_007));
    }

    #[test]
    fn test_freivalds_k_zero() {
        let a = vec![vec![1, 2], vec![3, 4]];
        let b = vec![vec![2, 0], vec![1, 2]];
        let c = vec![vec![4, 4], vec![10, 8]];
        // k=0: should always return true (no checks performed)
        #[cfg(feature = "simd")]
        assert!(freivalds_verify_simd::<u64, 4>(&a, &b, &c, 0));
        assert!(freivalds_verify_scalar_mod(&a, &b, &c, 0, 1_000_000_007));
    }

    #[test]
    fn test_freivalds_scalar_correct() {
        let a = vec![vec![1, 2], vec![3, 4]];
        let b = vec![vec![2, 0], vec![1, 2]];
        let c = vec![vec![4, 4], vec![10, 8]]; // c = a * b
        assert!(freivalds_verify_scalar(&a, &b, &c, 10));
    }

    #[test]
    fn test_freivalds_scalar_incorrect() {
        let a = vec![vec![1, 2], vec![3, 4]];
        let b = vec![vec![2, 0], vec![1, 2]];
        let c = vec![vec![4, 4], vec![10, 9]]; // incorrect
        assert!(!freivalds_verify_scalar(&a, &b, &c, 10));
    }

    #[test]
    fn test_freivalds_scalar_zero_matrix() {
        let a: Vec<Vec<i64>> = vec![];
        let b: Vec<Vec<i64>> = vec![];
        let c: Vec<Vec<i64>> = vec![];
        assert!(!freivalds_verify_scalar(&a, &b, &c, 10));
    }

    #[test]
    fn test_freivalds_scalar_mismatched_sizes() {
        let a = vec![vec![1, 2], vec![3, 4]];
        let b = vec![vec![2, 0]]; // wrong size
        let c = vec![vec![4, 4], vec![10, 8]];
        assert!(!freivalds_verify_scalar(&a, &b, &c, 10));
    }

    #[test]
    fn test_freivalds_scalar_negatives() {
        let a = vec![vec![-1, 2], vec![3, -4]];
        let b = vec![vec![2, 0], vec![1, -2]];
        let mut c = vec![vec![0i64; 2]; 2];
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    c[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        assert!(freivalds_verify_scalar(&a, &b, &c, 10));
    }

    #[test]
    fn test_freivalds_scalar_large_values() {
        let max = i64::MAX;
        let a = vec![vec![max, max], vec![max, max]];
        let b = vec![vec![1, 0], vec![0, 1]]; // identity
        let c = vec![vec![max, max], vec![max, max]];
        assert!(freivalds_verify_scalar(&a, &b, &c, 10));
    }

    #[test]
    fn test_freivalds_scalar_k_zero() {
        let a = vec![vec![1, 2], vec![3, 4]];
        let b = vec![vec![2, 0], vec![1, 2]];
        let c = vec![vec![4, 4], vec![10, 8]];
        // k=0: should always return true (no checks performed)
        assert!(freivalds_verify_scalar(&a, &b, &c, 0));
    }

    // SIMD-specific tests
    #[cfg(feature = "simd")]
    #[test]
    fn test_freivalds_simd_u32() {
        let a = vec![vec![1u32, 2], vec![3, 4]];
        let b = vec![vec![2u32, 0], vec![1, 2]];
        let c = vec![vec![4u32, 4], vec![10, 8]]; // c = a * b
        assert!(freivalds_verify_simd::<u32, 4>(&a, &b, &c, 10));
    }

    #[cfg(feature = "simd")]
    #[test]
    fn test_freivalds_simd_u8() {
        let a = vec![vec![1u8, 2], vec![3, 4]];
        let b = vec![vec![2u8, 0], vec![1, 2]];
        let c = vec![vec![4u8, 4], vec![10, 8]]; // c = a * b
        assert!(freivalds_verify_simd::<u8, 4>(&a, &b, &c, 10));
    }

    #[cfg(feature = "simd")]
    #[test]
    fn test_freivalds_simd_u16() {
        let a = vec![vec![1u16, 2], vec![3, 4]];
        let b = vec![vec![2u16, 0], vec![1, 2]];
        let c = vec![vec![4u16, 4], vec![10, 8]]; // c = a * b
        assert!(freivalds_verify_simd::<u16, 4>(&a, &b, &c, 10));
    }

    #[test]
    fn test_freivalds_signed_basic() {
        let a = vec![vec![-1, 2], vec![3, -4]];
        let b = vec![vec![2, 0], vec![1, -2]];
        let c = vec![vec![0, 7], vec![2, 8]]; // c = a * b
        let modulus = 11;
        assert!(freivalds_verify_scalar_mod(&a, &b, &c, 10, modulus));
    }

    #[test]
    fn test_freivalds_signed_negative_modulo() {
        let a = vec![vec![-5, -7], vec![6, -2]];
        let b = vec![vec![3, -1], vec![-4, 2]];
        // Compute c = a * b mod 13
        let mut c = vec![vec![0i64; 2]; 2];
        let modulus = 13;
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    c[i][j] = (c[i][j] + a[i][k] * b[k][j]).rem_euclid(modulus);
                }
            }
        }
        assert!(freivalds_verify_scalar_mod(&a, &b, &c, 10, modulus));
    }

    #[test]
    fn test_freivalds_signed_mixed_sign() {
        let a = vec![vec![1, -2, 3], vec![-4, 5, -6], vec![7, -8, 9]];
        let b = vec![vec![-1, 2, -3], vec![4, -5, 6], vec![-7, 8, -9]];
        let modulus = 17;
        // Compute c = a * b mod 17
        let mut c = vec![vec![0i64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    c[i][j] = (c[i][j] + a[i][k] * b[k][j]).rem_euclid(modulus);
                }
            }
        }
        assert!(freivalds_verify_scalar_mod(&a, &b, &c, 10, modulus));
    }

    #[test]
    fn test_freivalds_signed_all_zeros() {
        let a = vec![vec![0i64, 0], vec![0, 0]];
        let b = vec![vec![0i64, 0], vec![0, 0]];
        let c = vec![vec![0i64, 0], vec![0, 0]];
        let modulus = 5;
        assert!(freivalds_verify_scalar_mod(&a, &b, &c, 10, modulus));
    }

    #[test]
    fn test_freivalds_signed_identity() {
        let a = vec![vec![1i64, 0], vec![0, 1]];
        let b = vec![vec![1i64, 0], vec![0, 1]];
        let c = vec![vec![1i64, 0], vec![0, 1]];
        let modulus = 101;
        assert!(freivalds_verify_scalar_mod(&a, &b, &c, 10, modulus));
    }

    #[test]
    fn test_freivalds_signed_large_negative() {
        let a = vec![vec![-1000000007i64, 2], vec![3, -4]];
        let b = vec![vec![2, 0], vec![1, -2]];
        let modulus = 1_000_000_009i64;
        // Compute c = a * b mod modulus
        let mut c = vec![vec![0i64; 2]; 2];
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    c[i][j] = (c[i][j] + a[i][k] * b[k][j]).rem_euclid(modulus);
                }
            }
        }
        assert!(freivalds_verify_scalar_mod(&a, &b, &c, 10, modulus));
    }

    #[test]
    fn test_freivalds_signed_incorrect() {
        let a = vec![vec![1i64, -2], vec![3, 4]];
        let b = vec![vec![2, 0], vec![1, 2]];
        let c = vec![vec![4, 4], vec![10, 7]]; // incorrect
        let modulus = 13;
        assert!(!freivalds_verify_scalar_mod(&a, &b, &c, 10, modulus));
    }

    #[test]
    fn test_freivalds_signed_k_zero() {
        let a = vec![vec![1i64, 2], vec![3, 4]];
        let b = vec![vec![2, 0], vec![1, 2]];
        let c = vec![vec![4, 4], vec![10, 8]];
        let modulus = 17;
        // k=0: should always return true (no checks performed)
        assert!(freivalds_verify_scalar_mod(&a, &b, &c, 0, modulus));
    }
}
