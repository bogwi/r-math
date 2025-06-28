//! Chakravala method for solving Pell's equation x² - d*y² = 1.
//!
//! This module implements the Chakravala method for finding the fundamental
//! (smallest) integer solution to Pell's equation.
//!
//! # Example
//! ```
//! use rma::chakravala::chakravala;
//! use num_bigint::ToBigUint;
//! // d = 2, fundamental solution is (3, 2)
//! let d2 = 2u32.to_biguint().unwrap();
//! let (x2, y2) = chakravala(&d2).unwrap();
//! assert_eq!(x2, 3u32.to_biguint().unwrap());
//! assert_eq!(y2, 2u32.to_biguint().unwrap());
//! ``` 

use num_bigint::{BigUint, ToBigInt};
use num_traits::{One, Signed, Zero};

/// Solves Pell's equation x² - d*y² = 1 for a given non-square integer d.
///
/// The Chakravala method is an iterative algorithm for finding the fundamental
/// (smallest) integer solution to Pell's equation.
///
/// # Arguments
///
/// * `d` - A non-square positive integer.
///
/// # Returns
///
/// An `Option` containing a tuple `(x, y)` which is the fundamental solution,
/// or `None` if `d` is a perfect square.
///
/// # Examples
///
/// ```
/// use num_bigint::ToBigUint;
/// use rma::chakravala::chakravala;
///
/// // d = 13, fundamental solution is (649, 180)
/// let d = 13u32.to_biguint().unwrap();
/// let (x, y) = chakravala(&d).unwrap();
/// assert_eq!(x, 649u32.to_biguint().unwrap());
/// assert_eq!(y, 180u32.to_biguint().unwrap());
///
/// // d = 64 (perfect square), should return None
/// let d_sq = 64u32.to_biguint().unwrap();
/// assert!(chakravala(&d_sq).is_none());
/// ```
pub fn chakravala(d: &BigUint) -> Option<(BigUint, BigUint)> {
    let d_sqrt = d.sqrt();
    if &d_sqrt * &d_sqrt == *d {
        return None; // d is a perfect square, no non-trivial solution
    }

    let d_bigint = d.to_bigint().unwrap();

    let mut a = d_sqrt.clone();
    let mut b = BigUint::one();
    let mut k = (&a * &a).to_bigint().unwrap() - &d_bigint;

    if k.is_zero() {
        // This case should not be reached due to the perfect square check, but for safety:
        return None;
    }

    while !k.is_one() {
        let k_abs_bigint = k.abs();
        let k_abs_biguint = k_abs_bigint.to_biguint().unwrap();
        let m_congruence = (&a * &b.modinv(&k_abs_biguint).unwrap()).to_bigint().unwrap() % &k;
        let m_neg = -m_congruence;

        let mut m = m_neg.clone();
        if m.is_negative() || m.is_zero() {
            m += &k_abs_bigint;
        }

        let mut best_m = m.clone();
        let mut min_diff = (&best_m * &best_m - &d_bigint).abs();

        loop {
            let next_m = &m + &k_abs_bigint;
            let diff = (&next_m * &next_m - &d_bigint).abs();
            if diff < min_diff {
                min_diff = diff;
                best_m = next_m;
            } else {
                break;
            }
        }
        
        let m = best_m;
        
        let new_a = (&a.to_bigint().unwrap() * &m + &d_bigint * b.to_bigint().unwrap()) / &k_abs_bigint;
        let new_b = (&a.to_bigint().unwrap() + b.to_bigint().unwrap() * &m) / &k_abs_bigint;
        let new_k = (m.pow(2) - &d_bigint) / &k;

        a = new_a.to_biguint().unwrap();
        b = new_b.to_biguint().unwrap();
        k = new_k;
    }

    Some((a, b))
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::ToBigUint;

    #[test]
    fn test_chakravala_d_13() {
        let d = 13.to_biguint().unwrap();
        let (x, y) = chakravala(&d).unwrap();
        assert_eq!(x, 649.to_biguint().unwrap());
        assert_eq!(y, 180.to_biguint().unwrap());
    }

    #[test]
    fn test_chakravala_d_61() {
        let d = 61.to_biguint().unwrap();
        let (x, y) = chakravala(&d).unwrap();
        assert_eq!(x, 1766319049.to_biguint().unwrap());
        assert_eq!(y, 226153980.to_biguint().unwrap());
    }

    #[test]
    fn test_chakravala_d_67() {
        let d = 67.to_biguint().unwrap();
        let (x, y) = chakravala(&d).unwrap();
        assert_eq!(x, 48842.to_biguint().unwrap());
        assert_eq!(y, 5967.to_biguint().unwrap());
    }

    #[test]
    fn test_chakravala_perfect_square() {
        let d = 64.to_biguint().unwrap();
        assert!(chakravala(&d).is_none());
    }

    #[test]
    fn test_chakravala_d_2() {
        let d = 2.to_biguint().unwrap();
        let (x, y) = chakravala(&d).unwrap();
        assert_eq!(x, 3.to_biguint().unwrap());
        assert_eq!(y, 2.to_biguint().unwrap());
    }
} 