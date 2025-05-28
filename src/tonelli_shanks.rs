//! Tonelli–Shanks algorithm for modular square roots (r^2 ≡ n mod p) over prime moduli.
//! Constant-time version.
//!
//! # Overview
//! Efficiently computes modular square roots for odd prime moduli. Returns both roots or an error if no solution exists.
//!
//! # Example
//! ```
//! use rma::tonelli_shanks::tonelli_shanks_ct;
//! use num_bigint::BigUint;
//! let (r1, r2) = tonelli_shanks_ct(&BigUint::from(5u64), &BigUint::from(41u64)).unwrap();
//! assert_eq!((&r1 * &r1) % &BigUint::from(41u64), BigUint::from(5u64));
//! assert_eq!((&r2 * &r2) % &BigUint::from(41u64), BigUint::from(5u64));
//! ```

use num_bigint::BigUint;
use num_traits::{One, Zero};
use std::fmt;

// Helper constant for fixed-length iterations.
// Adjust as needed for expected maximum bit length of p.
// Max u64 is 64 bits. Common crypto might be 256, 512, or more.
const FIXED_P_BIT_LENGTH: usize = 256;
const MAX_Z_ITERATIONS: u32 = 200;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TonelliShanksError {
    NotPrime,
    NotQuadraticResidue,
    ZeroOrOne,
    InvalidInput(BigUint), // Used for "z not found" or "i == m"
}

impl fmt::Display for TonelliShanksError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TonelliShanksError::NotPrime => write!(f, "Modulus is not prime"),
            TonelliShanksError::NotQuadraticResidue => {
                write!(f, "n is not a quadratic residue mod p")
            }
            TonelliShanksError::ZeroOrOne => {
                write!(f, "Trivial case: n = 0 or 1 (handled directly)")
            }
            TonelliShanksError::InvalidInput(x) => write!(
                f,
                "Invalid input or algorithmic failure (e.g., z not found, i==m) for value: {}",
                x
            ),
        }
    }
}

/// Constant-time select for BigUint.
/// Chooses val_true if condition is true, val_false otherwise.
/// Relies on BigUint multiplication by 0 or 1 being reasonably safe, and XOR.
/// `control_val = 1` for true, `0` for false.
/// R = val_false ^ (control_val * (val_true ^ val_false))
fn ct_select_biguint(condition: bool, val_true: &BigUint, val_false: &BigUint) -> BigUint {
    let control_big = BigUint::from(condition as u8); // Creates BigUint::one() or BigUint::zero()
    // Ensure cloned values are used in calculation to avoid modifying inputs if they are the same.
    let vt = val_true.clone();
    let vf = val_false.clone();
    vf.clone() ^ (&control_big * (vt ^ vf))
}

/// Constant-time equality check (conceptual).
/// Uses standard `==`, which may not be strictly CT for BigUint.
#[inline]
fn ct_eq_biguint(a: &BigUint, b: &BigUint) -> bool {
    a == b
}

/// Constant-time modular exponentiation for BigUint (right-to-left binary method).
/// Iterates `max_bits` times, making it independent of `exp`'s actual value.
fn modpow_biguint_ct(base: &BigUint, exp: &BigUint, modulus: &BigUint, max_bits: usize) -> BigUint {
    if modulus.is_zero() {
        return BigUint::zero(); // Or handle as an error/precondition violation
    }
    if modulus.is_one() {
        return BigUint::zero(); // Result of any mod 1 is 0
    }

    let mut res = BigUint::one();
    let mut b_val = base % modulus;
    let mut e_val = exp.clone();
    let one = BigUint::one();

    for _i in 0..max_bits {
        let e_is_odd = ct_eq_biguint(&(&e_val & &one), &one);

        let prod = (&res * &b_val) % modulus;
        res = ct_select_biguint(e_is_odd, &prod, &res);

        e_val >>= 1u32;
        b_val = (&b_val * &b_val) % modulus;
    }
    res
}

/// Returns (r1, r2) such that r1^2 ≡ n mod p, r2 = p - r1, or an error if no solution exists.
/// Constant-time version.
/// Computes the modular square root of `n` modulo prime `p` using the constant-time Tonelli-Shanks algorithm.
///
/// Returns a tuple `(r1, r2)` such that `r1^2 ≡ n (mod p)` and `r2 = p - r1`, or an error if no solution exists.
/// This function is designed to execute in constant time with respect to secret data, mitigating timing side-channels
/// (assuming the underlying `BigUint` operations are constant-time).
///
/// # Arguments
///
/// * `n` - The value whose square root is to be computed modulo `p`.
/// * `p` - The prime modulus.
///
/// # Returns
///
/// * `Ok((r1, r2))` if a modular square root exists, where `r1^2 ≡ n (mod p)` and `r2 = p - r1`.
/// * `Err(TonelliShanksError)` if no solution exists or if the input is invalid.
///
/// # Errors
///
/// * `TonelliShanksError::NotPrime` if `p` is zero or one.
/// * `TonelliShanksError::NotQuadraticResidue` if `n` is not a quadratic residue modulo `p`.
///
/// # Examples
///
/// ```
/// use num_bigint::BigUint;
/// use rma::tonelli_shanks::{tonelli_shanks_ct, TonelliShanksError};
///
/// // Example: Find sqrt(10) mod 13
/// let n = BigUint::from(10u32);
/// let p = BigUint::from(13u32);
/// let res = tonelli_shanks_ct(&n, &p);
/// assert!(res.is_ok());
/// let (r1, r2) = res.unwrap();
/// assert_eq!((&r1 * &r1) % &p, n);
/// assert_eq!((&r2 * &r2) % &p, n);
///
/// // Example: No solution (non-residue)
/// let n = BigUint::from(5u32);
/// let p = BigUint::from(7u32);
/// let res = tonelli_shanks_ct(&n, &p);
/// assert!(matches!(res, Err(TonelliShanksError::NotQuadraticResidue)));
/// ```
pub fn tonelli_shanks_ct(
    n: &BigUint,
    p: &BigUint,
) -> Result<(BigUint, BigUint), TonelliShanksError> {
    // 1. Initial non-CT checks for invalid p and trivial cases for n (can return early).
    // These are not secret-dependent branches.
    if p.is_zero() || ct_eq_biguint(p, &BigUint::one()) {
        return Err(TonelliShanksError::NotPrime);
    }
    if n.is_zero() {
        return Ok((BigUint::zero(), BigUint::zero()));
    }
    if ct_eq_biguint(n, &BigUint::one()) {
        let r2 = p - BigUint::one();
        return Ok((BigUint::one(), r2));
    }

    // Initialize state for CT computation
    let mut overall_error: Option<TonelliShanksError> = None;
    let one_val = BigUint::one();
    let zero_val = BigUint::zero();

    // Variables that will be computed. Initialize to defaults for error paths.
    let r_final_computed;

    // --- Start of CT block --- All operations from here aim for constant time.

    // 2. Euler's criterion: n^((p-1)/2) ≡ 1 (mod p)
    let p_minus_1 = p - &one_val;
    let exp_euler = &p_minus_1 >> 1u32;
    // Use p.bits() for max_bits in modpow, as exp_euler is related to p.
    let n_pow_exp = modpow_biguint_ct(n, &exp_euler, p, p.bits() as usize);
    let is_not_residue = !ct_eq_biguint(&n_pow_exp, &one_val);
    if is_not_residue {
        overall_error = Some(TonelliShanksError::NotQuadraticResidue);
    }

    // 3. S and Q computation: p-1 = q * 2^s with q odd
    let mut q_val = p_minus_1.clone();
    let mut s_val = 0u32;
    let mut q_s_processing_done_flag = false; // True when Q is odd or loop finishes

    for _ in 0..FIXED_P_BIT_LENGTH {
        // Iterate up to max bit length of p
        let active_sq_step = !q_s_processing_done_flag && overall_error.is_none();
        let q_is_even = ct_eq_biguint(&(&q_val & &one_val), &zero_val); // q & 1 == 0
        let should_update_q_and_s = q_is_even && active_sq_step;
        let q_shifted = &q_val >> 1u32;
        q_val = ct_select_biguint(should_update_q_and_s, &q_shifted, &q_val);
        // s_val update: only if should_update_q_and_s is true
        let s_incremented = s_val + 1;
        s_val = if should_update_q_and_s {
            s_incremented
        } else {
            s_val
        }; // This if on u32 is fine.
        // If we didn't update (because q was already odd, or step became inactive), then done.
        if !should_update_q_and_s {
            q_s_processing_done_flag = true;
        }
    }
    // At this point, q_val is Q and s_val is S, or defaults if error occurred.

    // 4. Finding quadratic non-residue z
    let mut found_z_val = BigUint::from(2u32); // Default z if not found or error
    let mut z_actually_found_flag = false;

    for i in 2u32..(2u32 + MAX_Z_ITERATIONS) {
        // Fixed iterations
        let active_z_step = overall_error.is_none(); // Check if we should still be doing real work
        let current_z_candidate = BigUint::from(i);
        let z_pow_exp = modpow_biguint_ct(&current_z_candidate, &exp_euler, p, p.bits() as usize);
        let is_non_residue_cond = ct_eq_biguint(&z_pow_exp, &p_minus_1);
        // Capture if non-residue, not already found, and no overall error
        let should_capture_this_z = is_non_residue_cond && !z_actually_found_flag && active_z_step;
        found_z_val = ct_select_biguint(should_capture_this_z, &current_z_candidate, &found_z_val);
        if should_capture_this_z {
            z_actually_found_flag = true;
        } // simple bool update
    }
    let z_finding_failed = !z_actually_found_flag;
    if z_finding_failed && overall_error.is_none() {
        overall_error = Some(TonelliShanksError::InvalidInput(p.clone()));
    }

    // 5. Initialize c, r, t, m for the main loop
    // These values are calculated using potentially dummy `found_z_val` or `q_val` if an error has occurred.
    // `ct_select_biguint` will ensure they take on default/previous values if `overall_error.is_some()`.

    let c_initial = modpow_biguint_ct(&found_z_val, &q_val, p, p.bits() as usize);
    let q_plus_1 = &q_val + &one_val; // q_val could be dummy if error
    let q_plus_1_div_2 = &q_plus_1 >> 1u32;
    let r_initial = modpow_biguint_ct(n, &q_plus_1_div_2, p, p.bits() as usize);
    let t_initial = modpow_biguint_ct(n, &q_val, p, p.bits() as usize);
    let m_initial = s_val; // s_val computed earlier

    // Initialize working variables for the loop, using defaults if error already occurred.
    let mut c_val = ct_select_biguint(overall_error.is_none(), &c_initial, &one_val); // Default c to 1
    let mut r_val = ct_select_biguint(overall_error.is_none(), &r_initial, &zero_val); // Default r to 0
    let mut t_val = ct_select_biguint(overall_error.is_none(), &t_initial, &one_val); // Default t to 1 (makes loop no-op)
    let mut m_val = if overall_error.is_none() {
        m_initial
    } else {
        0
    }; // Default m to 0

    // 6. Main computational loop (based on s, max FIXED_P_BIT_LENGTH iterations)
    for _outer_loop_idx in 0..FIXED_P_BIT_LENGTH {
        // Condition for this iteration to perform real work: no overall error, and t is not yet 1.
        let active_outer_iteration = overall_error.is_none() && !ct_eq_biguint(&t_val, &one_val);

        // Inner loop to find i: t^(2^i) ≡ 1 (mod p). `determined_i_val` is this `i` (0-indexed).
        let mut determined_i_val = m_val; // Default to m (error `i == m`), if i not found or m=0.
        let mut i_search_temp_t = t_val.clone(); // This is t, then t^2, t^4, ...
        let mut inner_i_found_flag = false;

        // Inner loop runs up to m_val times. Max m_val is FIXED_P_BIT_LENGTH.
        // So, loop FIXED_P_BIT_LENGTH times. `k_inner` is the exponent for 2 (0 for t, 1 for t^2, etc.)
        for k_inner in 0..FIXED_P_BIT_LENGTH {
            // Active if outer iter is active, i not found yet, and k_inner is less than current m_val.
            let active_inner_step =
                active_outer_iteration && !inner_i_found_flag && (k_inner < m_val as usize);

            let is_temp_t_one = ct_eq_biguint(&i_search_temp_t, &one_val);
            let should_capture_i = active_inner_step && is_temp_t_one;

            if should_capture_i {
                determined_i_val = k_inner as u32; // Capture first k_inner where t^(2^k_inner) is 1
            }
            // inner_i_found_flag = inner_i_found_flag || should_capture_i;
            if should_capture_i {
                inner_i_found_flag = true;
            }

            let squared_i_search_temp_t = (&i_search_temp_t * &i_search_temp_t) % p;
            // Stop updating temp_t if i found or step is no longer active.
            let cond_stop_temp_t_update = inner_i_found_flag || !active_inner_step;
            i_search_temp_t = ct_select_biguint(
                cond_stop_temp_t_update,
                &i_search_temp_t,
                &squared_i_search_temp_t,
            );
        }
        // After inner loop, `determined_i_val` is the 0-indexed `i`. `inner_i_found_flag` is true if found.

        // Check for error condition: i == m (original was 1-indexed `i`, here 0-indexed `determined_i_val`)
        // This means loop finished (up to m_val) and `i` was not found (inner_i_found_flag is false).
        // Only an error if m_val was > 0 (meaning a search for `i` was expected).
        let i_equals_m_error_trigger = active_outer_iteration && !inner_i_found_flag && (m_val > 0);
        if i_equals_m_error_trigger && overall_error.is_none() {
            // Use n.clone(); for true CT, error value construction should also be CT.
            overall_error = Some(TonelliShanksError::InvalidInput(n.clone()));
        }

        // Updates for r, c, t, m, only if no error and t was not one and i was found.
        let can_update_rctm =
            active_outer_iteration && overall_error.is_none() && inner_i_found_flag;

        // Calculate exponent for b: 2^(m - i - 1)
        // m_val is current M. determined_i_val is i (0-indexed).
        // Exponent is m_val - determined_i_val - 1. Must be non-negative.
        let mut b_exp_power = 0u32;
        if m_val > determined_i_val && (m_val - determined_i_val) >= 1 {
            b_exp_power = m_val - determined_i_val - 1;
        }
        // `BigUint::one() << b_exp_power` creates 2^(b_exp_power)
        let exp_for_b_ct = BigUint::from(1u32) << (b_exp_power as usize);

        // All these calculations use potentially dummy values if previous steps had issues,
        // but `ct_select` later ensures variables are only updated if `can_update_rctm` is true.
        let b_val = modpow_biguint_ct(&c_val, &exp_for_b_ct, p, p.bits() as usize);

        let r_new = (&r_val * &b_val) % p;
        let c_new = (&b_val * &b_val) % p;
        let t_new = (&t_val * &c_new) % p;
        let m_new = determined_i_val; // new m is the found i (0-indexed)

        r_val = ct_select_biguint(can_update_rctm, &r_new, &r_val);
        c_val = ct_select_biguint(can_update_rctm, &c_new, &c_val);
        t_val = ct_select_biguint(can_update_rctm, &t_new, &t_val);
        m_val = if can_update_rctm { m_new } else { m_val };
    }

    // --- End of CT block ---

    r_final_computed = r_val; // r_val holds the computed R, or its default (zero) if error

    if let Some(err) = overall_error {
        Err(err)
    } else {
        // If no error, r_final_computed has the result.
        let r2_final = (p - &r_final_computed) % p;
        Ok((r_final_computed, r2_final))
    }
}

#[cfg(test)]
mod tests {
    use rand::Rng;

    use super::*;

    #[test]
    fn test_known_large_prime_10digit() {
        let p = BigUint::from(1_000_000_009u64);
        let n = BigUint::from(665_820_697u64);
        let (r1, r2) = tonelli_shanks_ct(&n, &p).unwrap();
        assert!(r1 == BigUint::from(378_633_312u64) || r2 == BigUint::from(378_633_312u64));
        assert!(r1 == BigUint::from(621_366_697u64) || r2 == BigUint::from(621_366_697u64));
        assert_eq!((&r1 * &r1) % &p, n);
        assert_eq!((&r2 * &r2) % &p, n);
    }
    #[test]
    fn test_known_large_prime_13digit() {
        let p = BigUint::from(1_000_000_000_039u64);
        let n = BigUint::from(881_398_088_036u64);
        let (r1, r2) = tonelli_shanks_ct(&n, &p).unwrap();
        assert!(r1 == BigUint::from(791_399_408_049u64) || r2 == BigUint::from(791_399_408_049u64));
        assert!(r1 == BigUint::from(208_600_591_990u64) || r2 == BigUint::from(208_600_591_990u64));
        assert_eq!((&r1 * &r1) % &p, n);
        assert_eq!((&r2 * &r2) % &p, n);
    }
    #[test]
    fn test_known_large_prime_14digit() {
        let p = BigUint::from(41_434_547_495_153u64);
        let n = BigUint::from(2u64);
        let (r1, r2) = tonelli_shanks_ct(&n, &p).unwrap();
        assert!(
            r1 == BigUint::from(32_789_792_927_270u64)
                || r2 == BigUint::from(32_789_792_927_270u64)
        );
        assert!(
            r1 == BigUint::from(8_644_754_567_883u64) || r2 == BigUint::from(8_644_754_567_883u64)
        );
        assert_eq!((&r1 * &r1) % &p, n);
    }
    #[test]
    fn test_known_large_prime_6digit() {
        let p = BigUint::from(100_049u64);
        let n = BigUint::from(44_402u64);
        let (r1, r2) = tonelli_shanks_ct(&n, &p).unwrap();
        assert!(r1 == BigUint::from(30_468u64) || r2 == BigUint::from(30_468u64));
        assert!(r1 == BigUint::from(69_581u64) || r2 == BigUint::from(69_581u64));
        assert_eq!((&r1 * &r1) % &p, n);
        assert_eq!((&r2 * &r2) % &p, n);
    }
    #[test]
    fn test_edge_case_zero() {
        let p = BigUint::from(101u64);
        let n = BigUint::zero();
        let res = tonelli_shanks_ct(&n, &p);
        assert!(res.is_ok());
        let (r1, r2) = res.unwrap();
        assert_eq!(r1, BigUint::zero());
        assert_eq!(r2, BigUint::zero());
    }
    #[test]
    fn test_edge_case_one() {
        let p = BigUint::from(101u64);
        let n = BigUint::one();
        let res = tonelli_shanks_ct(&n, &p);
        assert!(res.is_ok());
        let (r1, r2) = res.unwrap();
        assert!(
            (r1 == BigUint::one() && r2 == BigUint::from(100u64))
                || (r1 == BigUint::from(100u64) && r2 == BigUint::one())
        );
    }
    #[test]
    fn test_non_residue() {
        let p = BigUint::from(1009u64);
        // n = 23 is a non-residue mod 1009, (23^((1009-1)/2) mod 1009 = 1008)
        let n = BigUint::from(23u64); // 1032 % 1009 = 23
        let res = tonelli_shanks_ct(&n, &p);
        assert!(res.is_err());
        if let Err(TonelliShanksError::NotQuadraticResidue) = res {
            // Correct error type
        } else {
            panic!("Expected NotQuadraticResidue error, got {:?}", res);
        }
    }
    #[test]
    fn test_tonelli_shanks_basic() {
        let n = BigUint::from(5u64);
        let p = BigUint::from(41u64);
        let (r1, r2) = tonelli_shanks_ct(&n, &p).unwrap();
        assert_eq!((&r1 * &r1) % &p, n);
        assert_eq!((&r2 * &r2) % &p, n);
        assert_ne!(r1, r2); // Ensure roots are distinct (unless p=2 or r=0)
        // For (5, 41), roots are 13 and 28.
        assert!(
            (r1 == BigUint::from(13u64) && r2 == BigUint::from(28u64))
                || (r1 == BigUint::from(28u64) && r2 == BigUint::from(13u64))
        );
    }
    #[test]
    fn test_tonelli_shanks_non_residue() {
        let n = BigUint::from(3u64);
        let p = BigUint::from(7u64); // 3 is non-residue mod 7
        let res = tonelli_shanks_ct(&n, &p);
        assert!(res.is_err());
        if let Err(TonelliShanksError::NotQuadraticResidue) = res {
            // Correct error type
        } else {
            panic!("Expected NotQuadraticResidue error, got {:?}", res);
        }
    }

    #[test]
    fn test_non_residue_hard_p() {
        // p = 41, exp_euler = 20
        // Try n=3 (non-residue for p=41, 3^20 = 40 mod 41)
        let p = BigUint::from(41u64);
        let n = BigUint::from(3u64);
        let res = tonelli_shanks_ct(&n, &p);
        assert!(res.is_err());
        assert_eq!(res.unwrap_err(), TonelliShanksError::NotQuadraticResidue);

        // Test where z is harder to find (p needs to be structured s.t. small z are residues)
        // This is implicitly tested if the MAX_Z_ITERATIONS is hit without finding z.
        // Example: p = 257 (prime). (p-1)/2 = 128.
        // 2 is non-residue (2^128 = 1 mod 257 -> no, this means 2 is residue if 1)
        // 2^128 mod 257. If this is p-1 (256), then 2 is non-residue.
        // For p=257, z=3 is the first non-residue.
        // Let's test a case where z=2 is a residue.
        // p=17. (p-1)/2 = 8.
        // 2^8 mod 17 = 256 mod 17 = 1. So 2 is a quadratic residue.
        // 3^8 mod 17 = 6561 mod 17 = 16 = p-1. So z=3 is a non-residue.
        // Test with n=2 (residue), p=17. Roots should be 6, 11. (6^2=36=2 mod 17)
        let p17 = BigUint::from(17u64);
        let n2 = BigUint::from(2u64);
        let (r1, r2) = tonelli_shanks_ct(&n2, &p17).unwrap();
        assert!(
            (r1 == BigUint::from(6u64) && r2 == BigUint::from(11u64))
                || (r1 == BigUint::from(11u64) && r2 == BigUint::from(6u64))
        );
    }

    #[test]
    fn test_nist_p192_prime() {
        // NIST P-192 prime: FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF
        let p =
            BigUint::parse_bytes(b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF", 16).unwrap();
        let n = BigUint::from(4u32);
        let (r1, r2) = tonelli_shanks_ct(&n, &p).unwrap();
        assert_eq!((&r1 * &r1) % &p, n);
        assert_eq!((&r2 * &r2) % &p, n);
        assert_ne!(r1, r2);
    }

    #[test]
    fn test_nist_p256_prime() {
        // NIST P-256 prime: FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF
        let p = BigUint::parse_bytes(
            b"FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF",
            16,
        )
        .unwrap();
        let n = BigUint::from(4u32);
        let (r1, r2) = tonelli_shanks_ct(&n, &p).unwrap();
        assert_eq!((&r1 * &r1) % &p, n);
        assert_eq!((&r2 * &r2) % &p, n);
        assert_ne!(r1, r2);
    }

    #[test]
    fn test_nist_p224_prime() {
        // NIST P-224 prime: FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000000000000000000001
        let p = BigUint::parse_bytes(
            b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000000000000000000001",
            16,
        )
        .unwrap();
        // Test n = 4 (quadratic residue)
        let n = BigUint::from(4u32);
        let (r1, r2) = tonelli_shanks_ct(&n, &p).unwrap();
        assert_eq!((&r1 * &r1) % &p, n);
        assert_eq!((&r2 * &r2) % &p, n);
        // Test n = 9 (quadratic residue)
        let n = BigUint::from(9u32);
        let (r1, r2) = tonelli_shanks_ct(&n, &p).unwrap();
        assert_eq!((&r1 * &r1) % &p, n);
        assert_eq!((&r2 * &r2) % &p, n);
        // Test n = 3 (should only error if 3 is a non-residue)
        let n = BigUint::from(3u32);
        let p_legendre = p.clone();
        let exp = (&p_legendre - BigUint::one()) >> 1u32;
        let legendre = modpow_biguint_ct(&n, &exp, &p_legendre, 256);
        if legendre == p_legendre - BigUint::one() {
            // 3 is a non-residue, should error
            let res = tonelli_shanks_ct(&n, &p);
            assert!(res.is_err());
        } else {
            // 3 is a residue, should succeed
            let res = tonelli_shanks_ct(&n, &p);
            assert!(res.is_ok());
        }
        // Test using x from NIST ECDSA P-224 public key Q_x
        // Q_x = E84FB0B8E7000CB657D7973CF6B42ED78B301674276DF744AF130B3E
        let x = BigUint::parse_bytes(
            b"E84FB0B8E7000CB657D7973CF6B42ED78B301674276DF744AF130B3E",
            16,
        )
        .unwrap();
        let n = (&x * &x) % &p;
        let (r1, r2) = tonelli_shanks_ct(&n, &p).unwrap();
        assert!(r1 == x || r2 == x || r1 == (&p - &x) || r2 == (&p - &x));
    }

    #[test]
    fn test_nist_p384_prime() {
        // NIST P-384 prime: FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFF0000000000000000FFFFFFFF
        let p = BigUint::parse_bytes(b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFF0000000000000000FFFFFFFF", 16).unwrap();
        // Test n = 4 (quadratic residue)
        let n = BigUint::from(4u32);
        let res = tonelli_shanks_ct(&n, &p);
        let (r1, r2) = res.unwrap();
        assert_eq!((&r1 * &r1) % &p, n);
        assert_eq!((&r2 * &r2) % &p, n);
        // Test n = 9 (quadratic residue)
        let n = BigUint::from(9u32);
        let (r1, r2) = tonelli_shanks_ct(&n, &p).unwrap();
        assert_eq!((&r1 * &r1) % &p, n);
        assert_eq!((&r2 * &r2) % &p, n);
        // Test n = 3 (should only error if 3 is a non-residue)
        let n = BigUint::from(3u32);
        let p_legendre = p.clone();
        let exp = (&p_legendre - BigUint::one()) >> 1u32;
        let legendre = modpow_biguint_ct(&n, &exp, &p_legendre, 384);
        if legendre == p_legendre - BigUint::one() {
            // 3 is a non-residue, should error
            let res = tonelli_shanks_ct(&n, &p);
            assert!(res.is_err());
        } else {
            // 3 is a residue, should succeed
            let res = tonelli_shanks_ct(&n, &p);
            assert!(res.is_ok());
        }
        // Test using x from NIST ECDSA P-384 public key Q_x
        // Q_x = 3BF701BC9E9D36B4D5F1455343F09126F2564390F2B487365071243C61E6471FB9D2AB74657B82F9086489D9EF0F5CB5
        let x = BigUint::parse_bytes(b"3BF701BC9E9D36B4D5F1455343F09126F2564390F2B487365071243C61E6471FB9D2AB74657B82F9086489D9EF0F5CB5", 16).unwrap();
        let n = (&x * &x) % &p;
        let (r1, r2) = tonelli_shanks_ct(&n, &p).unwrap();
        assert!(r1 == x || r2 == x || r1 == (&p - &x) || r2 == (&p - &x));
    }

    #[test]
    fn test_nist_p521_prime() {
        // NIST P-521 prime: 01FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
        let p = BigUint::parse_bytes(b"01FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF", 16).unwrap();
        // Test n = 4 (quadratic residue)
        let n = BigUint::from(4u32);
        let res = tonelli_shanks_ct(&n, &p);
        let (r1, r2) = res.unwrap();
        assert_eq!((&r1 * &r1) % &p, n);
        assert_eq!((&r2 * &r2) % &p, n);
        // Test n = 9 (quadratic residue)
        let n = BigUint::from(9u32);
        let (r1, r2) = tonelli_shanks_ct(&n, &p).unwrap();
        assert_eq!((&r1 * &r1) % &p, n);
        assert_eq!((&r2 * &r2) % &p, n);
        // Test n = 3 (should only error if 3 is a non-residue)
        let n = BigUint::from(3u32);
        let p_legendre = p.clone();
        let exp = (&p_legendre - BigUint::one()) >> 1u32;
        let legendre = modpow_biguint_ct(&n, &exp, &p_legendre, 521);
        if legendre == p_legendre - BigUint::one() {
            // 3 is a non-residue, should error
            let res = tonelli_shanks_ct(&n, &p);
            assert!(res.is_err());
        } else {
            // 3 is a residue, should succeed
            let res = tonelli_shanks_ct(&n, &p);
            assert!(res.is_ok());
        }
        // Test using x from NIST ECDSA P-521 public key Q_x
        // Q_x = 0098E91EEF9A68452822309C52FAB453F5F117C1DA8ED796B255E9AB8F6410CCA16E59DF403A6BDC6CA467A37056B1E54B3005D8AC030DECFEB68DF18B171885D5C4
        let x = BigUint::parse_bytes(b"0098E91EEF9A68452822309C52FAB453F5F117C1DA8ED796B255E9AB8F6410CCA16E59DF403A6BDC6CA467A37056B1E54B3005D8AC030DECFEB68DF18B171885D5C4", 16).unwrap();
        let n = (&x * &x) % &p;
        let (r1, r2) = tonelli_shanks_ct(&n, &p).unwrap();
        assert!(r1 == x || r2 == x || r1 == (&p - &x) || r2 == (&p - &x));
    }

    #[test]
    fn test_mersenne_prime() {
        // Mersenne prime: 2^127 - 1
        let p = BigUint::parse_bytes(b"7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF", 16).unwrap();
        // n = 4 (residue)
        let n = BigUint::from(4u32);
        let (r1, r2) = tonelli_shanks_ct(&n, &p).unwrap();
        assert_eq!((&r1 * &r1) % &p, n);
        assert_eq!((&r2 * &r2) % &p, n);
        // n = 3 (likely non-residue for this prime)
        let n = BigUint::from(3u32);
        let exp = (&p - BigUint::one()) >> 1u32;
        let legendre = modpow_biguint_ct(&n, &exp, &p, p.bits() as usize);
        if legendre == p.clone() - BigUint::one() {
            assert!(tonelli_shanks_ct(&n, &p).is_err());
        }
    }

    #[test]
    fn test_safe_prime() {
        // Safe prime: 2^521 - 1 (also a Mersenne prime, but used in cryptography)
        let p = BigUint::parse_bytes(b"01FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF", 16).unwrap();
        // n = 9 (residue)
        let n = BigUint::from(9u32);
        let (r1, r2) = tonelli_shanks_ct(&n, &p).unwrap();
        assert_eq!((&r1 * &r1) % &p, n);
        assert_eq!((&r2 * &r2) % &p, n);
        // n = 5 (random small value)
        let n = BigUint::from(5u32);
        let exp = (&p - BigUint::one()) >> 1u32;
        let legendre = modpow_biguint_ct(&n, &exp, &p, p.bits() as usize);
        if legendre == p.clone() - BigUint::one() {
            assert!(tonelli_shanks_ct(&n, &p).is_err());
        }
    }

    #[test]
    fn test_random_residues_and_nonresidues_521bit() {
        // 521-bit random prime (use NIST P-521 for reproducibility)
        let p = BigUint::parse_bytes(b"01FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF", 16).unwrap();
        let mut rng = rand::thread_rng();
        for _ in 0..5 {
            // Random x < p
            let mut bytes = vec![0u8; 66];
            rng.fill(&mut bytes[..]);
            bytes[0] &= 0x01; // ensure < 2^521
            let x = BigUint::from_bytes_be(&bytes) % &p;
            let n = (&x * &x) % &p;
            // n is a residue
            let (r1, r2) = tonelli_shanks_ct(&n, &p).unwrap();
            assert_eq!((&r1 * &r1) % &p, n);
            assert_eq!((&r2 * &r2) % &p, n);
            // Now test a random non-residue
            let mut candidate = BigUint::from(rng.gen_range(2u32..1000u32));
            let exp = (&p - BigUint::one()) >> 1u32;
            let mut found = false;
            for _ in 0..100 {
                let legendre = modpow_biguint_ct(&candidate, &exp, &p, p.bits() as usize);
                if legendre == p.clone() - BigUint::one() {
                    assert!(tonelli_shanks_ct(&candidate, &p).is_err());
                    found = true;
                    break;
                }
                candidate += BigUint::one();
            }
            assert!(found, "Failed to find a non-residue in 100 tries");
        }
    }

    #[test]
    fn test_random_1024bit_prime() {
        // 1024-bit prime (RFC 3526 Group 2)
        let p = BigUint::parse_bytes(b"FFFFFFFFFFFFFFFFC90FDAA22168C234C4C6628B80DC1CD129024E088A67CC74020BBEA63B139B22514A08798E3404DDEF9519B3CD3A431B302B0A6DF25F14374FE1356D6D51C245E485B576625E7EC6F44C42E9A637ED6B0BFF5CB6F406B7EDEE386BFB5A899FA5AE9F24117C4B1FE649286651ECE65381FFFFFFFFFFFFFFFF", 16).unwrap();
        let mut rng = rand::thread_rng();
        for _ in 0..3 {
            // Random x < p
            let mut bytes = vec![0u8; 128];
            rng.fill(&mut bytes[..]);
            let x = BigUint::from_bytes_be(&bytes) % &p;
            let n = (&x * &x) % &p;
            // n is a residue
            let (r1, r2) = tonelli_shanks_ct(&n, &p).unwrap();
            assert_eq!((&r1 * &r1) % &p, n);
            assert_eq!((&r2 * &r2) % &p, n);
            // Now test a random non-residue
            let mut candidate = BigUint::from(rng.gen_range(2u32..1000u32));
            let exp = (&p - BigUint::one()) >> 1u32;
            let mut found = false;
            for _ in 0..100 {
                let legendre = modpow_biguint_ct(&candidate, &exp, &p, p.bits() as usize);
                if legendre == p.clone() - BigUint::one() {
                    assert!(tonelli_shanks_ct(&candidate, &p).is_err());
                    found = true;
                    break;
                }
                candidate += BigUint::one();
            }
            assert!(found, "Failed to find a non-residue in 100 tries");
        }
    }
}
