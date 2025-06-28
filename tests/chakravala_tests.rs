use num_bigint::{BigUint, ToBigUint};
use num_traits::One;
use num_integer::Integer;
use rma::chakravala::chakravala;

// ---------------- Property-based correctness ----------------
#[cfg(test)]
mod prop {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        /// For randomly chosen non-square d in [2, 1000), the algorithm must
        /// return (x,y) such that x^2 − d·y^2 = 1 and gcd(x,y) = 1.
        #[test]
        fn solves_pell_equation(d in 2u32..1000) {
            let sqrt = (d as f64).sqrt() as u32;
            // Skip perfect squares
            prop_assume!(sqrt * sqrt != d);

            let d_big = d.to_biguint().unwrap();
            let (x, y) = chakravala(&d_big).expect("Expected a solution");

            // Check Pell equation
            let lhs = &x * &x - &d_big * &y * &y;
            prop_assert_eq!(lhs, BigUint::one());

            // Fundamental solution should be primitive (coprime)
            prop_assert_eq!(x.gcd(&y), BigUint::one());
        }
    }
}

// ---------------- Exhaustive fundamental check for small d ----------------
#[test]
fn fundamental_solution_small_d() {
    for d in 2u64..=50 {
        let sqrt = (d as f64).sqrt() as u64;
        if sqrt * sqrt == d { continue; }

        let d_big = d.to_biguint().unwrap();
        let (x_found, y_found) = chakravala(&d_big).expect("Solution");

        let (x_brute, y_brute) = brute_force_pell(d);
        assert_eq!(x_found, x_brute.to_biguint().unwrap());
        assert_eq!(y_found, y_brute.to_biguint().unwrap());
    }
}

// ---------------- Perfect-square guard ----------------
#[test]
fn perfect_square_returns_none() {
    for s in 1u64..=31 {
        let d = s * s;
        let d_big = d.to_biguint().unwrap();
        assert!(chakravala(&d_big).is_none(), "d = {} should yield None", d);
    }
}

// ---------------- Helper ----------------
/// Brute-force search for the fundamental solution of x² − d·y² = 1.
/// Strictly for small d (d ≤ 50) so a naive search is fine.
fn brute_force_pell(d: u64) -> (u64, u64) {
    for y in 1u64.. {
        let rhs = d * y * y + 1;
        let x = (rhs as f64).sqrt() as u64;
        if x * x == rhs { return (x, y); }
    }
    unreachable!("No solution found for d={}", d)
} 