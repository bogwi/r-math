//! Berlekamp's Root Finding Algorithm for polynomials over finite fields (prime fields).
//!
//! This implementation supports polynomials over GF(p) for prime p, using u128 for all arithmetic.
//!
//! # Example
//!
//! ```
//! use rma::berlekamp::{Poly, berlekamp_roots};
//! let p = 7u128;
//! let f = Poly::new(vec![5, 0, 1], p); // x^2 + 5
//! let mut roots = berlekamp_roots(&f);
//! roots.sort();
//! assert_eq!(roots, vec![3, 4]);
//! ```
//!
//! # References
//! - https://encyclopedia.pub/entry/32373
//! - https://en.wikipedia.org/wiki/Berlekamp%27s_algorithm
//! - https://www.usc.es/regaca/mladra/docs/berlekamp.pdf

use std::cmp::min;

/// Represents a polynomial over GF(p) as coefficients in ascending order (c_0 + c_1 x + ...)
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Poly {
    /// Coefficients in ascending order (c_0 + c_1 x + ...)
    pub coeffs: Vec<u128>,
    /// Modulus (prime)
    pub p: u128,
}

impl Poly {
    /// Create a new polynomial, reducing all coefficients mod p and trimming leading zeros.
    ///
    /// # Arguments
    /// * `coeffs` - Coefficients in ascending order (c_0 + c_1 x + ...)
    /// * `p` - Prime modulus
    ///
    /// # Example
    ///
    /// ```
    /// use rma::berlekamp::Poly;
    /// let poly = Poly::new(vec![1, 2, 3], 7);
    /// assert_eq!(poly.coeffs, vec![1, 2, 3]);
    /// ```
    pub fn new(mut coeffs: Vec<u128>, p: u128) -> Self {
        coeffs.iter_mut().for_each(|c| *c %= p);
        while coeffs.last().map_or(false, |c| *c == 0) && coeffs.len() > 1 {
            coeffs.pop();
        }
        Poly {
            coeffs,
            p,
        }
    }

    /// Degree of the polynomial (deg 0 for constant, -1 for zero poly)
    pub fn degree(&self) -> isize {
        let mut d = self.coeffs.len() as isize - 1;
        while d > 0 && self.coeffs[d as usize] == 0 {
            d -= 1;
        }
        if d == 0 && self.coeffs[0] == 0 {
            -1
        } else {
            d
        }
    }

    /// Evaluate the polynomial at x (mod p).
    ///
    /// # Arguments
    /// * `x` - The value to evaluate at
    ///
    /// # Returns
    /// The result of the evaluation in GF(p)
    ///
    /// # Example
    /// ```
    /// use rma::berlekamp::Poly;
    /// let poly = Poly::new(vec![1, 2, 3], 7);
    /// assert_eq!(poly.eval(2), (1 + 2*2 + 3*4) % 7);
    /// ```
    pub fn eval(&self, x: u128) -> u128 {
        let mut res = 0u128;
        let mut pow = 1u128;
        for c in self.coeffs.iter() {
            res = (res + (*c * pow) % self.p) % self.p;
            pow = (pow * x) % self.p;
        }
        res
    }

    /// Addition of two polynomials (mod p).
    ///
    /// # Arguments
    /// * `other` - The other polynomial
    ///
    /// # Returns
    /// The sum as a new polynomial
    ///
    /// # Example
    /// ```
    /// use rma::berlekamp::Poly;
    /// let p = 7u128;
    /// let a = Poly::new(vec![1, 2, 3], p); // 1 + 2x + 3x^2
    /// let b = Poly::new(vec![6, 5], p);    // 6 + 5x
    /// let c = a.add(&b);
    /// assert_eq!(c.coeffs, vec![0, 0, 3]); // (1+6)%7, (2+5)%7, 3
    /// ```
    pub fn add(&self, other: &Poly) -> Poly {
        let n = min(self.coeffs.len(), other.coeffs.len());
        let mut coeffs = Vec::with_capacity(self.coeffs.len().max(other.coeffs.len()));
        for i in 0..n {
            coeffs.push((self.coeffs[i] + other.coeffs[i]) % self.p);
        }
        coeffs.extend_from_slice(&self.coeffs[n..]);
        coeffs.extend_from_slice(&other.coeffs[n..]);
        Poly::new(coeffs, self.p)
    }

    /// Subtraction of two polynomials (mod p).
    ///
    /// # Arguments
    /// * `other` - The other polynomial
    ///
    /// # Returns
    /// The difference as a new polynomial
    ///
    /// # Example
    /// ```
    /// use rma::berlekamp::Poly;
    /// let p = 7u128;
    /// let a = Poly::new(vec![1, 2, 3], p); // 1 + 2x + 3x^2
    /// let b = Poly::new(vec![6, 5], p);    // 6 + 5x
    /// let c = a.sub(&b);
    /// assert_eq!(c.coeffs, vec![2, 4, 3]); // (1-6)%7=2, (2-5)%7=4, 3
    /// ```
    pub fn sub(&self, other: &Poly) -> Poly {
        let n = min(self.coeffs.len(), other.coeffs.len());
        let mut coeffs = Vec::with_capacity(self.coeffs.len().max(other.coeffs.len()));
        for i in 0..n {
            coeffs.push((self.coeffs[i] + self.p - other.coeffs[i]) % self.p);
        }
        coeffs.extend_from_slice(&self.coeffs[n..]);
        for i in n..other.coeffs.len() {
            coeffs.push((self.p - other.coeffs[i]) % self.p);
        }
        Poly::new(coeffs, self.p)
    }

    /// Scalar multiplication (mod p).
    ///
    /// # Arguments
    /// * `k` - Scalar value
    ///
    /// # Returns
    /// The scaled polynomial
    ///
    /// # Example
    /// ```
    /// use rma::berlekamp::Poly;
    /// let p = 7u128;
    /// let a = Poly::new(vec![1, 2, 3], p); // 1 + 2x + 3x^2
    /// let b = a.scale(3);
    /// assert_eq!(b.coeffs, vec![3, 6, 2]); // (1*3)%7=3, (2*3)%7=6, (3*3)%7=2
    /// ```
    pub fn scale(&self, k: u128) -> Poly {
        let coeffs = self.coeffs.iter().map(|c| (*c * k) % self.p).collect();
        Poly::new(coeffs, self.p)
    }

    /// Multiplication of two polynomials (mod p).
    ///
    /// Uses Karatsuba multiplication for large degree polynomials.
    ///
    /// # Arguments
    /// * `other` - The other polynomial
    ///
    /// # Returns
    /// The product as a new polynomial
    ///
    /// # Example
    /// ```
    /// use rma::berlekamp::Poly;
    /// let a = Poly::new(vec![1, 2], 7); // 1 + 2x
    /// let b = Poly::new(vec![3, 4], 7); // 3 + 4x
    /// let c = a.mul(&b); // (1 + 2x) * (3 + 4x) = 3 + 4x + 6x + 8x^2 = 3 + 10x + 8x^2
    /// assert_eq!(c.coeffs, vec![3, 3, 1]); // mod 7: 3 + 3x + 1x^2
    /// ```
    pub fn mul(&self, other: &Poly) -> Poly {
        const KARATSUBA_THRESHOLD: isize = 4;
        if self.degree() < 0 || other.degree() < 0 {
            return Poly::new(vec![], self.p);
        }
        if self.degree() < KARATSUBA_THRESHOLD || other.degree() < KARATSUBA_THRESHOLD {
            return self.naive_mul(other);
        }
        self.karatsuba_mul(other)
    }

    /// Naive polynomial multiplication (mod p).
    ///
    /// Used internally for small degree polynomials.
    fn naive_mul(&self, other: &Poly) -> Poly {
        if self.coeffs.is_empty() || other.coeffs.is_empty() {
            return Poly::new(vec![], self.p);
        }
        let mut coeffs = vec![0u128; self.coeffs.len() + other.coeffs.len() - 1];
        if coeffs.is_empty() && (self.coeffs.len() == 1 && other.coeffs.len() == 1) {
            coeffs.push(0u128);
        }
        for (i, a) in self.coeffs.iter().enumerate() {
            if *a == 0 { continue; }
            for (j, b) in other.coeffs.iter().enumerate() {
                if *b == 0 { continue; }
                let term = (*a * *b) % self.p;
                coeffs[i + j] = (coeffs[i + j] + term) % self.p;
            }
        }
        Poly::new(coeffs, self.p)
    }

    /// Karatsuba multiplication for polynomials (mod p).
    ///
    /// Used internally for large degree polynomials.
    fn karatsuba_mul(&self, other: &Poly) -> Poly {
        let p = self.p;
        let n1 = self.coeffs.len();
        let n2 = other.coeffs.len();
        if n1 == 0 || n2 == 0 {
            return Poly::new(vec![], p);
        }
        if n1 == 1 { return other.scale(self.coeffs[0]); }
        if n2 == 1 { return self.scale(other.coeffs[0]); }
        let m = std::cmp::max(n1, n2) / 2;
        if m == 0 {
            return self.naive_mul(other);
        }
        let (a1_coeffs_slice, a0_coeffs_slice) = if n1 > m { self.coeffs.split_at(m) } else { (&[][..], &self.coeffs[..]) };
        let (b1_coeffs_slice, b0_coeffs_slice) = if n2 > m { other.coeffs.split_at(m) } else { (&[][..], &other.coeffs[..]) };
        let a0 = Poly::new(a0_coeffs_slice.to_vec(), p);
        let a1 = if a1_coeffs_slice.is_empty() { Poly::new(vec![], p) } else { Poly::new(a1_coeffs_slice.to_vec(), p) };
        let b0 = Poly::new(b0_coeffs_slice.to_vec(), p);
        let b1 = if b1_coeffs_slice.is_empty() { Poly::new(vec![], p) } else { Poly::new(b1_coeffs_slice.to_vec(), p) };
        let z0 = a0.mul(&b0);
        let z2 = a1.mul(&b1);
        let a0_plus_a1 = a0.add(&a1);
        let b0_plus_b1 = b0.add(&b1);
        let z1_full = a0_plus_a1.mul(&b0_plus_b1);
        let z1 = z1_full.sub(&z0).sub(&z2);
        let res_len = if self.degree() < 0 || other.degree() < 0 {
            0
        } else {
            n1 + n2 - 1
        };
        if res_len == 0 {
            return Poly::new(vec![], p);
        }
        let mut res_coeffs = vec![0u128; res_len];
        for (i, c) in z0.coeffs.iter().enumerate() {
            if i < res_coeffs.len() {
                res_coeffs[i] = (res_coeffs[i] + *c) % p;
            }
        }
        for (i, c) in z1.coeffs.iter().enumerate() {
            let target_idx = i + m;
            if target_idx < res_coeffs.len() {
                res_coeffs[target_idx] = (res_coeffs[target_idx] + *c) % p;
            }
        }
        for (i, c) in z2.coeffs.iter().enumerate() {
            let target_idx = i + 2 * m;
            if target_idx < res_coeffs.len() {
                res_coeffs[target_idx] = (res_coeffs[target_idx] + *c) % p;
            }
        }
        Poly::new(res_coeffs, p)
    }

    /// Polynomial division with remainder: self = q * d + r
    ///
    /// # Arguments
    /// * `d` - The divisor polynomial
    ///
    /// # Returns
    /// A tuple (q, r) where q is the quotient and r is the remainder
    pub fn div_rem(&self, d: &Poly) -> (Poly, Poly) {
        let n = self.coeffs.len();
        let m = d.coeffs.len();
        let p = self.p;
        if d.degree() < 0 {
            panic!("Division by zero polynomial");
        }
        if n < m {
            return (Poly::new(vec![0u128], p), self.clone());
        }
        let d_lead = d.coeffs[m - 1];
        let d_lead_inv = modinv(d_lead, p).expect("No inverse for leading coeff");
        let mut rem = self.coeffs.clone();
        let mut quot = vec![0u128; n - m + 1];
        for _ in (m - 1..n).rev() {
            let rem_deg = rem.len() - 1;
            let lead = rem[rem_deg];
            if lead == 0 {
                rem.pop();
                continue;
            }
            let q_coeff = (lead * d_lead_inv) % p;
            quot[rem_deg - (m - 1)] = q_coeff;
            for j in 0..m {
                let idx = rem_deg - (m - 1) + j;
                if idx < rem.len() {
                    let sub = (d.coeffs[j] * q_coeff) % p;
                    rem[idx] = (rem[idx] + p - sub) % p;
                }
            }
            rem.pop();
        }
        while rem.last().map_or(false, |c| *c == 0) && rem.len() > 1 {
            rem.pop();
        }
        (Poly::new(quot, p), Poly::new(rem, p))
    }

    /// Compute the greatest common divisor (GCD) of two polynomials (mod p).
    ///
    /// # Arguments
    /// * `a` - First polynomial
    /// * `b` - Second polynomial
    ///
    /// # Returns
    /// The monic GCD polynomial
    pub fn gcd(a: &Poly, b: &Poly) -> Poly {
        let mut a = a.clone();
        let mut b = b.clone();
        while b.degree() >= 0 {
            let (_, r) = a.div_rem(&b);
            a = b;
            b = r;
        }
        let lead = a.coeffs[a.degree() as usize];
        let inv = modinv(lead, a.p).unwrap();
        a.scale(inv)
    }

    /// Polynomial remainder: self mod d
    ///
    /// # Arguments
    /// * `d` - The divisor polynomial
    ///
    /// # Returns
    /// The remainder polynomial
    pub fn rem(&self, d: &Poly) -> Poly {
        let (_, r) = self.div_rem(d);
        r
    }

    /// Compute x^k mod f(x)
    ///
    /// # Arguments
    /// * `k_orig` - Exponent
    /// * `f` - Modulus polynomial
    ///
    /// # Returns
    /// The polynomial x^k mod f(x)
    pub fn x_pow_mod(k_orig: u128, f: &Poly) -> Poly {
        let mut res = Poly::new(vec![1u128], f.p);
        let mut x = Poly::new(vec![0u128, 1u128], f.p);
        let mut k = k_orig;
        if k == 0 {
            return Poly::new(vec![1u128], f.p);
        }
        while k > 0 {
            if k % 2 == 1 {
                res = res.mul(&x).rem(f);
            }
            x = x.mul(&x).rem(f);
            k /= 2;
        }
        res
    }

    /// Returns the Berlekamp matrix Q for f(x) over GF(p).
    ///
    /// # Arguments
    /// * `f` - The polynomial
    ///
    /// # Returns
    /// The Berlekamp matrix as a vector of vectors
    pub fn berlekamp_matrix(f: &Poly) -> Vec<Vec<u128>> {
        let deg = f.degree();
        if deg < 1 {
            return vec![];
        }
        let n = deg as usize;
        let p = f.p;
        let mut q_matrix = vec![vec![0u128; n]; n];
        for i_col in 0..n {
            let exponent = (i_col as u128) * p;
            let xp_poly = Poly::x_pow_mod(exponent, f);
            for j_row in 0..n {
                q_matrix[j_row][i_col] = if j_row < xp_poly.coeffs.len() {
                    xp_poly.coeffs[j_row]
                } else {
                    0u128
                };
            }
        }
        q_matrix
    }

    /// Compute the formal derivative of the polynomial.
    ///
    /// # Returns
    /// The derivative polynomial
    ///
    /// # Example
    /// ```
    /// use rma::berlekamp::Poly;
    /// let f = Poly::new(vec![3, 2, 0, 1], 5);
    /// let df = f.derivative();
    /// assert_eq!(df.coeffs, vec![2, 0, 3]);
    /// ```
    pub fn derivative(&self) -> Poly {
        if self.coeffs.len() <= 1 {
            return Poly::new(vec![0u128], self.p);
        }
        let mut deriv_coeffs = Vec::with_capacity(self.coeffs.len() - 1);
        for (i, coeff) in self.coeffs.iter().enumerate().skip(1) {
            deriv_coeffs.push((*coeff * (i as u128)) % self.p);
        }
        Poly::new(deriv_coeffs, self.p)
    }
}

/// Modular inverse using extended Euclidean algorithm
///
/// # Arguments
/// * `a_orig` - Value to invert
/// * `m` - Modulus
///
/// # Returns
/// Some(inverse) if exists, None otherwise
fn modinv(a_orig: u128, m: u128) -> Option<u128> {
    if m == 0 || m == 1 {
        return None;
    }
    let mut r_prev = m;
    let mut r_curr = a_orig % m;
    if r_curr == 0 {
        return None;
    }
    let mut t_prev = 0u128;
    let mut t_curr = 1u128;
    while r_curr != 0 {
        let q = r_prev / r_curr;
        let r_next = r_prev % r_curr;
        r_prev = r_curr;
        r_curr = r_next;
        let q_mul_t_curr = q * t_curr;
        let val_q_mul_t_curr_mod_m = q_mul_t_curr % m;
        let t_next = if t_prev >= val_q_mul_t_curr_mod_m {
            (t_prev - val_q_mul_t_curr_mod_m) % m
        } else {
            (m + t_prev - val_q_mul_t_curr_mod_m) % m
        };
        t_prev = t_curr;
        t_curr = t_next;
    }
    if r_prev == 1 {
        Some(t_prev % m)
    } else {
        None
    }
}

/// Helper function to factorize a monic, square-free polynomial using Berlekamp's matrix method.
///
/// # Arguments
/// * `f_sf_monic` - Monic, square-free polynomial
/// * `p` - Prime modulus
///
/// # Returns
/// Vector of roots in GF(p)
fn berlekamp_factorize_squarefree(f_sf_monic: &Poly, p: u128) -> Vec<u128> {
    let mut roots = Vec::new();
    let mut factors_to_process = vec![f_sf_monic.clone()];
    while let Some(mut poly_to_factor) = factors_to_process.pop() {
        let deg = poly_to_factor.degree();
        if deg < 1 {
            continue;
        }
        if deg == 1 {
            let b = poly_to_factor.coeffs[0];
            let root = (p - b % p) % p;
            roots.push(root);
            continue;
        }
        let lead_coeff_factor = poly_to_factor.coeffs[poly_to_factor.degree() as usize];
        if lead_coeff_factor != 1 {
            if let Some(inv_lead_factor) = modinv(lead_coeff_factor, p) {
                poly_to_factor = poly_to_factor.scale(inv_lead_factor);
            } else {
                continue;
            }
        }
        let current_deg_val = poly_to_factor.degree();
        if current_deg_val < 1 {
            continue;
        }
        if current_deg_val == 1 {
            let b = poly_to_factor.coeffs[0];
            let root = (p - b % p) % p;
            roots.push(root);
            continue;
        }
        let q_matrix = Poly::berlekamp_matrix(&poly_to_factor);
        if q_matrix.is_empty() || q_matrix[0].is_empty() {
            continue;
        }
        let null_space = nullspace_modp(&q_matrix, p);
        if null_space.len() <= 1 {
            continue;
        }
        let mut split_occurred = false;
        for v_basis_vec in null_space.iter() {
            let v_poly = Poly::new(v_basis_vec.clone(), p);
            if v_poly.degree() < 1 {
                continue;
            }
            let s_iteration_heuristic_limit = (current_deg_val as u64).max(2) * 2;
            let s_loop_upper_bound_u64: u64 = if p < 150 { p as u64 - 1 } else { s_iteration_heuristic_limit.min(p as u64 - 1) };
            for s_val_u64 in 0..=s_loop_upper_bound_u64 {
                let s_val = s_val_u64 as u128;
                if p >= 150 && s_val >= p {
                    break;
                }
                let h_poly = v_poly.sub(&Poly::new(vec![s_val], p));
                if h_poly.degree() < 0 {
                    continue;
                }
                if h_poly.degree() == 0 && h_poly.coeffs[0] != 0 {
                    continue;
                }
                let gcd_poly = Poly::gcd(&poly_to_factor, &h_poly);
                let gcd_deg = gcd_poly.degree();
                if gcd_deg >= 1 && gcd_deg < current_deg_val {
                    factors_to_process.push(gcd_poly.clone());
                    let (quotient_poly, _remainder_poly_split) = poly_to_factor.div_rem(&gcd_poly);
                    if quotient_poly.degree() >= 1 {
                        factors_to_process.push(quotient_poly);
                    }
                    split_occurred = true;
                    break;
                }
            }
            if split_occurred {
                break;
            }
        }
    }
    roots.sort();
    roots.dedup();
    roots
}

/// Compute the nullspace of (Q - I) mod p, where Q is n x n, returns list of basis vectors
///
/// # Arguments
/// * `q` - Matrix Q
/// * `p` - Prime modulus
///
/// # Returns
/// List of basis vectors for the nullspace
fn nullspace_modp(q: &Vec<Vec<u128>>, p: u128) -> Vec<Vec<u128>> {
    let n = q.len();
    if n == 0 {
        return vec![vec![]];
    }
    let mut mat = q.clone();
    for i in 0..n {
        if mat[i][i] < 1 {
            mat[i][i] += p;
            mat[i][i] -= 1;
        } else {
            mat[i][i] -= 1;
        }
    }
    let mut basis = vec![];
    let mut used_pivot_rows = vec![false; n];
    for col in 0..n {
        let mut pivot_row_opt = None;
        for r_idx in 0..n {
            if !used_pivot_rows[r_idx] && mat[r_idx][col] != 0 {
                pivot_row_opt = Some(r_idx);
                break;
            }
        }
        if let Some(pivot_row) = pivot_row_opt {
            used_pivot_rows[pivot_row] = true;
            let inv = modinv(mat[pivot_row][col], p).expect("Modular inverse failed for a pivot element");
            for j in 0..n {
                mat[pivot_row][j] = (mat[pivot_row][j] * inv) % p;
            }
            for i in 0..n {
                if i != pivot_row {
                    let factor = mat[i][col];
                    if factor != 0 {
                        for j in 0..n {
                            let val_to_sub = (mat[pivot_row][j] * factor) % p;
                            if mat[i][j] < val_to_sub {
                                mat[i][j] += p;
                                mat[i][j] -= val_to_sub;
                            } else {
                                mat[i][j] -= val_to_sub;
                            }
                        }
                    }
                }
            }
        }
    }
    let mut leading_col_for_row = vec![None; n];
    let mut is_pivot_col = vec![false; n];
    for r_idx in 0..n {
        for c_idx in 0..n {
            if mat[r_idx][c_idx] == 1 {
                let mut is_leading = true;
                for k_idx in 0..c_idx {
                    if mat[r_idx][k_idx] != 0 {
                        is_leading = false;
                        break;
                    }
                }
                if is_leading {
                    leading_col_for_row[r_idx] = Some(c_idx);
                    is_pivot_col[c_idx] = true;
                    break;
                }
            } else if mat[r_idx][c_idx] != 0 {
                break;
            }
        }
    }
    let mut free_cols = vec![];
    for c_idx in 0..n {
        if !is_pivot_col[c_idx] {
            free_cols.push(c_idx);
        }
    }
    for &free_col_idx in &free_cols {
        let mut v = vec![0u128; n];
        v[free_col_idx] = 1u128;
        for r_idx in 0..n {
            if let Some(pivot_col_idx) = leading_col_for_row[r_idx] {
                if pivot_col_idx != free_col_idx {
                    let val = mat[r_idx][free_col_idx];
                    if val != 0 {
                        v[pivot_col_idx] = p - val;
                    } else {
                        v[pivot_col_idx] = 0u128;
                    }
                }
            }
        }
        basis.push(v);
    }
    if basis.is_empty() && n > 0 {
        basis.push(vec![0u128; n]);
    }
    basis
}

/// Find all roots of f(x) in GF(p) using Berlekamp's algorithm.
///
/// # Arguments
/// * `original_f` - The polynomial to find roots of
///
/// # Returns
/// Vector of roots in GF(p)
///
/// # Example
/// ```
/// use rma::berlekamp::{Poly, berlekamp_roots};
/// let p = 7u128;
/// let f = Poly::new(vec![5, 0, 1], p); // x^2 + 5
/// let mut roots = berlekamp_roots(&f);
/// roots.sort();
/// assert_eq!(roots, vec![3, 4]);
/// ```
pub fn berlekamp_roots(original_f: &Poly) -> Vec<u128> {
    let p = original_f.p;
    let n = original_f.degree();
    if n < 0 {
        return vec![];
    }
    if n == 0 {
        return vec![];
    }
    if n == 1 {
        let a = original_f.coeffs[1];
        let b = original_f.coeffs[0];
        if a == 0 {
            return if b == 0 { vec![0u128] } else { vec![] };
        }
        let inv_a = modinv(a, p).expect("Inverse for linear term failed");
        let root = ((p - b % p) % p * inv_a) % p;
        return vec![root];
    }
    let f_derivative = original_f.derivative();
    if f_derivative.degree() < 0 {
        if p == 0 || p == 1 {
            return vec![];
        }
        let mut h_coeffs_vec: Vec<u128> = Vec::new();
        let mut max_h_deg_plus_1 = 0;
        let mut is_valid_form = true;
        for (i_idx, coeff) in original_f.coeffs.iter().enumerate() {
            if *coeff != 0 {
                if (i_idx as u128) % p != 0 {
                    is_valid_form = false;
                    break;
                }
                let h_idx = (i_idx as u128 / p) as usize;
                if h_idx + 1 > max_h_deg_plus_1 {
                    max_h_deg_plus_1 = h_idx + 1;
                    if h_coeffs_vec.len() < max_h_deg_plus_1 {
                        h_coeffs_vec.resize(max_h_deg_plus_1, 0u128);
                    }
                }
                h_coeffs_vec[h_idx] = *coeff;
            }
        }
        if !is_valid_form {
            return vec![];
        }
        while h_coeffs_vec.last().map_or(false, |c| *c == 0) && h_coeffs_vec.len() > 1 {
            h_coeffs_vec.pop();
        }
        if h_coeffs_vec.is_empty() {
            h_coeffs_vec.push(0u128);
        }
        let h_poly = Poly::new(h_coeffs_vec, p);
        return berlekamp_roots(&h_poly);
    }
    let mut f_monic = original_f.clone();
    let lead_coeff_orig = f_monic.coeffs[f_monic.degree() as usize];
    if lead_coeff_orig != 1 {
        if let Some(inv_lead) = modinv(lead_coeff_orig, p) {
            f_monic = f_monic.scale(inv_lead);
        } else {
            return vec![];
        }
    }
    let f_monic_derivative = f_monic.derivative();
    let common_divisor = Poly::gcd(&f_monic, &f_monic_derivative);
    let mut all_roots = Vec::new();
    let (square_free_part, _) = f_monic.div_rem(&common_divisor);
    if square_free_part.degree() >= 1 {
        let roots_from_sf_part = berlekamp_factorize_squarefree(&square_free_part, p);
        all_roots.extend(roots_from_sf_part);
    }
    if common_divisor.degree() >= 1 {
        let roots_from_gcd_part = berlekamp_roots(&common_divisor);
        all_roots.extend(roots_from_gcd_part);
    }
    all_roots.sort();
    all_roots.dedup();
    all_roots
}

#[cfg(test)]
mod tests {
    use super::*;
    fn assert_roots_match(found: &mut Vec<u128>, expected: Vec<u128>) {
        found.sort();
        let mut expected_sorted = expected;
        expected_sorted.sort();
        assert_eq!(*found, expected_sorted);
    }
    #[test]
    fn test_poly_add_mul() {
        let p = 7u128;
        let a = Poly::new(vec![1, 2, 3], p);
        let b = Poly::new(vec![6, 5], p);
        let c = a.add(&b);
        assert_eq!(c.coeffs, vec![0, 0, 3]);
        let d = a.mul(&b);
        assert_eq!(d.coeffs.len(), 4);
    }
    #[test]
    fn test_poly_gcd() {
        let p = 5u128;
        let a = Poly::new(vec![1, 0, 1], p);
        let b = Poly::new(vec![2, 0, 1], p);
        let g = Poly::gcd(&a, &b);
        assert!(g.degree() == 0);
        assert_eq!(g.coeffs, vec![1]);
    }
    #[test]
    fn test_linear_root() {
        let p = 7u128;
        let f = Poly::new(vec![2, 3], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![4]);
    }
    #[test]
    fn test_quadratic_roots() {
        let p = 7u128;
        let f = Poly::new(vec![5, 0, 1], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![3, 4]);
    }
    #[test]
    fn test_cubic_roots() {
        let p = 5u128;
        let f = Poly::new(vec![0, 4, 0, 1], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![0, 1, 4]);
    }
    #[test]
    fn test_poly_derivative() {
        let p = 5u128;
        let f = Poly::new(vec![3, 2, 0, 1], p);
        let df = f.derivative();
        assert_eq!(df.coeffs, vec![2, 0, 3]);
        let c = Poly::new(vec![4], p);
        let dc = c.derivative();
        assert_eq!(dc.coeffs, vec![0]);
        assert_eq!(dc.degree(), -1);
    }
    #[test]
    fn test_irreducible_over_gf2() {
        let p = 2u128;
        let f = Poly::new(vec![1, 1, 0, 1], p);
        let roots = berlekamp_roots(&f);
        assert_eq!(roots, vec![]);
    }
    #[test]
    fn test_repeated_root() {
        let p = 2u128;
        let f = Poly::new(vec![1, 0, 1], p);
        let roots = berlekamp_roots(&f);
        assert_eq!(roots, vec![1]);
    }
    #[test]
    fn test_high_degree_all_roots() {
        let p = 2u128;
        let f = Poly::new(vec![1, 0, 0, 0, 0, 1], p);
        let roots = berlekamp_roots(&f);
        assert_eq!(roots, vec![1]);
    }
    #[test]
    fn test_irreducible_over_gf3_reducible_over_gf9() {
        let p = 3u128;
        let f = Poly::new(vec![1, 0, 1], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![]);
    }
    #[test]
    fn test_zero_polynomial() {
        let p = 5u128;
        let f = Poly::new(vec![0], p);
        let roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots.clone(), vec![]);
    }
    #[test]
    fn test_constant_nonzero_polynomial() {
        let p = 5u128;
        let f = Poly::new(vec![3], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![]);
    }
    #[test]
    fn test_gf2_x_squared_plus_x_plus_1_irreducible() {
        let p = 2u128;
        let f = Poly::new(vec![1, 1, 1], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![]);
    }
    #[test]
    fn test_gf2_x_squared_repeated_root_0() {
        let p = 2u128;
        let f = Poly::new(vec![0, 0, 1], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![0]);
    }
    #[test]
    fn test_gf2_x_cubed_plus_x_squared_plus_x_plus_1_repeated_root_1_mult_3() {
        let p = 2u128;
        let f = Poly::new(vec![1, 1, 1, 1], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![1]);
    }
    #[test]
    fn test_gf2_x4_plus_x_plus_1_irreducible() {
        let p = 2u128;
        let f = Poly::new(vec![1, 1, 0, 0, 1], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![]);
    }
    #[test]
    fn test_gf2_x4_plus_x3_plus_x2_plus_x_plus_1_irreducible() {
        let p = 2u128;
        let f = Poly::new(vec![1, 1, 1, 1, 1], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![]);
    }
    #[test]
    fn test_gf2_x5_plus_x2_plus_1_irreducible() {
        let p = 2u128;
        let f = Poly::new(vec![1, 0, 1, 0, 0, 1], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![]);
    }
    #[test]
    fn test_gf3_x_squared_plus_x_plus_2_irreducible() {
        let p = 3u128;
        let f = Poly::new(vec![2, 1, 1], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![]);
    }
    #[test]
    fn test_gf3_x_cubed_plus_2x_plus_1_irreducible() {
        let p = 3u128;
        let f = Poly::new(vec![1, 2, 0, 1], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![]);
    }
    #[test]
    fn test_gf3_x_cubed_plus_2x_all_roots() {
        let p = 3u128;
        let f = Poly::new(vec![0, 2, 0, 1], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![0, 1, 2]);
    }
    #[test]
    fn test_gf5_x_squared_plus_2_irreducible() {
        let p = 5u128;
        let f = Poly::new(vec![2, 0, 1], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![]);
    }
    #[test]
    fn test_gf5_x_squared_plus_3_irreducible() {
        let p = 5u128;
        let f = Poly::new(vec![3, 0, 1], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![]);
    }
    #[test]
    fn test_gf5_x4_plus_1_product_of_irreducibles_no_roots() {
        let p = 5u128;
        let f = Poly::new(vec![1, 0, 0, 0, 1], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![]);
    }
    #[test]
    fn test_gf5_x5_minus_x_all_roots() {
        let p = 5u128;
        let f = Poly::new(vec![0, 4, 0, 0, 0, 1], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![0, 1, 2, 3, 4]);
    }
    #[test]
    fn test_gf13_x_squared_minus_3() {
        let p = 13u128;
        let f = Poly::new(vec![10, 0, 1], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![4, 9]);
    }
    #[test]
    fn test_gf2_x_cubed_plus_1_reducible_square_free_not_all_linear() {
        let p = 2u128;
        let f = Poly::new(vec![1, 0, 0, 1], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![1]);
    }
    #[test]
    fn test_gf3_x_cubed_plus_x_reducible_square_free_not_all_linear() {
        let p = 3u128;
        let f = Poly::new(vec![0, 1, 0, 1], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![0]);
    }
    #[test]
    fn test_gf3_repeated_factors_x_minus_1_sq_x_minus_2() {
        let p = 3u128;
        let f = Poly::new(vec![1, 2, 2, 1], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![1, 2]);
    }
    #[test]
    fn test_gf2_x_cubed_repeated_root_0_mult_3() {
        let p = 2u128;
        let f = Poly::new(vec![0, 0, 0, 1], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![0]);
    }
    #[test]
    fn test_gf5_x_minus_1_cubed_x_minus_2_squared() {
        let p = 5u128;
        let f1 = Poly::new(vec![4, 3, 2, 1], p);
        let f2 = Poly::new(vec![4, 1, 1], p);
        let f = f1.mul(&f2);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![1, 2]);
    }
    #[test]
    fn test_gf2_f_prime_is_zero_x4_x2_1() {
        let p = 2u128;
        let f = Poly::new(vec![1, 0, 1, 0, 1], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![]);
    }
    #[test]
    fn test_gf3_f_prime_is_zero_x6_x3_1() {
        let p = 3u128;
        let f = Poly::new(vec![1, 0, 0, 1, 0, 0, 1], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![1]);
    }
    #[test]
    fn test_gf3_f_prime_is_zero_x_cubed_minus_1() {
        let p = 3u128;
        let f = Poly::new(vec![2, 0, 0, 1], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![1]);
    }
    #[test]
    fn test_gf5_non_monic_becomes_monic_irreducible() {
        let p = 5u128;
        let f = Poly::new(vec![1, 0, 2], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![]);
    }
    #[test]
    fn test_gf11_xp_minus_x_all_roots() {
        let p_val = 11u128;
        let p = p_val;
        let mut coeffs = vec![0u128; (p_val + 1) as usize];
        coeffs[1] = (p - 1) % p;
        coeffs[p_val as usize] = 1u128;
        let f = Poly::new(coeffs, p);
        let expected_roots: Vec<u128> = (0..p_val).collect();
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, expected_roots);
    }
    #[test]
    fn test_gf2_x5_plus_x_plus_1_reducible() {
        let p = 2u128;
        let f = Poly::new(vec![1, 1, 0, 0, 0, 1], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![]);
    }
    #[test]
    fn test_gf7_x3_plus_2x2_plus_x_plus_2_some_roots() {
        let p = 7u128;
        let f = Poly::new(vec![2, 1, 2, 1], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![5]);
    }
    #[test]
    fn test_another_from_wikipedia_x4_plus_1_gf3() {
        let p = 3u128;
        let f = Poly::new(vec![1, 0, 0, 0, 1], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![]);
    }
    #[test]
    fn test_poly_from_paper_x8_plus_x6_plus_x4_plus_x3_plus_1_gf2() {
        let p = 2u128;
        let f = Poly::new(vec![1, 0, 0, 1, 1, 0, 1, 0, 1], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![]);
    }
    #[test]
    fn test_poly_x_to_p_minus_x_plus_1_gf_p() {
        let p_val = 7u128;
        let p = p_val;
        let mut coeffs = vec![0u128; (p_val + 1) as usize];
        coeffs[0] = 1u128;
        coeffs[1] = (p - 1) % p;
        coeffs[p_val as usize] = 1u128;
        let f = Poly::new(coeffs, p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![]);
    }
    #[test]
    fn test_product_of_all_linear_factors_gf3() {
        let p = 3u128;
        let f = Poly::new(vec![0, 2, 0, 1], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![0, 1, 2]);
    }
    #[test]
    fn test_reducible_large_prime_gf101_x2_minus_4() {
        let p = 101u128;
        let f = Poly::new(vec![(p - 4) % p, 0, 1], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![2, 99]);
    }
    #[test]
    fn test_irreducible_large_prime_gf101_x2_plus_1() {
        let p = 101u128;
        let f = Poly::new(vec![1, 0, 1], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![10, 91]);
    }
    #[test]
    fn test_gf2_x6_plus_x5_plus_x4_plus_x3_plus_x2_plus_x_plus_1() {
        let p = 2u128;
        let f = Poly::new(vec![1, 1, 1, 1, 1, 1, 1], p);
        let mut roots = berlekamp_roots(&f);
        assert_roots_match(&mut roots, vec![]);
    }
}
