// NOTE: Requires num-bigint with the `rand` feature enabled in Cargo.toml:
// num-bigint = { version = "0.4", features = ["rand"] }

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rand::SeedableRng;
use rand::rngs::StdRng;
use rand::prelude::SliceRandom;
use rand::Rng;
use rma::berlekamp::{Poly, berlekamp_roots};

/// Generate a random dense polynomial of degree `deg` over GF(p)
fn random_dense_poly(deg: usize, p: u128, rng: &mut StdRng) -> Poly {
    let mut coeffs: Vec<u128> = (0..deg).map(|_| rng.gen_range(0..p)).collect();
    // Ensure leading coefficient is nonzero
    let mut lead = 0u128;
    while lead == 0 {
        lead = rng.gen_range(0..p);
    }
    coeffs.push(lead);
    Poly::new(coeffs, p)
}

/// Generate a random sparse polynomial of degree `deg` over GF(p), with `nonzero` nonzero terms
fn random_sparse_poly(deg: usize, nonzero: usize, p: u128, rng: &mut StdRng) -> Poly {
    let mut coeffs = vec![0u128; deg + 1];
    let mut idxs: Vec<usize> = (0..deg).collect();
    idxs.shuffle(rng);
    for &i in idxs.iter().take(nonzero.saturating_sub(1)) {
        coeffs[i] = rng.gen_range(0..p);
    }
    // Ensure leading coefficient is nonzero
    let mut lead = 0u128;
    while lead == 0 {
        lead = rng.gen_range(0..p);
    }
    coeffs[deg] = lead;
    Poly::new(coeffs, p)
}

fn fixed_prime(bits: usize) -> u128 {
    match bits {
        // 8 => 251u128,
        // 16 => 65521u128,
        // 32 => 4_294_967_291u128,
        64 => 18_446_744_073_709_551_557u128,
        _ => 17u128,
    }
}

fn bench_berlekamp_roots(c: &mut Criterion) {
    let mut rng = StdRng::seed_from_u64(42);
    // Only the 256 degree polynomial is tested over the largest 64 bit prime.
    // This gives you an idea of the performance of the algorithm.
    let degrees = [256];
    let prime_bits = [64];
    for &deg in &degrees {
        for &bits in &prime_bits {
            let p = fixed_prime(bits);
            // Dense poly
            let poly = random_dense_poly(deg, p, &mut rng);
            let bench_name = format!("berlekamp_roots_dense_deg{}_p{}b", deg, bits);
            c.bench_function(&bench_name, |b| {
                b.iter(|| {
                    let roots = berlekamp_roots(black_box(&poly));
                    black_box(roots);
                })
            });
            // Sparse poly (10% nonzero terms, at least 2)
            let nonzero = (deg / 10).max(2);
            let poly_sparse = random_sparse_poly(deg, nonzero, p, &mut rng);
            let bench_name_sparse = format!("berlekamp_roots_sparse_deg{}_p{}b", deg, bits);
            c.bench_function(&bench_name_sparse, |b| {
                b.iter(|| {
                    let roots = berlekamp_roots(black_box(&poly_sparse));
                    black_box(roots);
                })
            });
        }
    }
}

criterion_group!(benches, bench_berlekamp_roots);
criterion_main!(benches); 