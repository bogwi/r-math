use criterion::{Criterion, black_box, criterion_group, criterion_main};
#[cfg(feature = "simd")]
use rma::freivalds_verify_simd;
use rma::{freivalds_verify_scalar, freivalds_verify_scalar_mod};

fn generate_matrix(n: usize) -> Vec<Vec<u64>> {
    (0..n)
        .map(|i| (0..n).map(|j| ((i * j + 1) % 100) as u64).collect())
        .collect()
}

fn bench_freivalds_scalar(c: &mut Criterion) {
    let n = 1024;
    let k = 10;
    let a = generate_matrix(n);
    let b = generate_matrix(n);
    let mut c_mat = vec![vec![0u64; n]; n];
    c_mat.iter_mut().enumerate().for_each(|(i, row)| {
        row.iter_mut().enumerate().for_each(|(j, cell)| {
            *cell = (0..n).fold(0, |acc, l| acc.wrapping_add(a[i][l].wrapping_mul(b[l][j])));
        });
    });
    c.bench_function(
        &format!("freivalds_verify scalar: {n} x {n} k={k}"),
        |bch| {
            bch.iter(|| {
                let result =
                    freivalds_verify_scalar(black_box(&a), black_box(&b), black_box(&c_mat), k);
                black_box(result)
            })
        },
    );
}

fn bench_freivalds_scalar_mod(c: &mut Criterion) {
    let n = 1024;
    let k = 10;
    let a = generate_matrix(n);
    let b = generate_matrix(n);
    let modulus = 1_000_000_007;
    let mut c_mat = vec![vec![0u64; n]; n];
    c_mat.iter_mut().enumerate().for_each(|(i, row)| {
        row.iter_mut().enumerate().for_each(|(j, cell)| {
            *cell = (0..n).fold(0, |acc, l| acc.wrapping_add(a[i][l].wrapping_mul(b[l][j])));
        });
    });
    c.bench_function(
        &format!("freivalds_verify scalar mod: {n} x {n} k={k} mod={modulus}"),
        |bch| {
            bch.iter(|| {
                let result = freivalds_verify_scalar_mod(
                    black_box(&a),
                    black_box(&b),
                    black_box(&c_mat),
                    k,
                    modulus,
                );
                black_box(result)
            })
        },
    );
}

fn custom_criterion() -> Criterion {
    Criterion::default().sample_size(50)
}

#[cfg(feature = "simd")]
fn bench_freivalds_simd(c: &mut Criterion) {
    let n = 1024;
    let k = 10;
    let a = generate_matrix(n);
    let b = generate_matrix(n);
    // Compute c = a * b
    let mut c_mat = vec![vec![0u64; n]; n];
    c_mat.iter_mut().enumerate().for_each(|(i, row)| {
        row.iter_mut().enumerate().for_each(|(j, cell)| {
            *cell = (0..n).fold(0, |acc, l| acc.wrapping_add(a[i][l].wrapping_mul(b[l][j])));
        });
    });
    c.bench_function(&format!("freivalds_verify simd: {n} x {n} k={k}"), |bch| {
        bch.iter(|| {
            let result = freivalds_verify_simd::<u64, 32>(
                black_box(&a),
                black_box(&b),
                black_box(&c_mat),
                k,
            );
            black_box(result)
        })
    });
}

#[cfg(feature = "simd")]
criterion_group!(
    benches,
    bench_freivalds_scalar,
    bench_freivalds_scalar_mod,
    bench_freivalds_simd
);
#[cfg(not(feature = "simd"))]
criterion_group! {
    name = benches;
    config = custom_criterion();
    targets = bench_freivalds_scalar, bench_freivalds_scalar_mod
}
criterion_main!(benches);
