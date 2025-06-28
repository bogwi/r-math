use criterion::{black_box, criterion_group, criterion_main, Criterion};
use num_bigint::ToBigUint;
use rma::chakravala::chakravala;

/// Benchmark Chakravala method for a variety of non-square discriminants.
fn bench_chakravala(c: &mut Criterion) {
    let mut group = c.benchmark_group("chakravala");

    // A representative set of non-square d values, increasing in size.
    let ds: [u64; 8] = [2, 13, 61, 67, 421, 991, 1021, 2027];

    for &d in &ds {
        group.bench_with_input(format!("d={}", d), &d, |b, &d_val| {
            // Pre-compute BigUint once per benchmark invocation.
            let d_big = d_val.to_biguint().unwrap();
            b.iter(|| {
                // Clone so each iteration works on a fresh copy.
                let result = chakravala(black_box(&d_big));
                black_box(result);
            });
        });
    }

    group.finish();
}

criterion_group!(benches, bench_chakravala);
criterion_main!(benches); 