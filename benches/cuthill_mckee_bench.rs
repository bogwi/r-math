use criterion::{Criterion, black_box, criterion_group, criterion_main};
use rand::prelude::*;
use rma::cuthill_mckee;

/// Generate a random sparse symmetric adjacency list with n nodes and avg_degree average degree.
fn random_sparse_symmetric_graph(n: usize, avg_degree: usize, seed: u64) -> Vec<Vec<usize>> {
    let mut adj = vec![Vec::new(); n];
    let mut rng = StdRng::seed_from_u64(seed);
    (0..n).into_iter().for_each(|i| {
        let mut neighbors = std::collections::HashSet::new();
        while neighbors.len() < avg_degree {
            let j = rng.gen_range(0..n);
            if j != i {
                neighbors.insert(j);
            }
        }
        neighbors.into_iter().for_each(|j| {
            adj[i].push(j);
            adj[j].push(i); // ensure symmetry
        });
    });
    // Remove duplicates
    adj.iter_mut().for_each(|nbrs| {
        nbrs.sort_unstable();
        nbrs.dedup();
    });
    adj
}

/// Generate a banded matrix adjacency list (bandwidth b)
fn banded_matrix_graph(n: usize, b: usize) -> Vec<Vec<usize>> {
    let mut adj = vec![Vec::new(); n];
    (0..n).into_iter().for_each(|i| {
        (i.saturating_sub(b)..=(i + b).min(n - 1))
            .into_iter()
            .for_each(|j| {
                if i != j {
                    adj[i].push(j);
                }
            });
    });
    adj
}

fn bench_cuthill_mckee(c: &mut Criterion) {
    // Benchmark on a random sparse symmetric graph
    let n = 100_000;
    let avg_degree = 10;
    let adj = random_sparse_symmetric_graph(n, avg_degree, 42);
    c.bench_function(&format!("cuthill_mckee_random_sparse_{}k", n / 1000), |b| {
        b.iter(|| {
            let order = cuthill_mckee(black_box(&adj)).unwrap();
            black_box(order);
        })
    });

    // Benchmark on a banded matrix graph
    let n = 100_000;
    let bandwidth = 5;
    let adj = banded_matrix_graph(n, bandwidth);
    c.bench_function(&format!("cuthill_mckee_banded_{}k_bw5", n / 1000), |b| {
        b.iter(|| {
            let order = cuthill_mckee(black_box(&adj)).unwrap();
            black_box(order);
        })
    });
}

criterion_group!(benches, bench_cuthill_mckee);
criterion_main!(benches);

// -----------------------------------------------------------------------------
// Performance Guidance:
//
// On a modern desktop CPU (e.g., Apple M1/M2, AMD Ryzen 5000+, Intel 12th+ gen):
//
// - For a random sparse symmetric graph with 10,000 nodes and avg_degree=10:
//     * Cutting-edge: < 100 ms per run (single-threaded), < 30 ms (multi-threaded)
//     * Good: < 200 ms per run
//
// - For a banded matrix graph with 10,000 nodes, bandwidth=5:
//     * Cutting-edge: < 50 ms per run (single-threaded), < 15 ms (multi-threaded)
//     * Good: < 100 ms per run
//
// These are rough targets; actual results depend on hardware, OS, and implementation details.
// -----------------------------------------------------------------------------
