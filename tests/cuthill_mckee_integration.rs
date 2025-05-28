use rma::cuthill_mckee::{
    MatrixMarketMatrix, MtxType, cuthill_mckee_csr, read_matrix_market, write_matrix_market,
};
use std::fs;
use std::io::Write;
use tempfile::NamedTempFile;

fn csr_eq(a: &rma::cuthill_mckee::CsrMatrix, b: &rma::cuthill_mckee::CsrMatrix) -> bool {
    a.nrows == b.nrows && a.ncols == b.ncols && a.indptr == b.indptr && a.indices == b.indices
}

fn csc_eq(a: &rma::cuthill_mckee::CscMatrix, b: &rma::cuthill_mckee::CscMatrix) -> bool {
    a.nrows == b.nrows && a.ncols == b.ncols && a.indptr == b.indptr && a.indices == b.indices
}

#[test]
fn test_cuthill_mckee_serialization_roundtrip_csr() {
    let orig_matrix = read_matrix_market("test_data/medium_matrix.mtx", MtxType::Csr).unwrap();
    let out_path = "test_data/output_matrix.mtx";

    // 1. Write to disk
    write_matrix_market(out_path, &orig_matrix).unwrap();

    // 2. Read back
    let roundtrip_matrix = read_matrix_market(out_path, MtxType::Csr).unwrap();

    // 3. Compare matrices
    match (&orig_matrix, &roundtrip_matrix) {
        (MatrixMarketMatrix::Csr(csr1), MatrixMarketMatrix::Csr(csr2)) => {
            assert!(csr_eq(csr1, csr2), "CSR matrices differ after roundtrip");
            // 4. Run cuthill_mckee_csr on both
            let order1 = cuthill_mckee_csr(csr1).expect("Cuthill-McKee failed on original");
            let order2 = cuthill_mckee_csr(csr2).expect("Cuthill-McKee failed on roundtrip");
            assert_eq!(
                order1, order2,
                "Cuthill-McKee orders differ after roundtrip"
            );

            // 5. Check it's a valid permutation
            let mut seen = vec![false; csr1.nrows];
            order1.iter().for_each(|&i| {
                assert!(i < csr1.nrows);
                assert!(!seen[i]);
                seen[i] = true;
            });
            seen.fill(false);
            order2.iter().for_each(|&i| {
                assert!(i < csr2.nrows);
                assert!(!seen[i]);
                seen[i] = true;
            });
        }
        _ => panic!("Expected CSR matrices for both original and roundtrip"),
    }

    // 5. Clean up
    fs::remove_file(out_path).unwrap();
}

#[test]
fn test_cuthill_mckee_serialization_roundtrip_csc() {
    let orig_matrix = read_matrix_market("test_data/large_csc_matrix.mtx", MtxType::Csc).unwrap();
    let out_path = "test_data/output_csc_matrix.mtx";

    // 1. Write to disk
    write_matrix_market(out_path, &orig_matrix).unwrap();

    // 2. Read back
    let roundtrip_matrix = read_matrix_market(out_path, MtxType::Csc).unwrap();

    // 3. Compare matrices
    match (&orig_matrix, &roundtrip_matrix) {
        (MatrixMarketMatrix::Csc(csc1), MatrixMarketMatrix::Csc(csc2)) => {
            assert!(csc_eq(csc1, csc2), "CSC matrices differ after roundtrip");
            // 4. Run cuthill_mckee_csc on both
            let order1 = rma::cuthill_mckee::cuthill_mckee_csc(csc1)
                .expect("Cuthill-McKee failed on original");
            let order2 = rma::cuthill_mckee::cuthill_mckee_csc(csc2)
                .expect("Cuthill-McKee failed on roundtrip");
            assert_eq!(
                order1, order2,
                "Cuthill-McKee orders differ after roundtrip"
            );

            // 5. Check it's a valid permutation
            let mut seen = vec![false; csc1.nrows];
            order1.iter().for_each(|&i| {
                assert!(i < csc1.nrows);
                assert!(!seen[i]);
                seen[i] = true;
            });
            seen.fill(false);
            order2.iter().for_each(|&i| {
                assert!(i < csc2.nrows);
                assert!(!seen[i]);
                seen[i] = true;
            });
        }
        _ => panic!("Expected CSC matrices for both original and roundtrip"),
    }

    // 5. Clean up
    fs::remove_file(out_path).unwrap();
}

#[test]
fn test_read_matrix_market_corrupted_files() {
    // Helper to write a string to a temp file and return its path
    fn write_temp(contents: &str) -> String {
        let mut file = NamedTempFile::new().unwrap();
        write!(file, "{}", contents).unwrap();
        file.into_temp_path().to_str().unwrap().to_string()
    }

    let cases = vec![
        // Missing header
        "% Just a comment\n1 2 3\n",
        // Wrong number of columns in header
        "23 23\n1 1 10\n",
        // Out-of-bounds index
        "23 23 1\n24 1 10\n",
        // Non-integer value in header
        "a b c\n1 1 10\n",
        // Non-integer value in entry
        "23 23 1\n1 1 x\n",
        // Fewer entries than nnz
        "3 3 5\n1 1 1\n2 2 2\n3 3 3\n",
        // Entry with only one column
        "3 3 1\n1\n",
        // Negative or zero index
        "3 3 1\n0 1 1\n",
        // Empty file
        "",
    ];

    for (i, contents) in cases.iter().enumerate() {
        let path = write_temp(contents);
        let mut result = read_matrix_market(&path, MtxType::Csr);
        assert!(result.is_err(), "Corrupted file case {} should error", i);
        result = read_matrix_market(&path, MtxType::Csc);
        assert!(result.is_err(), "Corrupted file case {} should error", i);
        // tempfile cleans up automatically
    }
}
