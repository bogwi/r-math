//! Cuthill–McKee and Reverse Cuthill–McKee (RCM) algorithm for bandwidth reduction of sparse symmetric matrices.
//!
//! # Overview
//! This module implements the Cuthill–McKee (CM) and Reverse Cuthill–McKee (RCM) algorithms for reordering
//! sparse symmetric matrices to reduce their bandwidth. RCM is particularly effective for minimizing
//! fill-in during matrix factorization.
//!
//! # Features
//! - Support for adjacency lists, CSR, and CSC matrix formats
//! - Conversion from `sprs` sparse matrix types
//! - Optimized for both bandwidth reduction and fill-in minimization
//!
//! # Example
//! ```
//! use rma::cuthill_mckee::{CscMatrix, reverse_cuthill_mckee_csc};
//!
//! // Create a sparse matrix in CSC format
//! let matrix = CscMatrix {
//!     nrows: 4,
//!     ncols: 4,
//!     indptr: vec![0, 2, 4, 6, 8],
//!     indices: vec![0, 1, 0, 2, 1, 3, 2, 3],
//!     data: Some(vec![1.0; 8]),
//! };
//!
//! // Apply RCM reordering
//! let perm = reverse_cuthill_mckee_csc(&matrix).unwrap();
//! ```

use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CuthillMckeeError {
    EmptyAdjacencyList,
    NeighborIndexOutOfBounds {
        node: usize,
        neighbor: usize,
        len: usize,
    },
    CsrIndptrLengthMismatch {
        expected: usize,
        actual: usize,
    },
    CsrIndexOutOfBounds {
        index: usize,
        ncols: usize,
    },
    CscIndptrLengthMismatch {
        expected: usize,
        actual: usize,
    },
    CscIndexOutOfBounds {
        index: usize,
        nrows: usize,
    },
    SprsMatrixNotCsc,
}

impl fmt::Display for CuthillMckeeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CuthillMckeeError::EmptyAdjacencyList => write!(f, "Adjacency list is empty"),
            CuthillMckeeError::NeighborIndexOutOfBounds {
                node,
                neighbor,
                len,
            } => write!(
                f,
                "Neighbor index {neighbor} out of bounds for node {node} (adjacency list length {len})"
            ),
            CuthillMckeeError::CsrIndptrLengthMismatch { expected, actual } => write!(
                f,
                "CSR indptr length {actual} does not match nrows+1 (expected {expected})"
            ),
            CuthillMckeeError::CsrIndexOutOfBounds { index, ncols } => write!(
                f,
                "CSR indices value {index} out of bounds for ncols {ncols}"
            ),
            CuthillMckeeError::CscIndptrLengthMismatch { expected, actual } => write!(
                f,
                "CSC indptr length {actual} does not match ncols+1 (expected {expected})"
            ),
            CuthillMckeeError::CscIndexOutOfBounds { index, nrows } => write!(
                f,
                "CSC indices value {index} out of bounds for nrows {nrows}"
            ),
            CuthillMckeeError::SprsMatrixNotCsc => write!(f, "Matrix must be in CSC storage"),
        }
    }
}

impl std::error::Error for CuthillMckeeError {}

/// Computes the Cuthill–McKee ordering for a symmetric sparse matrix represented as an adjacency list.
///
/// The Cuthill–McKee algorithm reorders the nodes of a sparse, symmetric graph (or matrix)
/// to reduce its bandwidth, which is useful for efficient numerical computations.
/// The input is an adjacency list, where each index represents a node and the vector at that index
/// lists its neighbors. The output is a permutation of node indices representing the new order.
///
/// # Arguments
///
/// * `adj_list` - A slice of vectors, where each vector contains the indices of neighbors for each node.
///
/// # Returns
///
/// A vector of node indices representing the Cuthill–McKee ordering.
///
/// # Example
///
/// ```
/// use rma::cuthill_mckee;
/// // Graph: 0-1, 0-2, 2-3, 2-4
/// let adj_list = vec![
///     vec![1, 2], // 0
///     vec![0],    // 1
///     vec![0, 3, 4], // 2
///     vec![2],    // 3
///     vec![2],    // 4
/// ];
/// let order = cuthill_mckee(&adj_list).unwrap();
/// assert_eq!(order.len(), adj_list.len());
/// assert!(order.iter().all(|&i| i < adj_list.len()));
/// ```
pub fn cuthill_mckee(adj_list: &[Vec<usize>]) -> Result<Vec<usize>, CuthillMckeeError> {
    let n = adj_list.len();
    if n == 0 {
        return Err(CuthillMckeeError::EmptyAdjacencyList);
    }

    // Validate adjacency list
    if let Some((i, nbr)) = adj_list
        .iter()
        .enumerate()
        .flat_map(|(i, nbrs)| nbrs.iter().map(move |&nbr| (i, nbr)))
        .find(|&(_, nbr)| nbr >= n)
    {
        return Err(CuthillMckeeError::NeighborIndexOutOfBounds {
            node: i,
            neighbor: nbr,
            len: n,
        });
    }

    let mut visited = vec![false; n];
    let mut order = Vec::with_capacity(n);

    // Avoid repeated allocations by reusing the queue
    let mut queue = std::collections::VecDeque::with_capacity(n);

    // Pre-calculate degrees for faster comparisons
    let degrees: Vec<_> = adj_list.iter().map(|neighbors| neighbors.len()).collect();

    // Buffer to reuse for storing neighbors
    let mut neighbors_buffer = Vec::with_capacity(n);

    while order.len() < n {
        // Find the next start node (unvisited node with minimum degree)
        let start = (0..n)
            .into_iter()
            .filter(|&i| !visited[i])
            .min_by_key(|&i| degrees[i]);

        let Some(start) = start else {
            break;
        };

        visited[start] = true;
        queue.push_back(start);

        while let Some(node) = queue.pop_front() {
            order.push(node);

            // Reuse the neighbors buffer instead of creating a new vector
            neighbors_buffer.clear();

            // Collect unvisited neighbors with their degrees, marking as visited to avoid duplicates
            neighbors_buffer.extend(adj_list[node].iter().filter_map(|&nbr| {
                if !visited[nbr] {
                    visited[nbr] = true; // Mark visited immediately to avoid duplicates in queue
                    Some((nbr, degrees[nbr]))
                } else {
                    None
                }
            }));

            // Sort by degree using faster unstable sort
            neighbors_buffer.sort_unstable_by_key(|&(_, degree)| degree);

            // Add to queue in sorted order
            queue.extend(neighbors_buffer.iter().map(|&(nbr, _)| nbr));
        }
    }
    Ok(order)
}

/// Computes the Reverse Cuthill–McKee (RCM) ordering for a symmetric sparse matrix represented as an adjacency list.
///
/// The RCM algorithm reverses the Cuthill–McKee ordering, which often further reduces fill-in during matrix factorization.
/// The input is an adjacency list, and the output is a permutation of node indices representing the new order.
///
/// # Arguments
///
/// * `adj_list` - A slice of vectors, where each vector contains the indices of neighbors for each node.
///
/// # Returns
///
/// A vector of node indices representing the Reverse Cuthill–McKee ordering.
///
/// # Example
///
/// ```
/// use rma::reverse_cuthill_mckee;
/// // Graph: 0-1, 0-2, 2-3, 2-4
/// let adj_list = vec![
///     vec![1, 2], // 0
///     vec![0],    // 1
///     vec![0, 3, 4], // 2
///     vec![2],    // 3
///     vec![2],    // 4
/// ];
/// let order = reverse_cuthill_mckee(&adj_list).unwrap();
/// assert_eq!(order.len(), adj_list.len());
/// assert!(order.iter().all(|&i| i < adj_list.len()));
/// ```
pub fn reverse_cuthill_mckee(adj_list: &[Vec<usize>]) -> Result<Vec<usize>, CuthillMckeeError> {
    Ok(cuthill_mckee(adj_list)?.into_iter().rev().collect())
}

/// A simple Compressed Sparse Row (CSR) matrix representation for symmetric matrices.
#[derive(Debug, Clone)]
pub struct CsrMatrix<T = f64> {
    pub nrows: usize,
    pub ncols: usize,
    pub indptr: Vec<usize>,   // Row pointer (len = nrows + 1)
    pub indices: Vec<usize>,  // Column indices (len = nnz)
    pub data: Option<Vec<T>>, // Optional values (len = nnz)
}

/// Convert a CSR matrix to an adjacency list (for symmetric matrices).
pub fn csr_to_adjlist(csr: &CsrMatrix) -> Result<Vec<Vec<usize>>, CuthillMckeeError> {
    if csr.indptr.len() != csr.nrows + 1 {
        return Err(CuthillMckeeError::CsrIndptrLengthMismatch {
            expected: csr.nrows + 1,
            actual: csr.indptr.len(),
        });
    }
    let mut adj = vec![Vec::new(); csr.nrows];
    for row in 0..csr.nrows {
        let start = csr.indptr[row];
        let end = csr.indptr[row + 1];
        for idx in start..end {
            if idx >= csr.indices.len() {
                return Err(CuthillMckeeError::CsrIndptrLengthMismatch {
                    expected: csr.indices.len(),
                    actual: idx + 1,
                });
            }
            let col = csr.indices[idx];
            if col >= csr.ncols {
                return Err(CuthillMckeeError::CsrIndexOutOfBounds {
                    index: col,
                    ncols: csr.ncols,
                });
            }
            adj[row].push(col);
            if row != col {
                adj[col].push(row);
            }
        }
    }
    // Deduplicate neighbors
    adj.iter_mut().for_each(|nbrs| {
        nbrs.sort_unstable();
        nbrs.dedup();
    });
    Ok(adj)
}

/// Computes the Cuthill–McKee ordering for a symmetric sparse matrix in CSR format.
///
/// # Arguments
/// * `csr` - A reference to a CsrMatrix.
///
/// # Returns
/// A vector of node indices representing the Cuthill–McKee ordering.
///
/// # Example
/// ```
/// use rma::{CsrMatrix, cuthill_mckee_csr};
/// let csr = CsrMatrix {
///     nrows: 4,
///     ncols: 4,
///     indptr: vec![0, 1, 2, 3, 6],
///     indices: vec![0, 1, 2, 0, 2, 3],
///     // Not needed for Cuthill-McKee, yet if present, the algorithm will ignore it
///     data: None,
/// };
/// let order = cuthill_mckee_csr(&csr).unwrap();
/// assert_eq!(order.len(), csr.nrows);
/// ```
pub fn cuthill_mckee_csr(csr: &CsrMatrix) -> Result<Vec<usize>, CuthillMckeeError> {
    let adj_list = csr_to_adjlist(csr)?;
    cuthill_mckee(&adj_list)
}

/// A simple Compressed Sparse Column (CSC) matrix representation for symmetric matrices.
#[derive(Debug, Clone)]
pub struct CscMatrix<T = f64> {
    pub nrows: usize,
    pub ncols: usize,
    pub indptr: Vec<usize>,   // Column pointer (len = ncols + 1)
    pub indices: Vec<usize>,  // Row indices (len = nnz)
    pub data: Option<Vec<T>>, // Optional values (len = nnz)
}

/// Convert a CSC matrix to an adjacency list (for symmetric matrices).
pub fn csc_to_adjlist(csc: &CscMatrix) -> Result<Vec<Vec<usize>>, CuthillMckeeError> {
    if csc.indptr.len() != csc.ncols + 1 {
        return Err(CuthillMckeeError::CscIndptrLengthMismatch {
            expected: csc.ncols + 1,
            actual: csc.indptr.len(),
        });
    }
    let mut adj = vec![Vec::new(); csc.nrows];
    for col in 0..csc.ncols {
        let start = csc.indptr[col];
        let end = csc.indptr[col + 1];
        for idx in start..end {
            if idx >= csc.indices.len() {
                return Err(CuthillMckeeError::CscIndptrLengthMismatch {
                    expected: csc.indices.len(),
                    actual: idx + 1,
                });
            }
            let row = csc.indices[idx];
            if row >= csc.nrows {
                return Err(CuthillMckeeError::CscIndexOutOfBounds {
                    index: row,
                    nrows: csc.nrows,
                });
            }
            adj[row].push(col);
            if row != col {
                adj[col].push(row);
            }
        }
    }
    // Deduplicate neighbors
    adj.iter_mut().for_each(|nbrs| {
        nbrs.sort_unstable();
        nbrs.dedup();
    });
    Ok(adj)
}

/// Computes the Cuthill–McKee ordering for a symmetric sparse matrix in CSC format.
///
/// # Arguments
/// * `csc` - A reference to a CscMatrix.
///
/// # Returns
/// A vector of node indices representing the Cuthill–McKee ordering.
///
/// # Example
/// ```
/// use rma::{CscMatrix, cuthill_mckee_csc};
/// let csc = CscMatrix {
///     nrows: 4,
///     ncols: 4,
///     indptr: vec![0, 1, 2, 3, 6],
///     indices: vec![0, 1, 2, 0, 2, 3],
///     data: None,
/// };
/// let order = cuthill_mckee_csc(&csc).unwrap();
/// assert_eq!(order.len(), csc.nrows);
/// // Check it's a valid permutation
/// assert!(order.iter().all(|&i| i < csc.nrows));
/// ```
pub fn cuthill_mckee_csc(csc: &CscMatrix) -> Result<Vec<usize>, CuthillMckeeError> {
    let adj_list = csc_to_adjlist(csc)?;
    cuthill_mckee(&adj_list)
}

#[cfg(feature = "sprs")]
use sprs::CsMat;

#[cfg(feature = "sprs")]
impl<'a, N: Clone + 'a> From<&'a CsMat<N>> for CsrMatrix<N> {
    fn from(mat: &'a CsMat<N>) -> Self {
        CsrMatrix {
            nrows: mat.rows(),
            ncols: mat.cols(),
            indptr: mat.indptr().to_vec(),
            indices: mat.indices().to_vec(),
            data: Some(mat.data().to_vec()),
        }
    }
}

/// Computes the Cuthill–McKee ordering for a symmetric sparse matrix in sprs::CsMat (CSR) format.
///
/// # Arguments
/// * `mat` - A reference to a sprs::CsMat (CSR) matrix.
///
/// # Returns
/// A vector of node indices representing the Cuthill–McKee ordering.
///
/// # Example (with sprs)
/// ```
/// use sprs::CsMat;
/// use rma::cuthill_mckee_sprs;
/// let indptr = vec![0, 1, 2, 3, 6];
/// let indices = vec![0, 1, 2, 0, 2, 3];
/// let data = vec![1.0; 6];
/// let mat = CsMat::new((4, 4), indptr, indices, data);
/// let order = cuthill_mckee_sprs(&mat).unwrap();
/// assert_eq!(order.len(), mat.rows());
/// ```
#[cfg(feature = "sprs")]
pub fn cuthill_mckee_sprs<N>(mat: &CsMat<N>) -> Result<Vec<usize>, CuthillMckeeError> {
    let csr = CsrMatrix::from(mat);
    cuthill_mckee_csr(&csr)
}

#[cfg(feature = "sprs")]
impl<'a, N: Clone + 'a> From<&'a CsMat<N>> for CscMatrix<N> {
    fn from(mat: &'a CsMat<N>) -> Self {
        assert!(mat.is_csc(), "Matrix must be in CSC storage");
        CscMatrix {
            nrows: mat.rows(),
            ncols: mat.cols(),
            indptr: mat.indptr().to_vec(),
            indices: mat.indices().to_vec(),
            data: Some(mat.data().to_vec()),
        }
    }
}

/// Computes the Cuthill–McKee ordering for a symmetric sparse matrix in sprs::CsMat (CSC) format.
///
/// # Arguments
/// * `mat` - A reference to a sprs::CsMat (CSC) matrix.
///
/// # Returns
/// A vector of node indices representing the Cuthill–McKee ordering.
///
/// # Example (with sprs)
/// ```
/// use sprs::CsMat;
/// use rma::cuthill_mckee_sprs_csc;
/// let indptr = vec![0, 1, 2, 3, 6];
/// let indices = vec![0, 1, 2, 0, 2, 3];
/// let data = vec![1.0; 6];
/// let mat = CsMat::new_csc((4, 4), indptr, indices, data);
/// let order = cuthill_mckee_sprs_csc(&mat).unwrap();
/// assert_eq!(order.len(), mat.rows());
/// ```
#[cfg(feature = "sprs")]
pub fn cuthill_mckee_sprs_csc<N>(mat: &CsMat<N>) -> Result<Vec<usize>, CuthillMckeeError> {
    if !mat.is_csc() {
        return Err(CuthillMckeeError::SprsMatrixNotCsc);
    }
    let csc = CscMatrix::from(mat);
    cuthill_mckee_csc(&csc)
}

/// Computes the Reverse Cuthill–McKee (RCM) ordering for a symmetric sparse matrix in CSR format.
///
/// # Arguments
/// * `csr` - A reference to a CsrMatrix.
///
/// # Returns
/// A vector of node indices representing the Reverse Cuthill–McKee ordering.
///
/// # Example
/// ```
/// use rma::{CsrMatrix, reverse_cuthill_mckee_csr};
/// let csr = CsrMatrix {
///     nrows: 4,
///     ncols: 4,
///     indptr: vec![0, 1, 2, 3, 6],
///     indices: vec![0, 1, 2, 0, 2, 3],
///     // Not needed for Cuthill-McKee, yet if present, the algorithm will ignore it
///     data: None,
/// };
/// let order = reverse_cuthill_mckee_csr(&csr).unwrap();
/// assert_eq!(order.len(), csr.nrows);
/// ```
pub fn reverse_cuthill_mckee_csr(csr: &CsrMatrix) -> Result<Vec<usize>, CuthillMckeeError> {
    let adj_list = csr_to_adjlist(csr)?;
    reverse_cuthill_mckee(&adj_list)
}

/// Computes the Reverse Cuthill–McKee (RCM) ordering for a symmetric sparse matrix in CSC format.
///
/// # Arguments
/// * `csc` - A reference to a CscMatrix.
///
/// # Returns
/// A vector of node indices representing the Reverse Cuthill–McKee ordering.
///
/// # Example
/// ```
/// use rma::{CscMatrix, reverse_cuthill_mckee_csc};
/// let csc = CscMatrix {
///     nrows: 4,
///     ncols: 4,
///     indptr: vec![0, 1, 2, 3, 6],
///     indices: vec![0, 1, 2, 0, 2, 3],
///     // Not needed for Cuthill-McKee, yet if present, the algorithm will ignore it
///     data: None,
/// };
/// let order = reverse_cuthill_mckee_csc(&csc).unwrap();
/// assert_eq!(order.len(), csc.nrows);
/// ```
pub fn reverse_cuthill_mckee_csc(csc: &CscMatrix) -> Result<Vec<usize>, CuthillMckeeError> {
    let adj_list = csc_to_adjlist(csc)?;
    reverse_cuthill_mckee(&adj_list)
}

/// Computes the Reverse Cuthill–McKee (RCM) ordering for a symmetric sparse matrix in sprs::CsMat (CSR) format.
///
/// # Arguments
/// * `mat` - A reference to a sprs::CsMat (CSR) matrix.
///
/// # Returns
/// A vector of node indices representing the Reverse Cuthill–McKee ordering.
///
/// # Example (with sprs)
/// ```
/// use sprs::CsMat;
/// use rma::reverse_cuthill_mckee_sprs;
/// let indptr = vec![0, 1, 2, 3, 6];
/// let indices = vec![0, 1, 2, 0, 2, 3];
/// let data = vec![1.0; 6];
/// let mat = CsMat::new((4, 4), indptr, indices, data);
/// let order = reverse_cuthill_mckee_sprs(&mat).unwrap();
/// assert_eq!(order.len(), mat.rows());
/// ```
#[cfg(feature = "sprs")]
pub fn reverse_cuthill_mckee_sprs<N>(mat: &CsMat<N>) -> Result<Vec<usize>, CuthillMckeeError> {
    let csr = CsrMatrix::from(mat);
    reverse_cuthill_mckee_csr(&csr)
}

/// Computes the Reverse Cuthill–McKee (RCM) ordering for a symmetric sparse matrix in sprs::CsMat (CSC) format.
///
/// # Arguments
/// * `mat` - A reference to a sprs::CsMat (CSC) matrix.
///
/// # Returns
/// A vector of node indices representing the Reverse Cuthill–McKee ordering.
///
/// # Example (with sprs)
/// ```
/// use sprs::CsMat;
/// use rma::reverse_cuthill_mckee_sprs_csc;
/// let indptr = vec![0, 1, 2, 3, 6];
/// let indices = vec![0, 1, 2, 0, 2, 3];
/// let data = vec![1.0; 6];
/// let mat = CsMat::new_csc((4, 4), indptr, indices, data);
/// let order = reverse_cuthill_mckee_sprs_csc(&mat)?;
/// assert_eq!(order.len(), mat.rows());
/// ```
#[cfg(feature = "sprs")]
pub fn reverse_cuthill_mckee_sprs_csc<N>(mat: &CsMat<N>) -> Result<Vec<usize>, CuthillMckeeError> {
    if !mat.is_csc() {
        return Err(CuthillMckeeError::SprsMatrixNotCsc);
    }
    let csc = CscMatrix::from(mat);
    reverse_cuthill_mckee_csc(&csc)
}

/// Enum to specify which matrix type to read from Matrix Market
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MtxType {
    Csr,
    Csc,
}

/// Enum to wrap either a CsrMatrix or CscMatrix
#[derive(Debug, Clone)]
pub enum MatrixMarketMatrix<T = f64> {
    Csr(CsrMatrix<T>),
    Csc(CscMatrix<T>),
}

#[derive(Debug)]
pub enum MatrixMarketError {
    Io(std::io::Error),
    Parse(String),
    Corrupted(String),
}

impl std::fmt::Display for MatrixMarketError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MatrixMarketError::Io(e) => write!(f, "IO error: {}", e),
            MatrixMarketError::Parse(msg) => write!(f, "Parse error: {}", msg),
            MatrixMarketError::Corrupted(msg) => write!(f, "Corrupted Matrix Market file: {}", msg),
        }
    }
}

impl std::error::Error for MatrixMarketError {}

impl From<std::io::Error> for MatrixMarketError {
    fn from(e: std::io::Error) -> Self {
        MatrixMarketError::Io(e)
    }
}

/// Reads a symmetric Matrix Market file (coordinate format, integer or float, symmetric) and returns either a CsrMatrix or CscMatrix.
/// If the file contains values, they are loaded as f64 and stored in data: Some(`Vec<f64>`). If only structure, data is None.
pub fn read_matrix_market(
    path: &str,
    mtx_type: MtxType,
) -> Result<MatrixMarketMatrix<f64>, MatrixMarketError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let mut nrows = 0;
    let mut ncols = 0;
    let mut nnz = 0;
    let mut entries = Vec::new(); // (row, col, Option<value>)
    let mut header_parsed = false;
    let mut has_values = false;
    let mut line_num = 0;

    for line in reader.lines() {
        let line = line.map_err(MatrixMarketError::Io)?;
        let line = line.trim();
        line_num += 1;
        if line.is_empty() || line.starts_with('%') {
            continue;
        }
        if !header_parsed {
            // Header: rows cols nnz
            let parts: Vec<_> = line.split_whitespace().collect();
            if parts.len() != 3 {
                return Err(MatrixMarketError::Parse(format!(
                    "Header line must have 3 columns (nrows ncols nnz), got {} on line {}",
                    parts.len(),
                    line_num
                )));
            }
            nrows = parts[0].parse().map_err(|_| {
                MatrixMarketError::Parse(format!("Failed to parse nrows on line {}", line_num))
            })?;
            ncols = parts[1].parse().map_err(|_| {
                MatrixMarketError::Parse(format!("Failed to parse ncols on line {}", line_num))
            })?;
            nnz = parts[2].parse().map_err(|_| {
                MatrixMarketError::Parse(format!("Failed to parse nnz on line {}", line_num))
            })?;
            if nrows == 0 || ncols == 0 {
                return Err(MatrixMarketError::Corrupted(
                    "Matrix dimensions must be positive".to_string(),
                ));
            }
            header_parsed = true;
            continue;
        }
        // Entry: row col [value]
        let parts: Vec<_> = line.split_whitespace().collect();
        if parts.len() < 2 {
            return Err(MatrixMarketError::Parse(format!(
                "Entry line must have at least 2 columns (row col), got {} on line {}",
                parts.len(),
                line_num
            )));
        }
        let row: usize = parts[0].parse::<usize>().map_err(|_| {
            MatrixMarketError::Parse(format!("Failed to parse row index on line {}", line_num))
        })?;
        let col: usize = parts[1].parse::<usize>().map_err(|_| {
            MatrixMarketError::Parse(format!("Failed to parse col index on line {}", line_num))
        })?;
        if row == 0 || row > nrows || col == 0 || col > ncols {
            return Err(MatrixMarketError::Corrupted(format!(
                "Row or column index out of bounds on line {}: row={}, col={}, nrows={}, ncols={}",
                line_num, row, col, nrows, ncols
            )));
        }
        let value = if parts.len() > 2 {
            has_values = true;
            Some(parts[2].parse::<f64>().map_err(|_| {
                MatrixMarketError::Parse(format!("Failed to parse value on line {}", line_num))
            })?)
        } else {
            None
        };
        // Store as 0-based
        entries.push((row - 1, col - 1, value));
        if row != col {
            entries.push((col - 1, row - 1, value)); // symmetric
        }
    }

    if !header_parsed {
        return Err(MatrixMarketError::Parse(
            "Matrix Market header not found".to_string(),
        ));
    }
    // The number of unique (row,col) pairs in the file should match nnz
    let mut unique_entries = std::collections::HashSet::new();
    for &(row, col, _) in &entries {
        if row > nrows - 1 || col > ncols - 1 {
            return Err(MatrixMarketError::Corrupted(format!(
                "Entry index out of bounds: row={}, col={}, nrows={}, ncols={}",
                row, col, nrows, ncols
            )));
        }
        unique_entries.insert((row, col));
    }
    if unique_entries.len() / 2 > nnz {
        return Err(MatrixMarketError::Corrupted(format!(
            "More entries than declared nnz: found {}, declared {}",
            unique_entries.len() / 2,
            nnz
        )));
    }
    // Sort for construction
    entries.sort_by_key(|&(row, col, _)| (row, col));

    match mtx_type {
        MtxType::Csr => {
            // Build CSR
            let mut indptr = vec![0];
            let mut indices = Vec::new();
            let mut data = if has_values { Some(Vec::new()) } else { None };
            let mut current_row = 0;
            for &(row, col, value) in &entries {
                if row > current_row {
                    indptr.extend(std::iter::repeat(indices.len()).take(row - current_row));
                    current_row = row;
                }
                indices.push(col);
                if let Some(ref mut d) = data {
                    d.push(value.unwrap_or(1.0));
                }
            }
            indptr.resize(nrows + 1, indices.len());
            Ok(MatrixMarketMatrix::Csr(CsrMatrix {
                nrows,
                ncols,
                indptr,
                indices,
                data,
            }))
        }
        MtxType::Csc => {
            // Build CSC
            let mut indptr = vec![0];
            let mut indices = Vec::new();
            let mut data = if has_values { Some(Vec::new()) } else { None };
            let mut current_col = 0;
            // Sort by col, then row
            let mut col_entries = entries.clone();
            col_entries.sort_by_key(|&(row, col, _)| (col, row));
            for &(row, col, value) in &col_entries {
                if col > current_col {
                    indptr.extend(std::iter::repeat(indices.len()).take(col - current_col));
                    current_col = col;
                }
                indices.push(row);
                if let Some(ref mut d) = data {
                    d.push(value.unwrap_or(1.0));
                }
            }
            indptr.resize(ncols + 1, indices.len());
            Ok(MatrixMarketMatrix::Csc(CscMatrix {
                nrows,
                ncols,
                indptr,
                indices,
                data,
            }))
        }
    }
}

/// Writes a CsrMatrix or CscMatrix to a Matrix Market (.mtx) file in coordinate format, symmetric.
/// If data is Some, writes values; if None, writes structure-only (all values 1).
pub fn write_matrix_market(
    path: &str,
    matrix: &MatrixMarketMatrix<f64>,
) -> Result<(), std::io::Error> {
    let mut file = File::create(path)?;
    match matrix {
        MatrixMarketMatrix::Csr(csr) => {
            // Collect entries (row, col, value) for symmetric matrix
            let mut entries = Vec::new();
            for row in 0..csr.nrows {
                for (_, idx) in (csr.indptr[row]..csr.indptr[row + 1]).enumerate() {
                    let col = csr.indices[idx];
                    if row <= col {
                        let value = csr.data.as_ref().map(|d| d[idx]).unwrap_or(1.0);
                        entries.push((row, col, value));
                    }
                }
            }
            writeln!(file, "%%MatrixMarket matrix coordinate real symmetric")?;
            writeln!(file, "% Exported by rma")?;
            writeln!(file, "{} {} {}", csr.nrows, csr.ncols, entries.len())?;
            for (row, col, value) in entries {
                writeln!(file, "{} {} {}", row + 1, col + 1, value)?;
            }
        }
        MatrixMarketMatrix::Csc(csc) => {
            // Collect entries (row, col, value) for symmetric matrix
            let mut entries = Vec::new();
            for col in 0..csc.ncols {
                for (_, idx) in (csc.indptr[col]..csc.indptr[col + 1]).enumerate() {
                    let row = csc.indices[idx];
                    if row <= col {
                        let value = csc.data.as_ref().map(|d| d[idx]).unwrap_or(1.0);
                        entries.push((row, col, value));
                    }
                }
            }
            writeln!(file, "%%MatrixMarket matrix coordinate real symmetric")?;
            writeln!(file, "% Exported by rma")?;
            writeln!(file, "{} {} {}", csc.nrows, csc.ncols, entries.len())?;
            for (row, col, value) in entries {
                writeln!(file, "{} {} {}", row + 1, col + 1, value)?;
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::VecDeque;

    // Helper to convert node labels to indices for the A-F example
    // A=0, B=1, C=2, D=3, E=4, F=5
    // Edges: A-B, A-C, C-D, C-F, D-F, F-E
    fn abcf_example_adj_list() -> Vec<Vec<usize>> {
        vec![
            vec![1, 2],    // A: B, C
            vec![0],       // B: A
            vec![0, 3, 5], // C: A, D, F
            vec![2, 5],    // D: C, F
            vec![5],       // E: F
            vec![2, 3, 4], // F: C, D, E
        ]
    }

    // Check if 'order' is a valid permutation of 0..n
    fn is_valid_permutation(order: &[usize], n: usize) -> bool {
        let mut seen = vec![false; n];
        for &i in order {
            if i >= n || seen[i] {
                return false;
            }
            seen[i] = true;
        }
        seen.into_iter().all(|x| x)
    }

    // Check if the order is a valid BFS traversal with degree sorting
    fn is_valid_cuthill_mckee_order(adj: &[Vec<usize>], order: &[usize]) -> bool {
        let n = adj.len();
        if order.len() != n {
            return false;
        }
        let mut pos = vec![0; n];
        for (i, &v) in order.iter().enumerate() {
            pos[v] = i;
        }
        let mut visited = vec![false; n];
        let mut queue = VecDeque::new();
        let mut idx = 0;
        while idx < n {
            let start = order[idx];
            if visited[start] {
                idx += 1;
                continue;
            }
            visited[start] = true;
            queue.push_back(start);
            while let Some(node) = queue.pop_front() {
                // Collect unvisited neighbors
                let neighbors: Vec<_> = adj[node]
                    .iter()
                    .copied()
                    .filter(|&nbr| !visited[nbr])
                    .collect();
                // Should be sorted by degree
                let mut sorted = neighbors.clone();
                sorted.sort_by_key(|&nbr| adj[nbr].len());
                if neighbors != sorted {
                    return false;
                }
                for nbr in neighbors {
                    visited[nbr] = true;
                    queue.push_back(nbr);
                }
            }
            idx += 1;
        }
        true
    }

    // Compute the bandwidth of the matrix given the ordering
    fn bandwidth(adj: &[Vec<usize>], order: &[usize]) -> usize {
        let n = adj.len();
        let mut pos = vec![0; n];
        for (i, &v) in order.iter().enumerate() {
            pos[v] = i;
        }
        let mut max_bw = 0;
        for (i, nbrs) in adj.iter().enumerate() {
            for &j in nbrs {
                let bw = (pos[i] as isize - pos[j] as isize).abs() as usize;
                if bw > max_bw {
                    max_bw = bw;
                }
            }
        }
        max_bw
    }

    #[test]
    fn test_cuthill_mckee_abcf_properties() {
        let adj = abcf_example_adj_list();
        let order = cuthill_mckee(&adj).unwrap();
        assert!(is_valid_permutation(&order, adj.len()));
        assert!(is_valid_cuthill_mckee_order(&adj, &order));
        // Optionally: check bandwidth is reduced compared to identity
        let orig_bw = bandwidth(&adj, &(0..adj.len()).collect::<Vec<_>>());
        let new_bw = bandwidth(&adj, &order);
        assert!(new_bw <= orig_bw);
    }

    #[test]
    fn test_reverse_cuthill_mckee_abcf_properties() {
        let adj = abcf_example_adj_list();
        let order = reverse_cuthill_mckee(&adj).unwrap();
        assert!(is_valid_permutation(&order, adj.len()));
        // RCM is just the reverse, so it should also be a valid permutation
        // Optionally: check bandwidth is reduced compared to identity
        let orig_bw = bandwidth(&adj, &(0..adj.len()).collect::<Vec<_>>());
        let new_bw = bandwidth(&adj, &order);
        assert!(new_bw <= orig_bw);
    }

    #[test]
    fn test_cuthill_mckee_disconnected() {
        // Two disconnected components: 0-1, 2-3
        let adj = vec![
            vec![1], // 0
            vec![0], // 1
            vec![3], // 2
            vec![2], // 3
        ];
        let order = cuthill_mckee(&adj).unwrap();
        assert!(is_valid_permutation(&order, adj.len()));
        assert!(is_valid_cuthill_mckee_order(&adj, &order));
    }

    #[test]
    fn test_cuthill_mckee_single_node() {
        let adj = vec![vec![]];
        let order = cuthill_mckee(&adj).unwrap();
        assert!(is_valid_permutation(&order, adj.len()));
        assert!(is_valid_cuthill_mckee_order(&adj, &order));
    }

    #[test]
    fn test_cuthill_mckee_csc_equivalence() {
        // Equivalent to the CSR test, but using CSC
        let csc = CscMatrix {
            nrows: 4,
            ncols: 4,
            indptr: vec![0, 1, 2, 3, 6],
            indices: vec![0, 1, 2, 0, 2, 3],
            data: None,
        };
        let order = cuthill_mckee_csc(&csc).unwrap();
        assert_eq!(order.len(), csc.nrows);
        // Check it's a valid permutation
        assert!(order.iter().all(|&i| i < csc.nrows));
    }

    #[test]
    fn test_cuthill_mckee_invalid_adjlist() {
        // Neighbor index out of bounds
        let adj = vec![vec![1], vec![2]]; // node 1 has neighbor 2, but only 2 nodes
        let err = cuthill_mckee(&adj).unwrap_err();
        assert_eq!(
            err,
            CuthillMckeeError::NeighborIndexOutOfBounds {
                node: 1,
                neighbor: 2,
                len: 2
            }
        );
    }

    #[test]
    fn test_cuthill_mckee_empty_adjlist() {
        let adj: Vec<Vec<usize>> = vec![];
        let err = cuthill_mckee(&adj).unwrap_err();
        assert_eq!(err, CuthillMckeeError::EmptyAdjacencyList);
    }

    #[test]
    fn test_csr_to_adjlist_invalid_indptr() {
        let csr = CsrMatrix {
            nrows: 3,
            ncols: 3,
            indptr: vec![0, 1, 2], // should be len 4
            indices: vec![0, 1],
            data: None,
        };
        let err = csr_to_adjlist(&csr).unwrap_err();
        assert_eq!(
            err,
            CuthillMckeeError::CsrIndptrLengthMismatch {
                expected: 4,
                actual: 3
            }
        );
    }

    #[test]
    fn test_csr_to_adjlist_index_out_of_bounds() {
        let csr = CsrMatrix {
            nrows: 2,
            ncols: 2,
            indptr: vec![0, 1, 3],
            indices: vec![0, 2], // 2 is out of bounds
            data: None,
        };
        let err = csr_to_adjlist(&csr).unwrap_err();
        assert_eq!(
            err,
            CuthillMckeeError::CsrIndexOutOfBounds { index: 2, ncols: 2 }
        );
    }

    #[test]
    fn test_csc_to_adjlist_invalid_indptr() {
        let csc = CscMatrix {
            nrows: 3,
            ncols: 3,
            indptr: vec![0, 1, 2], // should be len 4
            indices: vec![0, 1],
            data: None,
        };
        let err = csc_to_adjlist(&csc).unwrap_err();
        assert_eq!(
            err,
            CuthillMckeeError::CscIndptrLengthMismatch {
                expected: 4,
                actual: 3
            }
        );
    }

    #[test]
    fn test_csc_to_adjlist_index_out_of_bounds() {
        let csc = CscMatrix {
            nrows: 2,
            ncols: 2,
            indptr: vec![0, 1, 3],
            indices: vec![0, 2], // 2 is out of bounds
            data: None,
        };
        let err = csc_to_adjlist(&csc).unwrap_err();
        assert_eq!(
            err,
            CuthillMckeeError::CscIndexOutOfBounds { index: 2, nrows: 2 }
        );
    }

    #[cfg(feature = "sprs")]
    #[test]
    fn test_cuthill_mckee_sprs_csc_not_csc() {
        use sprs::CsMat;
        let indptr = vec![0, 1, 2, 3, 6];
        let indices = vec![0, 1, 2, 0, 2, 3];
        let data = vec![1.0; 6];
        let mat = CsMat::new((4, 4), indptr, indices, data); // CSR, not CSC
        let err = cuthill_mckee_sprs_csc(&mat).unwrap_err();
        assert_eq!(err, CuthillMckeeError::SprsMatrixNotCsc);
    }
}
