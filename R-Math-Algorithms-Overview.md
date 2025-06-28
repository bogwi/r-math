# R-Math Algorithms Overview

## Table of Contents
- [The Chakravala Method](#the-chakravala-method)
  - [Overview](#overview)
  - [Purpose](#purpose)
  - [Algorithm Steps](#algorithm-steps)
  - [Example](#example)
  - [Key Features](#key-features)
  - [Applications](#applications)
  - [Sources](#sources)
- [Berlekamp's Root Finding Algorithm](#berlekamps-root-finding-algorithm)
  - [Overview](#overview-1)
  - [Key Features](#key-features-1)
  - [Applications](#applications-1)
  - [Summary Table](#summary-table)
  - [Sources](#sources-1)
- [Tonelli-Shanks Algorithm](#tonelli-shanks-algorithm)
  - [Overview](#overview-2)
  - [Algorithm Overview](#algorithm-overview)
  - [Steps of the Algorithm](#steps-of-the-algorithm)
  - [Example](#example-1)
  - [Applications](#applications-2)
  - [Sources](#sources-2)
- [Cuthill-McKee Algorithm](#cuthill-mckee-algorithm)
  - [Overview](#overview-3)
  - [Purpose](#purpose-1)
  - [Algorithm Steps](#algorithm-steps-1)
  - [Example](#example-2)
  - [Key Features](#key-features-2)
  - [Applications](#applications-3)
  - [Sources](#sources-3)
- [Freivalds' Algorithm](#freivalds-algorithm)
  - [Overview](#overview-4)
  - [How it works](#how-it-works)
  - [Error and Probability](#error-and-probability)
  - [Key Features](#key-features-3)
  - [Applications](#applications-4)
  - [Sources](#sources-4)

## The Chakravala method

### Overview
The Chakravala method is an ancient Indian cyclic algorithm for solving indeterminate quadratic equations, particularly Pell's equation ($x^2 - Ny^2 = 1$). Developed by Jayadeva (c. 950–1000 CE) and refined by Bhāskara II (c. 1114–1185 CE), it efficiently finds integer solutions through iterative composition and minimization&nbsp;[1][4].

### Purpose
The method solves equations of the form $x^2 - Ny^2 = k$ (where $N$ is a non-square integer and $k = \pm 1, \pm 2, \pm 4$), with Pell's equation ($k=1$) being a primary focus. It addresses the challenge of finding minimal integer solutions $(x, y)$&nbsp;[1][6].

### Algorithm Steps
1. **Initialization**: Start with a trivial triple $(a_0, b_0, k_0) = (1, 0, 1)$ satisfying $a_0^2 - N b_0^2 = k_0$&nbsp;[6][7].
2. **Iteration**: For each triple $(a, b, k)$:
   - **Select $m$**: Choose a positive integer $m$ such that:
     - $m \equiv -a \cdot b^{-1} \pmod{k}$ (using modular inverses or linear congruences)&nbsp;[6][7].
     - $|m^2 - N|$ is minimized&nbsp;[1][6].
   - **Update Triples**: Compute new values:
     $$
     a_{\text{next}} = \frac{a m + N b}{|k|}, \quad
     b_{\text{next}} = \frac{a + b m}{|k|}, \quad
     k_{\text{next}} = \frac{m^2 - N}{k}.
     $$
     This yields $a_{\text{next}}^2 - N b_{\text{next}}^2 = k_{\text{next}}$&nbsp;[1][6][7].
3. **Termination**: Repeat until $|k_{\text{next}}| = 1, 2,$ or $4$. If $k_{\text{next}} = \pm 1$, the solution is $(a_{\text{next}}, b_{\text{next}})$. For $|k_{\text{next}}| = 2$ or $4$, apply Brahmagupta's composition (Bhavana) to derive $(x, y)$ satisfying $x^2 - Ny^2 = 1$&nbsp;[1][3][6].

### Example 
Solving $x^2 - 67y^2 = 1$
1. **Initialization**: $(a_0, b_0, k_0) = (1, 0, 1)$.
2. **First Iteration**:
   - Find $m > 0$ minimizing $|m^2 - 67|$. Choose $m = 8$ (since $8^2 = 64 < 67$, $9^2 = 81 > 67$).
   - Update:
     $$
     a_1 = \frac{1 \cdot 8 + 67 \cdot 0}{1} = 8, \quad
     b_1 = \frac{1 + 0 \cdot 8}{1} = 1, \quad
     k_1 = \frac{8^2 - 67}{1} = -3.
     $$
     Thus, $8^2 - 67 \cdot 1^2 = -3$&nbsp;[1][6].
3. **Second Iteration**:
   - Solve $m \equiv -8 \cdot 1^{-1} \pmod{3} \implies m \equiv 1 \pmod{3}$. Minimize $|m^2 - 67|$ for $m = 1, 4, 7, \ldots$. Choose $m = 7$ (since $|49 - 67| = 18$ is minimal).
   - Update:
     $$
     a_2 = \frac{8 \cdot 7 + 67 \cdot 1}{3} = 41, \quad
     b_2 = \frac{8 + 1 \cdot 7}{3} = 5, \quad
     k_2 = \frac{7^2 - 67}{-3} = 6.
     $$
     Thus, $41^2 - 67 \cdot 5^2 = 6$&nbsp;[1][6].
4. **Final Steps**:
   - After 5 iterations, reach $(a_5, b_5, k_5) = (221, 27, -2)$.
   - Compose with itself via Brahmagupta's method:
     $$
     x = \frac{221^2 + 67 \cdot 27^2}{2} = 48842, \quad
     y = 221 \cdot 27 = 5967.
     $$
     Verify: $48842^2 - 67 \cdot 5967^2 = 1$&nbsp;[1][6].

### Key Features
- **Efficiency**: Minimizes $|m^2 - N|$ at each step, avoiding large numbers and ensuring rapid convergence&nbsp;[1][4].
- **Cyclic Structure**: The name "Chakravala" (Sanskrit for "wheel") reflects its iterative nature&nbsp;[1][3].
- **Historical Significance**: Solved $x^2 - 61y^2 = 1$ (solution: $x = 1766319049, y = 226153980$) centuries before European methods&nbsp;[1][4].
- **Minimality**: Generates solutions with minimal length and effort, proven to halt&nbsp;[6][7].

### Applications
- **Number Theory**: Solves Pell's equation, crucial for algebraic number theory and Diophantine analysis&nbsp;[1][6].
- **Approximation**: Delivers rational approximations to $\sqrt{N}$ (e.g., $\frac{48842}{5967} \approx \sqrt{67}$ with error $\sim 2 \times 10^{-9}$)&nbsp;[1].
- **Cryptography**: Underpins primality tests and integer factorization algorithms&nbsp;[6]

### Sources  
1. [Chakravala method - Wikipedia](https://en.wikipedia.org/wiki/Chakravala_method)  
2. [Pell equation: A generalization of continued fraction... (arXiv)](https://arxiv.org/html/2308.02742v3)  
3. [Exposition of the Chakrawal Method (PDF)](https://instavm.org/wp-content/uploads/2021/05/H17.pdf)  
4. [Rationale of the chakravala process (PDF)](https://scispace.com/pdf/rationale-of-the-chakravala-process-of-jayadeva-and-bhaskara-9rrtdhb1n1.pdf)  
5. [The Chakravala Method (PDF)](https://publications.azimpremjiuniversity.edu.in/1630/1/3_The%20Chakravala%20Method.pdf)  
6. [The Cakravāla Method for Solving Quadratic... (PDF)](https://www.imsc.res.in/~knr/acadimsc13/mds.pdf)  
7. [An Elementary Proof of the Halting Property... (arXiv)](https://arxiv.org/pdf/1406.6809.pdf)


## Berlekamp's Root Finding Algorithm

### Overview

Berlekamp's algorithm efficiently finds all roots of a univariate polynomial over a finite field GF(p), where p is a prime. The implementation in `berlekamp.rs` supports polynomials with coefficients and arithmetic in `u128`, and works for any prime modulus. The algorithm handles all cases, including repeated roots and non-monic polynomials, and is robust for both small and large primes.

The core idea is to use the Berlekamp matrix to find a basis for the nullspace of a certain linear operator, which reveals nontrivial factors of the polynomial. The algorithm recursively factors the polynomial into irreducible components, ultimately extracting all roots.

### Key Features

- **Supports any prime modulus**: Works for all prime fields GF(p), with p up to 2^128-1.
- **Handles all polynomial degrees**: Works for linear, quadratic, cubic, and higher-degree polynomials.
- **Robust to repeated roots**: Correctly finds roots even when the polynomial has repeated factors.
- **Automatic monic conversion**: Internally converts polynomials to monic form as needed.
- **Efficient arithmetic**: Uses Karatsuba multiplication for large polynomials, and modular arithmetic throughout.
- **Formal derivative and square-free factorization**: Handles cases where the formal derivative is zero (i.e., inseparable polynomials).
- **Comprehensive test suite**: Extensively tested on all edge cases, including irreducible, reducible, and repeated-root polynomials.
- **Example usage**:
  ```rust
  use rma::berlekamp::{Poly, berlekamp_roots};
  let p = 7u128;
  let f = Poly::new(vec![5, 0, 1], p); // x^2 + 5
  let mut roots = berlekamp_roots(&f);
  roots.sort();
  assert_eq!(roots, vec![3, 4]);
  ```

### Applications

- **Coding theory**: Decoding BCH and Reed–Solomon codes.
- **Cryptography**: Finding roots in cryptographic protocols over finite fields.
- **Computational algebra**: Factoring polynomials and finding solutions to equations in GF(p).
- **Number theory**: Studying properties of polynomials over finite fields.

### Summary Table

| Feature                | Supported |
|------------------------|-----------|
| Prime fields (GF(p))   | Yes       |
| Arbitrary degree       | Yes       |
| Repeated roots         | Yes       |
| Non-monic polynomials  | Yes       |
| Large primes           | Yes       |
| Efficient multiplication | Yes    |
| Handles f'(x) = 0      | Yes       |

### Sources

- [Encyclopedia: Berlekamp's Algorithm](https://encyclopedia.pub/entry/32373)
- [Wikipedia: Berlekamp's algorithm](https://en.wikipedia.org/wiki/Berlekamp%27s_algorithm)
- [USC lecture notes (PDF)](https://www.usc.es/regaca/mladra/docs/berlekamp.pdf)

## Tonelli-Shanks Algorithm

### Overview
The Tonelli–Shanks algorithm efficiently computes square roots modulo a prime $p$, solving equations of the form $r^2 \equiv n \mod p$. It is particularly useful in computational number theory and cryptography, where modular arithmetic is fundamental[1][4][5].

### Algorithm Overview
1. **Purpose**: Find $r$ such that $r^2 \equiv n \mod p$, where $p$ is an odd prime and $n$ is a quadratic residue modulo $p$.
2. **Key Insight**: Factor $p-1$ into $Q \cdot 2^S$, where $Q$ is odd. Use a quadratic non-residue $z$ to iteratively adjust the solution[1][4].
3. **Limitations**: Only works for primes. For composite moduli, the problem is equivalent to integer factorization[1][5].

---

### Steps of the Algorithm
1. **Check Legitimacy**:
   - Verify $n$ is a quadratic residue using Euler's criterion: $n^{(p-1)/2} \equiv 1 \mod p$. If not, no solution exists[4][5].

2. **Factor $p-1$**:
   - Write $p-1 = Q \cdot 2^S$, with $Q$ odd[4][6].

3. **Find a Quadratic Non-Residue $z$**:
   - Randomly test $z$ until $z^{(p-1)/2} \equiv -1 \mod p$. On average, this takes 2 trials[1][4].

4. **Initialize Variables**:
   - $m = S$, $c = z^Q \mod p$, $t = n^Q \mod p$, $r = n^{(Q+1)/2} \mod p$[4][6].

5. **Iterative Adjustment**:
   - While $t \neq 1$:
     - Find the smallest $j$ where $t^{2^j} \equiv 1 \mod p$.
     - Update $b = c^{2^{m-j-1}} \mod p$.
     - Set $m = j$, $c = b^2 \mod p$, $t = t \cdot c \mod p$, $r = r \cdot b \mod p$[4][6].

6. **Output**:
   - The solutions are $r$ and $p - r$[4][6].

---

### Example
Let $p = 41$, $n = 5$:
1. Factor $p-1 = 40 = 5 \cdot 2^3$ ($Q=5$, $S=3$).
2. Find $z=2$ (a non-residue since $2^{20} \equiv -1 \mod 41$).
3. Initialize:
   - $m=3$, $c=2^5 \mod 41 = 32$,
   - $t=5^5 \mod 41 = 40$,
   - $r=5^{(5+1)/2} \mod 41 = 28$.
4. Iterate until $t=1$:
   - $t=40 \neq 1$, find $j=1$ (since $40^2 \equiv 1 \mod 41$).
   - $b = 32^{2^{3-1-1}} = 32^2 \mod 41 = 40$.
   - Update $m=1$, $c=40^2 \mod 41 = 1$, $t=40 \cdot 1 \mod 41 = 40$, $r=28 \cdot 40 \mod 41 = 28$.
   - Next iteration: $t=40 \neq 1$, find $j=0$.
   - $b=1^{2^{1-0-1}}=1$, no change. Loop exits[4][6].
5. Solutions: $r=28$ and $41-28=13$.

---

### Applications
- Cryptography (e.g., elliptic curve point compression).
- Solving congruences in primality testing[5][7].

The Tonelli–Shanks algorithm balances efficiency and generality, making it a cornerstone in computational number theory.

### Sources
[1] Tonelli–Shanks algorithm - Wikipedia https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm
[2] Tonelli-Shanks algorithm - Rosetta Code https://rosettacode.org/wiki/Tonelli-Shanks_algorithm
[3] Tonelli-Shanks アルゴリズムの理解と実装(1) - Qiita https://qiita.com/taiyaki8926/items/f62f534d43ff006129f7
[4] Computing square roots modulo a prime with the Tonelli-Shanks ... https://zerobone.net/blog/math/tonelli-shanks/
[5] Tonelli–Shanks algorithm - HandWiki https://handwiki.org/wiki/Tonelli%E2%80%93Shanks_algorithm
[6] Tonelli-Shanks アルゴリズムの理解と実装(2) #Python - Qiita https://qiita.com/taiyaki8926/items/6720a1c19e9c9afbfd12
[7] Finding Mod-p Square Roots with the Tonelli-Shanks Algorithm https://www.youtube.com/watch?v=d7ZFCf95MAQ
[8] mod_sqrtについて.md - GitHub Gist https://gist.github.com/kirika-comp/2562ffecbd4db344d66bf7c08a40b548

## Cuthill-McKee Algorithm

### Overview
The Cuthill–McKee (CM) algorithm is a graph-based method to reorder rows and columns of symmetric sparse matrices, reducing their bandwidth for efficient numerical computations like Gaussian elimination. Developed by Elizabeth Cuthill and James McKee, it operates by renumbering nodes in a breadth-first search (BFS) pattern while prioritizing neighbors with lower degrees [[2]](https://en.wikipedia.org/wiki/Cuthill%E2%80%93McKee_algorithm) [[5]](https://student.cs.uwaterloo.ca/~cs475/CS475-Lecture06.pdf). The **reverse Cuthill–McKee (RCM)** variant, which reverses the CM ordering, often achieves even better results by minimizing fill-in during matrix factorization [[2]](https://en.wikipedia.org/wiki/Cuthill%E2%80%93McKee_algorithm) [[8]](https://epubs.siam.org/doi/10.1137/0713020).

### Purpose
- **Bandwidth reduction**: Converts sparse matrices into banded forms with narrower non-zero diagonals, optimizing memory and computational efficiency [[2]](https://en.wikipedia.org/wiki/Cuthill%E2%80%93McKee_algorithm) [[6]](http://juliafem.org/examples/2017-08-29-reordering-nodes-with-the-RCM-algorithm).
- **Fill-in minimization**: RCM specifically reduces new non-zero entries (fill-in) during matrix decomposition, critical for solving linear systems [[8]](https://epubs.siam.org/doi/10.1137/0713020) [[5]](https://student.cs.uwaterloo.ca/~cs475/CS475-Lecture06.pdf).

### Algorithm Steps
1. **Starting node selection**:  
   Choose a peripheral node (e.g., one with the lowest degree) as the root [[2]](https://en.wikipedia.org/wiki/Cuthill%E2%80%93McKee_algorithm) [[5]](https://student.cs.uwaterloo.ca/~cs475/CS475-Lecture06.pdf).
2. **Breadth-first traversal**:  
   Visit nodes level by level, appending neighbors of the current node to the ordering list [[2]](https://en.wikipedia.org/wiki/Cuthill%E2%80%93McKee_algorithm) [[4]](https://www.youtube.com/watch?v=E6yIXP1OYeQ).
3. **Ordering within levels**:  
   Sort adjacent nodes by:
   - **Predecessor order**: Nodes connected to earlier-numbered predecessors first.
   - **Degree**: Break ties using ascending node degree [[2]](https://en.wikipedia.org/wiki/Cuthill%E2%80%93McKee_algorithm) [[5]](https://student.cs.uwaterloo.ca/~cs475/CS475-Lecture06.pdf).
4. **Repeat**:  
   Continue until all nodes are ordered. For disconnected graphs, apply the algorithm separately to each component [[5]](https://student.cs.uwaterloo.ca/~cs475/CS475-Lecture06.pdf).

### Example
Consider a graph with nodes **A–F** (as in [[5]](https://student.cs.uwaterloo.ca/~cs475/CS475-Lecture06.pdf)):
- **CM ordering**: Starts at **A**, visits neighbors **B** (degree 1) and **C** (degree 3). Next, **C**'s neighbors **D** (degree 2) and **F** (degree 3), followed by **F**'s neighbor **E**. Result: `A → B → C → D → F → E`.
- **RCM ordering**: Reverse the CM sequence: `E → F → D → C → B → A`. This reduces fill-in to zero in elimination [[5]](https://student.cs.uwaterloo.ca/~cs475/CS475-Lecture06.pdf).

### Key Features
- **Reverse Cuthill–McKee (RCM)**:  
  Superior to CM for Gaussian elimination, as it often minimizes fill-in [[8]](https://epubs.siam.org/doi/10.1137/0713020) [[5]](https://student.cs.uwaterloo.ca/~cs475/CS475-Lecture06.pdf).
- **Optimality on trees**: RCM guarantees no fill-in when applied to tree structures [[5]](https://student.cs.uwaterloo.ca/~cs475/CS475-Lecture06.pdf).
- **Heuristic nature**: Performance depends on the starting node; peripheral nodes yield better results [[2]](https://en.wikipedia.org/wiki/Cuthill%E2%80%93McKee_algorithm) [[6]](http://juliafem.org/examples/2017-08-29-reordering-nodes-with-the-RCM-algorithm).

### Applications
- **Finite element methods**: Mesh node reordering for sparse matrix solvers [[6]](http://juliafem.org/examples/2017-08-29-reordering-nodes-with-the-RCM-algorithm).
- **Circuit simulation**: Reducing interconnection complexity [[7]](https://www.kurims.kyoto-u.ac.jp/~kyodo/kokyuroku/contents/pdf/1375-27.pdf).
- **Linear algebra**: Bandwidth-sensitive algorithms like Cholesky decomposition [[8]](https://epubs.siam.org/doi/10.1137/0713020).

The algorithm's efficiency makes it a cornerstone in computational methods requiring sparse matrix optimizations.


### Sources
[1] カットヒル・マキー法 - Wikipedia https://ja.wikipedia.org/wiki/%E3%82%AB%E3%83%83%E3%83%88%E3%83%92%E3%83%AB%E3%83%BB%E3%83%9E%E3%82%AD%E3%83%BC%E6%B3%95
[2] Cuthill–McKee algorithm - Wikipedia https://en.wikipedia.org/wiki/Cuthill%E2%80%93McKee_algorithm
[3] cuthill_mckee_ordering - boostjp https://boostjp.github.io/archive/boost_docs/libs/graph/cuthill_mckee_ordering.html
[4] Cuthill Mckee Algorithm by Dr. Aarthi - YouTube https://www.youtube.com/watch?v=E6yIXP1OYeQ
[5] [PDF] Lecture 06 - Matrix Re-Ordering https://student.cs.uwaterloo.ca/~cs475/CS475-Lecture06.pdf
[6] Reordering nodes with the Reverse Cuthill-McKee algorithm http://juliafem.org/examples/2017-08-29-reordering-nodes-with-the-RCM-algorithm
[7] [PDF] バンド幅問題に対する遺伝的アルゴリズムの交叉方法の提案 https://www.kurims.kyoto-u.ac.jp/~kyodo/kokyuroku/contents/pdf/1375-27.pdf
[8] Comparative Analysis of the Cuthill–McKee and the ... - SIAM.org https://epubs.siam.org/doi/10.1137/0713020

## Freivalds' Algorithm

### Overview
Freivalds' algorithm is a probabilistic randomized method to verify if the product of two $n \times n$ matrices $A$ and $B$ equals a third matrix $C$ without computing the full product $A \times B$ explicitly [[1]](https://en.wikipedia.org/wiki/Freivalds'_algorithm) [[6]](http://ravi-bhide.blogspot.com/2012/05/freivalds-algorithm.html) [[8]](https://iq.opengenus.org/freivalds-algorithm/).

### How it works
- **Input:** Three $n \times n$ matrices $A, B, C$.
- Choose a random $n \times 1$ vector $r$ with entries 0 or 1, selected uniformly at random.
- Compute the vectors $A \times (B \times r)$ and $C \times r$.
- If these two vectors are different, output **No** (i.e., $A \times B \neq C$).
- If they are equal, output **Yes** (i.e., $A \times B = C$ with high probability).

This verification takes $O(n^2)$ time, which is faster than the $O(n^{2.3729})$ or $O(n^3)$ time needed for full matrix multiplication [[1]](https://en.wikipedia.org/wiki/Freivalds'_algorithm) [[6]](http://ravi-bhide.blogspot.com/2012/05/freivalds-algorithm.html).

### Error and Probability
- If $A \times B = C$, the algorithm always returns **Yes** (no false negatives).
- If $A \times B \neq C$, the probability that the algorithm incorrectly returns **Yes** is at most $1/2$ per run (one-sided error).

### Key Features
- **Randomized verification:** Monte Carlo algorithm.
- **Linear time** in the size of the matrices.

### Applications
- **Fast verification of matrix products**, especially in large-scale computations.

Freivalds' algorithm efficiently verifies matrix multiplication correctness with high probability by using randomization and vector-matrix multiplications, making it much faster than recomputing the full product [[1]](https://en.wikipedia.org/wiki/Freivalds'_algorithm) [[6]](http://ravi-bhide.blogspot.com/2012/05/freivalds-algorithm.html) [[8]](https://iq.opengenus.org/freivalds-algorithm/).


### Sources
[1] Freivalds' algorithm - Wikipedia https://en.wikipedia.org/wiki/Freivalds'_algorithm
[2] [PDF] 1 Verifying Matrix Product (Freivalds' algorithm) - Avah Banerjee https://www.avahbanerjee.com/files/Teaching/CS6200/CS_6200_Lec1-3.pdf
[3] Freivalds' Algorithm - YouTube https://www.youtube.com/watch?v=3E5CZ9ytYAY
[4] Freivalds' algorithm - HandWiki https://handwiki.org/wiki/Freivalds'_algorithm
[5] Freivalds' Algorithm | PDF | Areas Of Computer Science - Scribd https://www.scribd.com/document/380296415/Freivalds-Algorithm
[6] Freivalds' Algorithm - Theory behind the Technology http://ravi-bhide.blogspot.com/2012/05/freivalds-algorithm.html
[7] [PDF] Lecture 13: Dimensionality Reduction https://www.cl.cam.ac.uk/teaching/1920/Probablty/materials/lec13.pdf
[8] Freivalds' algorithm for verifying Matrix Multiplication - OpenGenus IQ https://iq.opengenus.org/freivalds-algorithm/
