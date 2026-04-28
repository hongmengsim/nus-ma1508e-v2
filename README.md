# MATLAB for Linear Algebra (V2)

## What's New in V2
This toolkit has been upgraded for the AY25/26 semester to prioritize mathematical accuracy and ease of use for NUS Engineering students.
* **Symbolic Math Integration:** Methods now natively support and return `sym` objects, providing exact fractional results (e.g., `1/3`) instead of messy floating-point decimals.
* **Step-by-Step Outputs:** Functions like `gramSchmidt` now provide the exact algebraic working, intermediate projections, and inner products for exam verification.
* **Subspace Extraction:** Added a new utility to automatically find a linearly independent basis for the span of any given set of vectors.
* **Chapter 7 Visualizer:** Added a Phase Portrait plotter to easily identify Sinks, Sources, Saddles, and Spirals geometrically.
* **Rigorous QR Decomposition:** Upgraded `gramSchmidt` to the Modified Gram-Schmidt algorithm, ensuring the $R$ matrix is strictly upper-triangular and $V = QR$ reconstruction is successful.
* **Parametric LSS Solver:** The Least Squares solver now identifies infinite solution spaces and formats them as $\mathbf{x} = \mathbf{x}_p + s_1 \mathbf{v}_1 + \dots$

## Setup 
*Note: This toolkit requires the **Symbolic Math Toolbox** to be installed in your MATLAB Add-Ons.*

```sh
git clone https://github.com/hongmengsim/nus-ma1508e-v2.git
cd nus-ma1508e-v2
```
Usage
Initialise the class instance in your Command Window:
```MATLAB
m = MA1508E();
```
To ensure exact fractional outputs, it is highly recommended to pass your matrices as symbolic objects using the `sym()` function:
```MATLAB
A = sym([1 0 -2; 1 1 1; 1 -1 1]);
m.getPD(A);
```
**Example Verification**

$$
S = \begin{Bmatrix} 
\begin{pmatrix} 1 \cr 1 \cr 1 \end{pmatrix}, 
\begin{pmatrix} 0 \cr 1 \cr -1 \end{pmatrix}, 
\begin{pmatrix} -2 \cr 1 \cr 1 \end{pmatrix} 
\end{Bmatrix}
$$

To check if the set of vectors S is orthogonal:
```MATLAB
A = sym([1 0 -2; 1 1 1; 1 -1 1]);
m.isOrthogonalSet(A);
```
Output:
```MATLAB
  The set is orthogonal
```
---
## 📚 Detailed API Reference
Click on any method below to expand its full documentation, parameters, and conditions.

## ✅ How to Verify Your Results
The V2.1 engine is designed for easy verification. Run these checks in your Command Window:

1. **Check Orthonormality:** `simplify(e' * e)` should yield the **Identity Matrix**.
2. **Check QR Reconstruction:** `simplify(e * r)` should yield your original matrix **V**.
3. **Check LSS Accuracy:** `simplify(transpose(A) * (b - A*x_p))` should yield a **Zero Vector**.

### Chapter 1: Linear Systems
<details>
<summary><code>isValidERO(str)</code></summary>

- **`str`**: String containing the elementary row operation.
- **Description**: Validates the syntax of an ERO string.
- **V2.1 Upgrades**: 
  - **Space-Insensitive**: The parser now accepts both standard (`R2 + 2R1`) and compact (`R2+2R1`) formats for faster input during timed labs.
  - **Rules**: Row prefixes must remain an uppercase **R**. Supports exact fractions (e.g., `1/2R1`).
</details>

<details>
<summary><code>generateElemAdd(data)</code> / <code>generateElemSwap(data)</code></summary>

- **`data`**: Struct containing the elementary row operation.
- **Description**: Returns the elementary matrix corresponding to the elementary row operation.
</details>

<details>
<summary><code>generateElemMatrix(str, n)</code></summary>

- **`str`**: Elementary row operation to be performed (e.g., `R2+0.5R3`, `3R1`, `R4-1R2`, `R2SR3`).
- **`n`**: The number of rows of the target matrix.
- **Description**: Generates the corresponding elementary matrix used to perform the row operation via pre-multiplication.
- **V2.1 Upgrades**: Fully supports space-insensitive strings and seamlessly parses symbolic fractions into exact matrix elements.
</details>

<details>
<summary><code>performERO(A)</code></summary>
  
### `performERO(A)`
**Description:** Launches an interactive, step-by-step console loop that allows the user to apply Elementary Row Operations (EROs) to a given matrix. The matrix is updated and printed to the console after every successful operation. Type `exit` or `quit` to break the loop.

> **Note for V2:**
  - Input strings no longer require quotation marks in the console. Furthermore, passing a symbolic matrix (`sym()`) ensures all intermediate fractions remain exact.
  - Now utilizes the space-insensitive parser for rapid data entry. 
  - **Tip:** Pass your matrix as `A = sym(A)` before running this method to ensure all intermediate row additions and scalar multiplications remain in exact fractional form instead of decimals.

#### **Parameters**
- **`A`** *(Matrix | sym)*: The target matrix you want to perform operations on. 

#### **Supported Operation Syntax**
The console input parser requires specific formatting. Use the letter `R` followed by the row number.

| Operation Type | Syntax | Description | Example |
| :--- | :--- | :--- | :--- |
| **Row Swapping** | `Ri S Rj` | Swaps Row *i* with Row *j* | `R1 S R2` |
| **Scalar Multiplication** | `cRi` | Multiplies Row *i* by a scalar value *c* | `3R1` or `1/2R3` |
| **Row Addition** | `Rj + cRi` | Adds *c* times Row *i* to Row *j* | `R2 + 2R1` |
| **Row Subtraction** | `Rj - cRi` | Subtracts *c* times Row *i* from Row *j* | `R4 - R2` |

#### **Example Usage Session**

**1. Initialize in Command Window:**
```matlab
>> m = MA1508E();
>> A = sym([1 2 3; 2 4 8]);
>> m.performERO(A);
```
</details>

### Chapter 2: Matrices
<details>
<summary><code>leftInverse(A)</code></summary>

- **`A`**: $m \times n$ matrix to check.
- **Description**: Returns the left inverse of $A$ if it exists. 
- **Condition**: $m > n$ for left inverse to exist.
</details>

<details>
<summary><code>rightInverse(A)</code></summary>

- **`A`**: $m \times n$ matrix to check.
- **Description**: Returns the right inverse of $A$ if it exists. 
- **Condition**: $m < n$ for right inverse to exist.
</details>

### Chapter 5: Orthogonal Projection
<details>
<summary><code>isOrthogonalSet(S)</code> / <code>isOrthonormalSet(S)</code></summary>

- **`S`**: Set of vectors in matrix form.
- **Description**: Prints whether the given set of vectors is strictly orthogonal or orthonormal.
</details>

<details>
<summary><code>isOrthogonalTo(S, target)</code></summary>

- **`S`**: Set of vectors in matrix form.
- **`target`**: The vector to be checked.
- **Description**: Prints whether `target` is orthogonal to the subspace spanned by the set `S`.
</details>

<details>
<summary><code>toOrthonormalSet(OG)</code></summary>

- **`OG`**: Orthogonal set of vectors in matrix form.
- **Description**: Normalises the set of vectors. Returns `ON`, the corresponding orthonormal set.
</details>

<details>
<summary><code>dotWithSet(S, v)</code></summary>

- **`S`**: Set of vectors in matrix form.
- **`v`**: Column vector with the same number of rows as `S`.
- **Description**: Returns the dot product between each column of `S` and `v`.
</details>

<details>
<summary><code>orthogonalProj(S, w)</code></summary>

- **`S`**: Set of vectors in matrix form.
- **`w`**: The vector to be projected.
- **Description**: Returns the exact orthogonal projection of `w` onto the span of `S`.
</details>

<details>
<summary><code>gramSchmidt(v, showSteps)</code></summary>

- **`v`**: The basis to be converted into an orthonormal basis.
- **`showSteps`**: Boolean (Default = `true`). Set to false to hide algebraic working.
- **Description**: Performs **Modified Gram-Schmidt** orthonormalization.
- **V2.1 Upgrades**:
  - **Surd Support**: Returns exact symbolic square roots (e.g., $\frac{\sqrt{6}}{3}$) instead of decimal approximations.
  - **QR Decomposition**: Returns `[e, r]` where `e` is the orthonormal $Q$ matrix and `r` is the upper-triangular $R$ matrix. Guaranteed to satisfy $V = QR$ and $Q^T Q = I$.
</details>

<details>
<summary><code>calcLSS(A, b)</code></summary>

- **`A`**: The matrix to calculate the least squares solution.
- **`b`**: The column vector.
- **Description**: Solves the Normal Equations $(A^T A)\mathbf{x} = A^T \mathbf{b}$ symbolically.
- **V2.1 Upgrades**:
  - **Parametric Vector Form**: If the system is underdetermined (infinite solutions), the output automatically identifies the null space and formats the general solution as a geometric translation: $\mathbf{x} = \mathbf{x}_p + s_1 \mathbf{v}_1 + \dots$
  - Returns `[x_gen, basis, x_p]` for immediate workspace use.
</details>

### Chapter 6: Diagonalisation
<details>
<summary><code>getEigenvalues(A)</code></summary>

- **`A`**: The matrix corresponding to the differential system.
- **Description**: Returns both real and complex eigenvalues corresponding to `A`.
</details>

<details>
<summary><code>getEigenvector(A, lambda, output)</code></summary>

- **`A`**: The matrix corresponding to the differential system.
- **`lambda`**: An eigenvalue of `A`.
- **`output`**: Boolean (Default = `true`). Show steps to obtain the eigenvector.
- **Description**: Returns the basis for the eigenspace associated to `lambda`.
- **V2.1 Upgrades**: 
  - **Symbolic Null Space Analysis**: By leveraging the symbolic engine, this method now reliably extracts exact basis vectors even for complex eigenvalues ($\lambda = a \pm bi$) without falling victim to floating-point truncation errors common in standard numerical solvers.
</details>

<details>
<summary><code>isAssociatedEigenvalue(A, lambda)</code></summary>

- **`A`**: The matrix corresponding to the differential system.
- **`lambda`**: Eigenvalue to verify.
- **Description**: Returns `true` if lambda is an eigenvalue associated to `A`, `false` otherwise.
</details>

<details>
<summary><code>getGeneralisedEigenvector(A, lambda)</code></summary>

- **`A`**: The matrix corresponding to the differential system.
- **`lambda`**: An eigenvalue of `A`.
- **Description**: Returns the matrix which gives you the valid generalised eigenvector for defective matrices.
</details>

### Chapter 7: System of Linear Differential Equations
<details>
<summary><code>generateInitialConditions(n)</code></summary>

- **`n`**: The number of initial conditions of a differential system.
- **Description**: Returns formatted string to be input as initial conditions in MATLAB's `dsolve`.
</details>

<details>
<summary><code>solveDifferentialSystem(A, isInitial)</code></summary>

- **`A`**: The matrix corresponding to the differential system.
- **`isInitial`**: Boolean. `false` to obtain general solution of the differential system, `true` for particular solution.
- **Description**: Returns a struct containing the symbolic solution to the system.
</details>

<details>
<summary><code>plotPhasePortrait(A)</code></summary>
  
Use the `m.plotPhasePortrait(A)` function to instantly visualize the stability of 2x2 differential systems. The thick red lines indicate the exact straight-line solutions (eigenvectors).

### 1. The Saddle Point (Unstable)

One positive eigenvalue, one negative eigenvalue. The arrows flow inward along one eigenvector and outward along the other.
```MATLAB
A = sym([1 2; 2 1]);
m.plotPhasePortrait(A);
```
![Saddle](saddle.svg)

### 2. The Sink / Stable Node

Both eigenvalues are negative. All trajectories are pulled directly into the origin over time.
```MATLAB
A = sym([-2 1; 1 -2]);
m.plotPhasePortrait(A);
```
![Sink](sink.svg)

### 3. The Source / Unstable Node

Both eigenvalues are positive. All trajectories are pushed outward and away from the origin.
```MATLAB
A = sym([2 1; 1 2]);
m.plotPhasePortrait(A);
```
![Source](source.svg)

### 4. The Spiral

Complex eigenvalues ($\lambda = a \pm bi$). If the real part a is negative, it spirals inward (Stable). If a is positive, it spirals outward (Unstable).
```MATLAB
A = sym([-1 1; -5 -1]);
m.plotPhasePortrait(A);
```
![Spiral](spiral.svg)
</details>

