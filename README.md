# MATLAB for Linear Algebra (V2)

## What's New in V2
This toolkit has been upgraded for the AY25/26 semester to prioritize mathematical accuracy and ease of use for NUS Engineering students.
* **Symbolic Math Integration:** Methods now natively support and return `sym` objects, providing exact fractional results (e.g., `1/3`) instead of messy floating-point decimals.
* **Step-by-Step Outputs:** Functions like `gramSchmidt` now provide the exact algebraic working, intermediate projections, and inner products for exam verification.
* **Subspace Extraction:** Added a new utility to automatically find a linearly independent basis for the span of any given set of vectors.
* **Chapter 7 Visualizer:** Added a Phase Portrait plotter to easily identify Sinks, Sources, Saddles, and Spirals geometrically.

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
  The set is othogonal
```
---

## 📈 Chapter 7 Phase Portrait Gallery

Use the `m.plotPhasePortrait(A)` function to instantly visualize the stability of 2x2 differential systems. The thick red lines indicate the exact straight-line solutions (eigenvectors).

### 1. The Saddle Point (Unstable)
One positive eigenvalue, one negative eigenvalue. The arrows flow inward along one eigenvector and outward along the other.
```MATLAB
A = sym([1 2; 2 1]);
m.plotPhasePortrait(A);
```
### 2. The Sink / Stable Node

Both eigenvalues are negative. All trajectories are pulled directly into the origin over time.
```MATLAB
A = sym([-2 1; 1 -2]);
m.plotPhasePortrait(A);
```
### 3. The Source / Unstable Node

Both eigenvalues are positive. All trajectories are pushed outward and away from the origin.
```MATLAB
A = sym([2 1; 1 2]);
m.plotPhasePortrait(A);
```
### 4. The Spiral

Complex eigenvalues (λ=a±bi). If the real part a is negative, it spirals inward (Stable). If a is positive, it spirals outward (Unstable).
```MATLAB
A = sym([-1 1; -5 -1]);
m.plotPhasePortrait(A);
```
