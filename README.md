# Implementation of Numerical Methods 

Hey everyone! This repository contains my progress on university tasks related to numerical methods. I hope you'll find something useful here!

## Task 1: Condition Number, LU-decomposition, Inverse Matrix

Write a program to calculate the condition number of a square matrix, using norms:

```math
||x||_{\alpha}, \qquad where \quad \alpha = 1, \infty, E
```

The calculation of an inverse matrix should include LU-decomposition. Find the condition number of matrices following this pattern:

```math
a_{i,j} = min(i,j), \quad i=1,2,3,...,n, \quad j=1,2,3,...,n
```

Ensure that the inverse matrix is tridiagonal, with non-diagonal elements:
```math
a^{-1}_{i, i \pm 1}=-1
```
and diagonal elements:
```math
a^{-1}_{i,i} =
\begin{cases}
2, & i=1,2,...,n-1 \
1, & i=n
\end{cases}
```
Plot the dependence of the calculation time of the inverse matrix on the size of the matrix.

### Description of the Solution

The `Task 1` directory contains 4 files:

- `external.py`
- `cond_num.py`
- `results.csv`
- `one_minute.png`

`external.py` contains 8 functions that run the entire process:

- `generate_matrix`: generates the matrix according to the task
- `matrix_multiplication`: multiplies two matrices using the basic algorithm
- `lu_decompose`: decomposes a matrix according to the Gaussian algorithm of LU-decomposition
- `invert`: returns the inverse matrix
- `norm_1`, `norm_inf`, `norm_euclid`: calculate 3 kinds of norms of a matrix
- `condition_numbers`: returns three condition numbers for each kind of norm

`cond_num.py` is the main program that accomplishes the task. It retrieves the condition numbers and plots the dependency.

`results.csv` is the file with those particular condition numbers.

`one_minute.png` is the actual plot of the time and size relation.

> The execution time of the program is limited by the number of minutes mentioned at the beginning (this number could be a float).

`Plato.png` is a screenshot of a plateau during a one-minute test.

## Task 2: Iterative Methods, Seidel's Method, Alternating-triangular Method

### Part 1

Write a program that implements the approximate solution of a system of linear algebraic equations using Seidel's method. Use the program to find the solution of the following system:

```math
Ax = f
```
where elements of A are:

```math
a_{i i}=2, \quad a_{i i+1} = -1 - \alpha, \quad a_{i i-1}
```
and elements of f are:

```math
f_{1} = 1 - \alpha, \quad f_{i} = 0, \quad f_{n} = 1 + \alpha
```
The accurate solution is:

```math
x_{i} = 1, \quad i = 10, 20, 30, ..., 100
```
Explore the dependence of the iteration count on n and

```math
\alpha, \quad \text{as} \quad 0 \leq \alpha \leq 1
```
and plot the convergence graphs.
