# Implementation of numerical methods 


Hey, everyone. That is a repository for fixing progress on my university tasks. Hope you'll find something useful for yourself!

<br />

## The first task: Condition number. LU-decomposition. Inverse matrix

Write a program to calculate a condition number of a square matrix, using norms:
```math
||x||_{\alpha}, \qquad where \quad \alpha = 1, \infty, E
```
Calculation of an inverse matrix should include *LU*-decomposition.
Find the condition number of matrices of this pattern:
```math
a_{i,j} = min(i,j), \quad  i=1,2,3,...,n, \quad  j=1,2,3,...,n 
```


Make sure that the inverse matrix is tridiagonal,
its non-diagonal elements are 
```math
a^{-1}_{i, i \pm 1}=-1
``` 
and its diagonal elements are 
```math
a^{-1}_{i,i} =
\begin{cases}
2, i=1,2,...,n-1
\\
1, i=n
\end{cases}
```

Plot the dependence of the calculation time of the inverse matrix on the size
of the matrix.

<br />

### Description of the solution
So the `Task 1` directory contains 4 files:
- external.py
- cond_num.py
- results.csv
- one_minute.png

`external.py` is a place where i described 8 functions that run the whole process:
- generate_matrix - generates the matrix according to the task
- matrix_multiplication - multiplies two matrices using the basic algorithm
- lu_decompose - decomposes a matrix according to the Gaussian algorithm of LU-decomposition
- invert - returns the inverse matrix
- norm_1, norm_inf, norm_euclid - calculate 3 kinds of norms of a matrix
- condition_numbers - returns three condition numbers for each kind of norm

`cond_num.py` is the main program, that accomplishes the task. It gets the condition numbers and plots the dependency.

`_results.csv` is the file with those particular condition numbers.

`_minute.png` is the actual plot of the time and size relation.

> the execution time of a program is limited by the number of minutes mentioned at the beginning
> (the number could be float).

`Plato.png` is a screenshot of a plato during one-minute test.