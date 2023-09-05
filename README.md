## The first task

Write a program to calculate a condition number of a square matrix, using norms: *1, inf, E*.
Calculation of an inverse matrix should include *LU*-decomposition.
Find the condition number of matrices of this pattern:
```math
a_{i,j} = min(i,j)   i=1,2,3,...,n   j=1,2,3,...,n 
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
\end{cases}```