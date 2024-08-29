## Class: `Matrix`
### Constructor
```javascript
@constructor
@param {Number[][]} data - Two-dimensional array representing the matrix data.
```
### Methods:
**transpose()**
Time complexity: $O(i\cdot j)$ where $i$ is the number of rows and $j$ is the number of columns.

A matrix transpose is an operation on a matrix in which its rows and columns are swapped $A_{ij}^{T}=A_{ji}$. An example:
$\quad \\ A = \begin{vmatrix}1&2\\3&4\\5&6\end{vmatrix};\quad
A^{T} = \begin{vmatrix}1&3&5\\2&4&6\end{vmatrix};\\\quad$
```javascript
let A = new Matrix([[1,2],[3,4],[5,6]]);
let AT = A.transpose();
console.log(AT)
```
```bash
Matrix {data:[[1,3,5],[2,4,6]]}
```
Properties of the transposed matrix:

$(A^{T})^{T}=A;\\\quad
(k\cdot A)^{T}=k\cdot A^{T};\\\quad
(A+B)^{T}=A^{T}+B^{T};\\\quad
(A\cdot B)^{T}=B^{T}\cdot A^{T};\\\quad$
determinant  $A$ = determinant $A^{T};$
