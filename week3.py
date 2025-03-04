# MATH2089 NM Tutorial Week 3 - Python


#%% QUESTION 1

from numpy import array
from scipy.linalg import norm, cholesky, solve_triangular, eig

A = array([[2,-1,2], [-1,1,-1], [2,-1,3]])
b = array([2, 4, 1])
print('\nA =', A )
print('\nb =', b )

# (a) Is A symmetric?
ndiff=norm(A-A.T)
nsymm = ndiff/norm(A,2)
print('\ndiff =', ndiff )
print('\nsymm =', nsymm )

# (b) Does Cholesky factorization exist?
R = cholesky(A)
Rt = R.T
print('\nR =', R )

check = norm(A-R.T@R, 2)
print('\ncheck =', check) 

# (c) Is A positive definite?
eigvals,eigvecs = eig(A)

print('\nEigenvalues of A = ', eigvals)
print('\nMinimum eigenvalue of A = ', min(eigvals))

# (d) solve Ax = b

y=solve_triangular(R.T,b,lower=True)
x=solve_triangular(R,y,lower=False)

r=A@x-b
rnorm=norm(r,2) #abs error
relrnorm=norm(r,2)/norm(b,2) #rel error

print('\ny =',y ) 
print('\nx =', x ) 
print('\nr =', r )
print('\nrnorm =', relrnorm) 


#%% QUESTION 2

from scipy.io import mmread
from scipy.sparse import issparse
from scipy.sparse.linalg import norm as spnorm 
import matplotlib.pyplot as plt

# (a) load west0156 into sparse matrix A
A = mmread('west0156.mtx.gz').tocsc()

# (b) issparse, size, nnz
print('\nIs A sparse? =', issparse(A) )
print('\nsize of A =', A.shape )
print('\nnnz of A =', A.nnz )
print('\nA =', A )

# (c) full
F = A.todense()  # turns A into dense matrix format
print('\nIs F sparse? =', issparse(F) )
print('\nsize of F =', F.shape )
print('\nF =', F )

# (d) spy plot of A

# (e) Is A symmetric?
print('\nsymchkA ='  )

# (f) submatrix A_ij for i = 146,...,156 and j = 1,...,5 (index from 1)


#%% QUESTION 3

from numpy import prod, count_nonzero
from scipy.io import mmread
from scipy.linalg import cholesky, bandwidth
from scipy.sparse import eye as speye
from scipy.sparse.linalg import norm as spnorm 
from scipy.sparse.csgraph import reverse_cuthill_mckee
import matplotlib.pyplot as plt

# (a) B = A.T A
A = mmread('west0156.mtx.gz').tocsc()
B = A.T @ A

# (b) spy plot of B

plt.spy(B)

# (c) Is B symmetric?
print('\nsymchkB =' )

# (d) nonzero elements of B
print('\nnnz of B =' )

# (e) Is B postive definite?

# (f) Is C postive definite?



# (g) symmetric reverse Cuthill McKee permutation
p = reverse_cuthill_mckee(C)
C2 = C[p,:][:,p]    # moves the nonzero entries closer to diagonal
                    # Warning: C[p,p] does not work

# (h) spy plots for C, R, C2, R2

# (i) density of C, R, C2, R2
print('\ndensity of C =' )
print('\ndensity of R =' )
print('\ndensity of C2 =' )
print('\ndensity of R2 =' )