from numpy import array
from scipy.linalg import cholesky, solve_triangular

A = array([[1,2,-1], [2,8,4], [-1,4,11]]) 
b = array([1,4,3]) 
R = cholesky(A) # A = R.T@R 
y = solve_triangular(R.T, b, lower=True) 
x = solve_triangular(R, y) 
print("y = ", y) 
print("x = ", x)

