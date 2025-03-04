# MATH2089 NM Tutorial Week 2 - Python


#%% QUESTION 1

from numpy import exp, linspace
import matplotlib.pyplot as plt

# f(x) = e^{-x^2/2}


#%% QUESTION 2

a = 0.01; b = 2000; c = -0.001; x = -2.000000000005e5;
print('A1 =', a*x**2 + b*x + c)
print('A2 =', x*(a*x + b) + c)
print('A3 =', c + b*x + a*x**2)


#%% QUESTION 3

# numerically preferable quadratic formula
# b > 0 => diamond, heart
# b < 0 => spade, club

from numpy import sqrt

def quadsolve(a, b, c):
    D = b**2 - 4*a*c # D = discriminant
    if D > 0: # real roots
        r1 = 0   # <replace with appropriate code>
        r2 = 0   # <replace with appropriate code>
        
    elif D < 0: # complex conjugate roots
        sqrtmD = sqrt(-D)
        r1 = complex(-b,  sqrtmD) / (2*a)
        r2 = complex(-b, -sqrtmD) / (2*a)
    else: # D=0 so equal roots
        r1 = r2 = -b / (2*a)
    return r1, r2

a, b, c = 1, 3, 2              # case (i)
#a, b, c = 0.01, 2000, -0.001  # case (ii)

r1, r2 = quadsolve(a, b, c)
d1 = a*r1**2 + b*r1 + c
d2 = a*r2**2 + b*r2 + c
print(f'Coefficients: {a:g}, {b:g}, {c:g}')
print(f'       Roots: {r1:g}, {r2:g}')
print(f'      Values: {d1:.8e}, {d2:.8e}')   


#%% QUESTION 4

from numpy import array, inf
from numpy.linalg import norm

# vector norm
x = array([5, -4, 0, -6])


#%% QUESTION 5

from numpy import array, inf, sqrt, max
from numpy.linalg import norm, eig

# matrix norm
A = array([[3, -4, 2], [0, 4, -5], [2, -2, 3]])


#%% QUESTION 6

from numpy import array, eye, inf
from numpy.linalg import cond
from scipy.linalg import norm, eig, solve_triangular

A = array([[0, 3, -2], [-1, -4, 2], [5, 14, 26]])
B = array([[-11/8, -53/48, -1/48], [3/8, 5/48, 1/48], [1/16, 5/32, 1/32]])


# (a) B is the inverse of A

 
# (b) condition number measured in 1-norm
a1 = norm(A,1)   
b1 = norm(B,1)
c1 = cond(A,1)   
print('\na1  =', a1 )        
print('b1  =', b1 )        
print('c1  =', c1 ) 
print('c11 =', a1*b1 )

# (c) condition number measured in inf-norm
ainf = norm(A,inf)    
binf = norm(B,inf)
cinf = cond(A,inf)    
print('\nainf  =', ainf )       
print('binf  =', binf )        
print('cinf  =', cinf ) 
print('cinff =', ainf*binf )

# (d) condition number measured in 2-norm
a2 = norm(A,2)    
b2 = norm(B,2)
c2 = cond(A,2)  
print('\na2  =', a2 )       
print('b2  =', b2 )        
print('c2  =', c2 ) 
print('c22 =', a2*b2 )

print('\nEigenvalues of A =' )
print('Eigenvalues of B ='  )
print('      reciprical ='  )
print('           ratio ='  )


#%% QUESTION 7

from numpy import array
from scipy.linalg import norm, lu, solve_triangular

# (a) LU factorization PA = LU
A = array([[0, 3, -2], [-1, -4, 2], [5, 14, 26]])


# (b) check LU factorization


# (c) solve linear system using LU factorization
#     Ax = b => PAx = Pb => LUx = Pb => Ly = Pb
b = array([1, 2, 3])
