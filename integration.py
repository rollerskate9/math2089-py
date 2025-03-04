import scipy.integrate as integrate

# define the function
def f(x):
    return (x)

# attatch the integration result to the variable result
result = integrate.quad(f, 0, 1)
print(result)