from math import *
import numpy as np
from matplotlib.pyplot import plot, legend, show
N = 10000
x = sqrt(2)
x_exact = 1./3

h = np.linspace(10**-20,0.5,N)
f2cs = np.zeros(N)
f3cs = np.zeros(N)

def f(x):
	return np.arctan(x)
def f2c(x,h):
	return ( f(x+h)-f(x) ) /h #+ O(h)
def f3c(x,h):
	return ( f(x+h)-f(x-h) )/(2*h) #+ O(h**2)

for i in range(N):
	f2cs[i] = f2c(x,h[i])
	f3cs[i] = f3c(x,h[i])

h_best = h[min_element(abs(x_exact-f2cs))] # Step length corresponding to low err
f2c_best = f2c(x,h_best) 
f3c_best = f3c(x,h_best)

print "Exact: %f" % (x_exact)
print "2c: %f with error: %g" % (f2c_best, min(abs(x_exact-f2cs)))
print "3c: %f with error: %g" % (f3c_best, min(abs(x_exact-f3cs)))
print "Best h: %g" % (h_best)

eps2c = np.log10(abs((f2cs-x_exact)/(x_exact)))
eps3c = np.log10(abs((f3cs-x_exact)/(x_exact)))
plot(h, eps2c, h, eps3c)
legend(['2c','3c'])
show()
