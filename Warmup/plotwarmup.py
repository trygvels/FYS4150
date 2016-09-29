import matplotlib.pyplot as plt
import numpy as np

x_exact = 1./3.
N=1000
f2cs = np.zeros(N)
f3cs = np.zeros(N)
h = np.zeros(N)

f = open('data.txt','r')
k = 0
for line in f:
	f2cs[k] = line.split()[0]
	f3cs[k] = line.split()[1]
	h[k] = line.split()[2]
	k += 1

# Old code

h_best = h[np.argmin(abs(x_exact-f2cs))] # Step length corresponding to low err
f2c_best = f2cs[np.argmin(abs(x_exact-f2cs))]
f3c_best = f3cs[np.argmin(abs(x_exact-f2cs))]

print "Exact: %f" % (x_exact)
print "2c: %f with error: %g" % (f2c_best, min(abs(x_exact-f2cs)))
print "3c: %f with error: %g" % (f3c_best, min(abs(x_exact-f3cs)))
print "Best h: %g" % (h_best)

eps2c = np.log10(abs((f2cs-x_exact)/(x_exact)))
eps3c = np.log10(abs((f3cs-x_exact)/(x_exact)))
plt.plot(h, eps2c, h, eps3c)
plt.legend(['2c','3c'])
plt.show()
