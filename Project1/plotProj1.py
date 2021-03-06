import sys
import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18

name = str(sys.argv[1])
f = open('dataproj'+name+'.txt','r')
n = int(f.readline())
x = np.zeros(n+2)
u = np.zeros(n+2)
v = np.zeros(n+2)

k = 0
next(f)
for line in f:

	x[k] = line.split()[0]
	u[k] = line.split()[1]
	v[k] = line.split()[2]
	k += 1

tall = str(n)
plt.plot(x,u,'-b',x,v,'--r')
plt.title('n ='+tall, fontsize=22)
plt.xlabel(r'$x_i$', fontsize=20)
plt.ylabel(r'$u(x_i) \quad \& \quad v_i$', fontsize=20)
plt.legend(['Exact solution','Numerical solution'],fontsize=18)
plt.show()

