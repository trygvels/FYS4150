import sys
import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18

f1 = open('data2Iw1.txt','r')
f2 = open('data2Ow1.txt','r') 

Omega_r = float(f.readline().split()[1])
rho_min = float(f.readline().split()[1])
rho_max = float(f.readline().split()[1])
n 	= int(f.readline().split()[1])
psi = np.zeros(n)
rho = np.linspace(rho_min,rho_max,n)


allLines = f.readlines()
allLines.pop(0)
allLines.pop(-1)

i=0
for line in allLines:
	psi[i] = float(line)
	i+=1

f.close()

plt.plot(rho,abs(psi)**2,'-b')
plt.title(r'$\omega_i=$'+str(Omega_r), fontsize=22)
plt.xlabel(r'$\rho_i$', fontsize=20)
plt.ylabel(r'$|u_i|^2$', fontsize=20)
plt.legend(['Numerical solution'],fontsize=18)
plt.show()

