import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma
from matplotlib import rc

m1 = np.loadtxt("data/5d_0-0.5-0_s1e3.dat",unpack=True)
m2 = np.loadtxt("data/5d_0-1-0_s1e3.dat",unpack=True)
m3 = np.loadtxt("data/5ab/5a_000_s1e3.dat",unpack=True)
var = np.loadtxt("data/var.dat",unpack=True)
"""
m3 = np.loadtxt("data/5ab/5a_0.5-00_s1e3.dat",unpack=True)
m4 = np.loadtxt("data/5ab/5a_0.9-00_s1e3.dat",unpack=True)
"""

# ---------------- New color scheme -----------------
# These are the "Tableau 20" colors as RGB.    
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]    
  
# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.) 


#--------------------Configuration------------------------------
fig = plt.figure()
# Remove the plot frame lines. They are unnecessary chartjunk.    
ax = plt.subplot(111)  
ax.yaxis.grid() 
#ax.xaxis.grid()  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(False)    
ax.spines["right"].set_visible(False)    
ax.spines["left"].set_visible(False)   
plt.xticks(fontsize=10)    
plt.yticks(fontsize=10)    
# Remove the tick marks; they are unnecessary with the tick lines we just plotted.    
plt.tick_params(axis="both", which="both", bottom="off", top="off",    
                labelbottom="on", left="off", right="off", labelleft="on")

plt.grid(b=True, which='minor', alpha=0.2)

#-------------------------------------------------------------
N=len(m1)
"""
beta = 1/np.mean(m1)
print np.mean(m1)
omega1 = beta*np.exp(beta*m1)
omega2 = beta*np.exp(beta*m2)
omega3 = beta*np.exp(beta*m3)
omega4 = beta*np.exp(beta*m4)
plt.semilogy(m1,omega1,  color=tableau20[0],label=r'$\lambda = 0.0$')
plt.semilogy(m2,omega2,  color=tableau20[2],label=r'$\lambda = 0.25$')
plt.semilogy(m3,omega3,  color=tableau20[4],label=r'$\lambda = 0.5$')
plt.semilogy(m4,omega4,  color=tableau20[6],label=r'$\lambda = 0.9$')
plt.legend(loc="center right")
"""
plt.title(r'Wealth distribution -  N = '+str(N)+r' with 1e3 simulations ',fontsize=12)
plt.ylabel(r'Percent of total agents',fontsize=14)
plt.xlabel(r'wealth m ($m_0 = 100$)',fontsize=14)

#---------PLOT---------------------------------------

#Specify bin size
binsize = 10
N1=int(max(m1)/binsize)
N2=int(max(m2)/binsize)
N3=int(max(m3)/binsize)
#N4=int(max(m4)/binsize)

m = np.linspace(min(m3),max(m3),N3)
lamb = 0
n = 1+ (3*lamb)/(1-lamb)
x = m*n/np.mean(m3)
P = x**(n-1)*np.exp(-x)/(gamma(n))
#plt.loglog(m,m**-2.8,label='Trygve parameteriserer')
#plt.loglog(x,P, label='P(x)')

"""
data_hist1, binEdges = np.histogram(m1,bins=N1)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1]) #Center bin data
plt.loglog(bincenters, data_hist1/float(N), color=tableau20[0],label=r'$\alpha = 0.5$')

data_hist1, binEdges = np.histogram(m2,bins=N2)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1]) #Center bin data
plt.loglog(bincenters, data_hist1/float(N), color=tableau20[2],label=r'$\alpha = 1.0$')



data_hist1, binEdges = np.histogram(m3,bins=N3)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1]) #Center bin data
bincenters = bincenters * n / 100. #Rescale bincenters
dbins = bincenters[1]-bincenters[0] #Width of bins
plt.loglog(bincenters, data_hist1/(float(len(m3))*dbins), color=tableau20[0],label=r'$\alpha = 0$')
"""
plt.plot(var[0,:], var[1,:])
plt.legend(loc='upper right')
#plt.xlim(10,1000)
#plt.ylim(1e-2,0.2)



"""
# Simple histogram
plt.hist(m1,bins=N1, color=tableau20[0],label=r'$\lambda = 0.0$')
plt.hist(m2, bins=N2,color=tableau20[2],label=r'$\lambda = 0.25$')
plt.legend()
"""
#------------------------------------------------------------
#fig.savefig('5cLOGLOG_N500_varSavings.pdf', bbox_inches='tight',pad_inches=0.106)
plt.show()
