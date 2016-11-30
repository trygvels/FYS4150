import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

#m1 = np.loadtxt("data/5a_N5e4_trans1e7.dat", unpack=True)
m1 = np.loadtxt("build-Project5-Desktop-Debug/mony.dat",unpack=True)
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
beta = 1/np.mean(m1)
omega = beta*np.exp(beta*m1)

plt.title(r'Money distribution with N = '+str(N)+r' agents with savings',fontsize=12)

plt.ylabel(r'log($\omega_m$)',fontsize=14)
plt.xlabel(r'Money m',fontsize=14)
plt.hist(m1,bins=50, color=tableau20[4])
#plt.semilogy(m1,omega, color=tableau20[4]) #Shows straight line; omega is exponential

#fig.savefig('hist5c.pdf', bbox_inches='tight',pad_inches=0.106)
plt.show()
