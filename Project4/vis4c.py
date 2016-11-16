import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy import special
#T (temperature), Cv (specific heat), X (susceptibility) 

c= 60 #Number of cores!
T,E,Cv,M,X,Mabs,Nacc = np.loadtxt("data/20_1_1_0.05_1e6_c60.dat",usecols=(0,1,2,3,4,5,6), unpack=True)
Ts,Es,Cvs,Ms,Xs,Mabss,Naccs = np.loadtxt("data/20_24_24_0.05_1e6_c60.dat",usecols=(0,1,2,3,4,5,6), unpack=True)
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
#-------------------------------------------------------

fig = plt.figure()
plt.subplot(111)
plt.plot(T, Nacc,color=tableau20[2],label='T=1')
plt.plot(T, Naccs,color=tableau20[4],label='T=2.4')

# Sjekk plott. Gaar fra 0 til 3500, 500 av gangen. 
# Gang med antall kjerner (32 i dette eksempelet).
plt.xticks(np.arange(0,len(T),2500),np.arange(0,len(T)*c,2500*c))
plt.title('Number of accepted states T=1 and T=2.4')
plt.ylabel('Number of accepted states',fontsize=10)

plt.xlabel("Number of cycles",fontsize=10)
#--------------------Configuration------------------------------
# Remove the plot frame lines. They are unnecessary chartjunk.    
ax = plt.subplot(111)  
ax.yaxis.grid()  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(False)    
ax.spines["right"].set_visible(False)    
ax.spines["left"].set_visible(False)   
plt.xticks(fontsize=10)    
plt.yticks(fontsize=10)    
plt.legend(loc="upper center")
#plt.xlim(0,100)
#plt.grid(b=True, which='major')
#plt.grid(b=True, which='minor', alpha=0.2)
#plt.minorticks_on()
# Remove the tick marks; they are unnecessary with the tick lines we just plotted.    
plt.tick_params(axis="both", which="both", bottom="off", top="off",    
                labelbottom="on", left="off", right="off", labelleft="on")  
#------------------------------------------------------------
#-------------------------------------------------------------
# if you want/need to save the plot in some format, you can use
# (bbox and pad make the figure to be tighten to the plot-box)
fig.savefig('20NaccofNT1and24.pdf', bbox_inches='tight',pad_inches=0.106)
plt.show()

