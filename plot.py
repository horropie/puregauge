import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.integrate


O=[]
err=[]
path="data/"
data2 = np.genfromtxt(path+"correl_poly.txt",delimiter=";",skip_header=1,usecols=1)
for i in range(1,100):
	#einlesen der Daten

	t=i*0.05
	
	data = np.genfromtxt(path+"Metropolis"+str("{:1.2f}".format(t))+"0000.txt",delimiter=";",skip_header=1,usecols=1,max_rows=1000)
	
	
	O.append(abs(np.mean(data))*1000)
	err.append(np.std(data)*1000*data2[i-1])

x=np.arange(0.05,5,0.05)
plt.figure()
plt.subplot(211)
plt.errorbar(x,O,yerr=err)

plt.ylabel(r"Polyakov Line $\cdot 10^3$  ")
plt.subplot(212)

plt.xlabel(r"$\beta$")
plt.ylabel(r"$\tau_{int}\cdot 10^4$")
plt.plot(x,data2*10000,color=(0,0.18,0.36))
#Plot Autocorrel

#data = np.genfromtxt(path+"correl_O1.txt",delimiter=";",skip_header=1,usecols=1)
#x=np.arange(99)
#plt.plot(x,data)
plt.xlim(0,5)
plt.show()
 

