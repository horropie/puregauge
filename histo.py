import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.integrate

path="data/"

for i in [1,20,50,100]:
	#einlesen der Daten

	t=i*0.05
	
	data = np.genfromtxt(path+"Metropolis"+str("{:1.2f}".format(t))+"0000.txt",delimiter=";",skip_header=1,usecols=0,max_rows=1000)
	plt.hist(data,label=str(t))
plt.legend()
plt.show()