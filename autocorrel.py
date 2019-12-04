import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.integrate

#Schreibe Autokorrelationszeiten und Temperatur in folgende Datei:

path="data20x20/"
autocorrelzeit = open(path+"correl.txt","w")
autocorrelzeit.write("Temperatur ; t_int\n" )

tint=0.0
#Nutze folgende Funktion für den fit:
def func(x, a, t):
 	  return a*np.exp(-x/t)

for i in range(0,100):
	#einlesen der Daten
	t=i*0.05
	data = np.genfromtxt(path+"Wolff"+str("{:1.2f}".format(t))+"0000.txt",delimiter=";",skip_header=1,usecols=0,max_rows=1000)

	mag= np.mean(data)
	gamma= np.zeros_like(data)

	#print(mag)
	#fitten der exponentialfunktion

	#berechnen der  exponentiellen Autokorrelationszeit t für jede Temperatur
	
	x=np.arange(len(gamma))	
	
	for n in range(0,100):
		for i in range(0,len(data)-n):
		#print(gamma[n])
			gamma[n]=gamma[n]+((data[n+i]-mag)*(data[i]-mag))
		gamma[n]=gamma[n]/(len(data)-n)
	
	popt, pcov = curve_fit(func, x, abs(gamma))
	#yy=func(x,*popt)
	f= lambda x:popt[0]*np.exp(-x/(popt[1]))
	tint=scipy.integrate.quad(f,0,len(data))[0]
	#in Datei schreiben
	autocorrelzeit.write(str("{:1.2f}".format(t))+";"+str(tint) +"\n")
	#print(tint)'''


#print(gamma[n])
##Der Versuch eines exponential fits





	#print(popt)
##Plotten der Autokorrelationsfkt	

'''plt.plot(x,yy)
plt.plot(x,abs(gamma))
plt.xlim(0,100)
plt.show()'''
 