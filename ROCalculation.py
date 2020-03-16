# Es wird versucht, nach Vorbild von MA Werner eine Fitting Curve fuer das Ramberg-Osgood Modell zu erstellen, allerdings nicht in Matlab, sondern in Numpy.
# Um die Parameter a und n berechnen zu koennen, werden die Schubmodule benoetigt.
# Es werden zunaechst 2 Arrays erstellt, einmal werden die Temperaturen gespeichert, zum anderen die gemittelten Schubmodule aus den realen Messungen in N/mm^2.
#
# Material: IM7/8552

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Fitting functions for Ramber-Osgood
#-----------------------------------#
def func_ro_fit20(y, a, n):
	g12=4513.196315
	return (y+a*np.sign(y)*((np.abs(y)))**n)/g12

def func_ro_fit40(y, a, n):
	g12=4141.575438
	return y/g12+a*np.sign(y)*((np.abs(y))**n)/g12

def func_ro_fit60(y, a, n):
	g12=3770.993697
	return y/g12+a*np.sign(y)*((np.abs(y))**n)/g12

def func_ro_fit80(y, a, n):
	g12=3852.001111
	return y/g12+a*np.sign(y)*((np.abs(y))**n)/g12

def func_ro_fit100(y, a, n):
	g12=3234.411975
	return y/g12+a*np.sign(y)*((np.abs(y))**n)/g12

#----------------------------------------------------------#

# Sheae Modulus and Shear strength of IM7/8552 Real Tests
#-------------------------------------------------------#

g12=[4513.196315, 4141.575438, 3770.993697, 3852.001111, 3234.411975]
sl=[63.24905705, 58.31070307, 53.56012142, 52.44344755, 48.21212345]


# Ramberg-Osgood-Parameters of the CompDam-Model
#----------------------------------------------#
aPL_old=4.06e-9
nPL_old=5.4


# Filename
#--------#
name='80-1'
ending='.txt'
filename=name+ending


# Load xy Data from File; x: Strain, y: shear stress in MPa
#---------------------------------------------------------#
x, y = np.loadtxt(filename, delimiter='	', usecols=(0,1), unpack=True)
#y = np.loadtxt(filename, delimiter='	', usecols=(1), unpack=True)


# Plot of Dataset
plt.plot(x,y, '-b', label='Dataset '+name)

# Fitting RO curve
if name[0:2]=="20":
	popt, pcov = curve_fit(func_ro_fit20, y, x, maxfev=10000)
	plt.plot(func_ro_fit20(y, *popt), y, '-r', label='RO-Parameters, aPL = %.3E, nPL = %3.2f' % tuple(popt))
elif name[0:2]=="40":
	popt, pcov = curve_fit(func_ro_fit40, y, x, maxfev=10000)
	plt.plot(func_ro_fit40(y, *popt), y, '-r', label='RO-Parameters, aPL = %.3E, nPL = %3.2f' % tuple(popt))
elif name[0:2]=="60":
	popt, pcov = curve_fit(func_ro_fit60, y, x, maxfev=10000)
	plt.plot(func_ro_fit60(y, *popt), y, '-r', label='RO-Parameters, aPL = %.3E, nPL = %3.2f' % tuple(popt))
elif name[0:2]=="80":
	popt, pcov = curve_fit(func_ro_fit80, y, x, maxfev=10000)
	plt.plot(func_ro_fit80(y, *popt), y, '-r', label='RO-Parameters, aPL = %.3E, nPL = %3.2f' % tuple(popt))
elif name[0:3]=="100":
	popt, pcov = curve_fit(func_ro_fit100, y, x, maxfev=10000)
	plt.plot(func_ro_fit100(y, *popt), y, '-r', label='RO-Parameters, aPL = %.3E, nPL = %3.2f' % tuple(popt))
else:
	print "No Fitting Function used."




print popt[0], popt[1]

#plt.plot(func_ro_fit(y, *popt), y, '-r', label='RO-Parameters, aPL = %.3E, nPL = %3.2f' % tuple(popt))

plt.xlabel('$\\gamma_{12}$')
plt.ylabel('$\\tau_{12}$ [MPa]')
plt.legend()
plt.show()


