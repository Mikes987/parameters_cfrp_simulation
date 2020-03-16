# Load necessary libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from numpy import nan as NaN

# Load rop.txt
df=pd.read_csv('rop.txt', header = None, delimiter=" ")
df=df.dropna(axis=1, how='all')

# Make Arrays for fitting curves
t=df[0].values
a=df[1].values
n=df[8].values

# Create functions for curve fit
def eatan(x, a, b, c, d, e):
	return a*np.exp(b*np.arctan(c*x+d)+e)

def atan(x, a, b, c, d):
	return a*np.arctan(b*x+c)+d

# do curvefit
popta, pcova = curve_fit(eatan, t, a, maxfev=1000000)
poptn, pcovn = curve_fit(atan, t, n, maxfev=1000000)

# Create massive number of datapoints
xfit = np.linspace(20,100,1000)
afit = eatan(xfit, *popta)
nfit = atan(xfit, *poptn)

# Create figure
fig, ax = plt.subplots(1, 2, figsize=(11,5))

ax[0].set_title('$\\alpha_\mathrm{PL}$')
ax[0].set_xlabel('$T$  [' + u'\xb0C]')
ax[0].set_ylabel('$\\alpha_\mathrm{PL}$')
ax[0].set_yscale('log')

ax[0].scatter(t, a)
ax[0].plot(xfit, afit)

ax[1].set_title('$n_\mathrm{PL}$')
ax[1].set_xlabel('$T$  [' + u'\xb0C]')
ax[1].set_ylabel('$n_\mathrm{PL}$')

ax[1].scatter(t, n)
ax[1].plot(xfit, nfit)

#plt.show()
