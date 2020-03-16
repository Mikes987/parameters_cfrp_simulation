"""
I wrote this Python script for my master thesis in order to calculate sufficient parameters for simulating
temperature dependent progressive damage analysis of carbon fiber reinforced polymers (CFRP). For a successful
simulation, I need to find a way to transform the dataset and plot of shear stress analysis into the mathematical
Ramberg-Osgood-Equation.

This script is the first one of two in order to calculate sufficient Ramberg-Osgood Parameters by calculatin these
according to the operating temperature.

Further Information:
I wrote this script before I came into closer contact with libraries like numpy and pandas, so I mostly used the regular Python arrays and lists.
I am sure there are couple of methods to shorten the code and I'll check if I will be able to do further improvements.
"""

# Important Info
#--------------#

# Goal of this script:
# Try to make a mean curve of the shear tests according to their temperature. In
# order for this script to work, it is necessary, that txt-files with data follow a
# a certain pattern. Begin the name of txt files with their temperature. E.g. if you
# did shear tests at 20 Degrees Celsius then begin the txt-filenames with 20-1, 20-2,
# 20-x ... .txt

# The necessary Parameter Input will be the temperature.
# The Script works that way, that temperaturevalues have to be put in as integers.



# Import necessary packages
#-------------------------#

import numpy as np
import os
import glob
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit



# Fitting Functions for Ramberg-Osgood
# g12 is the mean of the measured shear moduli in MPa according to the operating temperature
#------------------------------------------------------------------------------------------#

def func_ro_fit20(y, a, n):
	g12=4436.89625250778
	return (y+a*np.sign(y)*((np.abs(y)))**n)/g12

def func_ro_fit40(y, a, n):
	g12=4157.30052490382
	return y/g12+a*np.sign(y)*((np.abs(y))**n)/g12

def func_ro_fit60(y, a, n):
	g12=3671.28353638837
	return y/g12+a*np.sign(y)*((np.abs(y))**n)/g12

def func_ro_fit80(y, a, n):
	g12=3365.50230552121
	return y/g12+a*np.sign(y)*((np.abs(y))**n)/g12

def func_ro_fit100(y, a, n):
	g12=3280.10977583825
	return y/g12+a*np.sign(y)*((np.abs(y))**n)/g12



# Main Function to load files and prepare datasets
#------------------------------------------------#

def ro(temp):

	# Sign for degrees celsius:
	cs = u'\xb0'

	# Write Log to follow data treatment.
	q=open("Logdata.txt", 'a+')

	stl= 'Plotting and RO-Fitting for ' + str(temp) + ' Degrees Celsius.\n'
	q.write(stl)
	q.write((len(stl)-1)*'-'+'\n')
	q.write("\n")

	# Store all txtdata from directory into an array to choose the specific ones.
	txtfiles = glob.glob('*.txt')

	# Arrays for x and y values, x is shear strain without any unit, y is shear stress in MPa. Most probably, they
	# will be 2d-arrays.
	x=[]
	y=[]
	l=[] # Storing the lengths of each array

	# Colorset for curves.
	colors=['b', 'g', 'r', 'c', 'm']

	# This script shall create 2 diagrams side by side.
	# Left:  Plot all datasets and its corresponding fitting curve.
	# Right: Plot the mean curve and its corresponding Ramberg-Osgood fit.

	fig, cvs = plt.subplots(1,2, sharex=True, sharey=True, figsize=(12,5))
	c1 = cvs[0]
	c2 = cvs[1]

	# Store x and y values of the essential datasets.
	filename=str(temp) + "-"
	for txt in txtfiles:
		if len(filename)==2 or len(filename)==3 or len(filename)==4:
			if txt[0:2]==filename or txt[0:3]==filename or txt[0:4]==filename:
				x.append(np.loadtxt(txt, delimiter='	', usecols=(0), unpack = True))
				y.append(np.loadtxt(txt, delimiter='	', usecols=(1), unpack = True))

	# If no Data was found, then the script will terminate, otherwise it will continue.
	if len(x)==0:
		q.write('No Data Found. Program will terminate here. Please check filename and input temperature.\n')
	else:
		q.write(str(len(x)) + ' Datasets were found and will be used for plotting and calculation.\n')

		# Store lengths of each dataset array, needed for checking to match array lengths.
		for i in range(len(x)):
			l.append(len(x[i]))
		q.write(str(l) + '\n')

		# Correct the values if there is a potential offset ==> Set the initial values of the shear stress as well as the shear strain to 0.
		q.write('Doing offset correction if needed ==> Initial values of Strain and Stress will be set to Zero.\n')
		for i in range(len(l)):
			xinit=x[i][0]
			for j in range(len(x[i])):
				x[i][j]-=xinit
		for i in range(len(l)):
			yinit=y[i][0]
			for j in range(len(y[i])):
				y[i][j]-=yinit

		# All datasets or subarrays for x and y need to have the same length. If they are not equal in the first place, a correction is needed.
		boolean=True
		i=1
		while i<len(x) and boolean==True:
			if len(x[0])!=len(x[i]):
				boolean=False
			else:
				i+=1
		if boolean==True:
			q.write("All datasets or arrays have the same length. No correction needed.\n")
		else:
			q.write("Lengths of datasets differ, correction is needed.\n")

			# Check which dataset or array is the smallest and which is the biggest for further calculation.
			minlen=min(l)
			minpos=l.index(min(l))
			maxlen=max(l)

			q.write('x[' + str(minpos) + '] is dataset with smallest amount of data with ' + str(minlen) + ' datapoints\n')

			# Since all values should be nearby: Delete all values in the bigger arrays, starting from position pos = minlen, if these values are higher than.
			i=0
			while i<len(x):
				if i!=minpos:
					if len(x[i])>minlen:
						boolean=x[i][minlen]>x[minpos][minlen-1]
						if boolean==True:
							q.write('Dataset of x[' + str(i) + '] at position ' + str(minlen) + ' is higher than final value of smallest dataset   ==> Removing all the higher values. ')

							# Using slice function.
							x[i]=x[i][np.s_[:minlen]]
							y[i]=y[i][np.s_[:minlen]]
							l[i]=len(x[i])
							maxlen=max(l)
						else:
							q.write('Dataset of x[' + str(i) + '] at position ' + str(minlen) + ' isnt higher than final value of smallest dataset ==> Search higher values to remove. ')
							delmin=minlen
							while delmin<len(x[i]) and x[minpos][minlen-1]>x[i][delmin]:
								delmin+=1
							if delmin<len(x[i]):
								x[i]=x[i][np.s_[:delmin]]
								y[i]=y[i][np.s_[:delmin]]
								l[i]=len(x[i])
						q.write(str(l) + '\n')
				i+=1

			maxlen=max(l)

			# If lengths of certain datasets are still bigger than the shortest one: check if
			# the dataset with more values is the minimum at that specific row ==> delete if yes.

			while maxlen!=minlen:
				for i in range(len(x)):
					if len(x[i])>minlen:
						q.write('Still adjustment needed for x[' + str(i) + ']. ')
						val=[]
						for j in range(len(x)):
							val.append(x[i][1])
						j=1
						while j<minlen and len(x[i])>minlen:
							if val[i]==min(val):
								x[i]=np.delete(x[i],j)
								y[i]=np.delete(y[i],j)
								val[i]=x[i][j]
								l[i]=len(x[i])
							else:
								j+=1
								if j<minlen:
									for k in range(len(x)):
										val[k]=x[k][j]
						q.write(str(l) + '\n')
				maxlen=max(l)


		# Now the lengths of all arrays and datasets should be equal.
		# If we can be sure that all arrays have the same length, plot and calculate
		# the mean values for x and y.

		xmean=np.mean((x),0)
		ymean=np.mean((y),0)
		
		# Plotting left diagram:
		# Remark: The thesis had been written in German, so titles and axis labels are in German, too.
		for i in range(len(x)):
			c1.plot(x[i],y[i], colors[i], label=filename+str(i+1))
		c1.plot(xmean, ymean, '-k', label='Mittel')

		c1.set_title('Kurvenschar bei ' + str(temp) + u'\xb0' +'C')
		c1.set_xlabel('$\\gamma_{12}$')
		c1.set_ylabel('$\\tau_{12}$ [MPa]')
		c1.legend()

		# Plotting right diagram and fitting Ramberg-Osgood values
		c2.plot(xmean, ymean, '-k', label='Mittel')

		if filename[0:2]=="20":
			popt, pcov = curve_fit(func_ro_fit20, ymean, xmean, maxfev=10000)
			c2.plot(func_ro_fit20(ymean, *popt), ymean, '-r', label='RO-Fit, $\\alpha_{\mathrm{PL}}$ = %.3E, $n_{\mathrm{PL}}$ = %3.2f' % tuple(popt))
		elif filename[0:2]=="40":
			popt, pcov = curve_fit(func_ro_fit40, ymean, xmean, maxfev=10000)
			c2.plot(func_ro_fit40(ymean, *popt), ymean, '-r', label='RO-Fit, $\\alpha_{\mathrm{PL}}$ = %.3E, $n_{\mathrm{PL}}$ = %3.2f' % tuple(popt))
		elif filename[0:2]=="60":
			popt, pcov = curve_fit(func_ro_fit60, ymean, xmean, maxfev=10000)
			c2.plot(func_ro_fit60(ymean, *popt), ymean, '-r', label='RO-Fit, $\\alpha_{\mathrm{PL}}$ = %.3E, $n_{\mathrm{PL}}$ = %3.2f' % tuple(popt))
		elif filename[0:2]=="80":
			popt, pcov = curve_fit(func_ro_fit80, ymean, xmean, maxfev=10000)
			c2.plot(func_ro_fit80(ymean, *popt), ymean, '-r', label='RO-Fit, $\\alpha_{\mathrm{PL}}$ = %.3E, $n_{\mathrm{PL}}$ = %3.2f' % tuple(popt))
		elif filename[0:3]=="100":
			popt, pcov = curve_fit(func_ro_fit100, ymean, xmean, maxfev=10000)
			c2.plot(func_ro_fit100(ymean, *popt), ymean, '-r', label='RO-Fit, $\\alpha_{\mathrm{PL}}$ = %.3E, $n_{\mathrm{PL}}$ = %3.2f' % tuple(popt))
		else:
			q.write("No Fitting Function used.")

		c2.set_xlabel('$\\gamma_{12}$')
		c2.set_ylabel('$\\tau_{12}$ [MPa]')
		c2.set_title('Ramberg-Osgood Fit durch gemittelte Kurve')
		c2.legend()

		# As a final step: Save temperature and RO-Parameters into a file and show and save plot.
		f=open("rop.txt", 'a+')
		f.write('%d %.6E %10.2f \n' %(temp, popt[0], popt[1]))
		f.close

		plt.savefig(str(temp)+'.png')
		#plt.show()

		end='|Done!|'
		q.write('\n' + str(len(end)*'-') + '\n' + end + '\n' + str(len(end)*'-') + '\n')
		q.write('\n')
		q.close
		


# Parameters to put/use manually

t=[20, 40, 60, 80, 100]

for temp in t:
	ro(temp)

q=open("Logdata.txt", 'a+')
q.write("Calculation Complete. Check Directory, you should see:\n")
for temp in t:
	q.write(' - ' + str(temp) + '.png\n')
q.write(' - rop.txt\n')
q.write(' - And of course this Log-File.\n')
q.close

print('Done. Check Logfile and further data in the directory of this python script.')
