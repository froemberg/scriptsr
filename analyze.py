#!/usr/bin/python

import numpy as np
import os
import matplotlib.pyplot as plt

data_dir = "/group/ag_compstatphys/data/froemberg/enzyme_fbm/sim05q"
file_path = os.path.join(data_dir, "msd1.dat")
data = np.loadtxt(file_path, skiprows=1)
data = np.transpose(data)
t = data[0]
tc = t[2:]
msd = data[1]

qfile_path = os.path.join(data_dir, "mqd1.dat")
qdata = np.loadtxt(qfile_path, skiprows=1)
qdata = np.transpose(qdata)
tq = qdata[0]
tqc = t[2:]
mqd = qdata[1]

t1 = t[::5]
msd1 = msd[::5]
tc1 = t1[2:]
eqcorr = np.fabs(np.diff(msd,2)/(t[2]-t[1])**2)


#plt.xscale('log')
#plt.yscale('log')

#plt.plot(tq, mqd, 'r', label="mqd")
#plt.xlabel("time")
#plt.ylabel("mqd")
#plt.show()

#non-Gaussian parameter:
NGP = []
for j in range(len(t)):
	#lgngp = np.log(3*(mqd[j]-3*msd[j]*msd[j]))-np.log(5*(msd[j]))-np.log(msd[j])
	#ngp  =  np.exp(lgngp) - 1
	ngp = 3.0*(mqd[j]-3.0*msd[j]**2)/(5.0*msd[j]**2) - 1
	NGP.append(ngp)
	

plt.xscale('log')
plt.yscale('log')
plt.plot(tq, np.abs(NGP), 'r', label="NGP")
plt.xlabel("time")
plt.ylabel("NGP")
plt.show()

#plt.plot(tc, eqcorr, 'b', label="Incremental Autocorrelation")
#plt.xlabel("time")
#plt.ylabel("corr")
#plt.show()
