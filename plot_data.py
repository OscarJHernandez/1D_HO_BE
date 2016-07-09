import matplotlib
import numpy as np
from numpy.linalg import inv
import pylab
import math
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.integrate
from scipy import integrate
from scipy.optimize import curve_fit
from scipy import interpolate
from scipy import stats
import sys
import random

# constants
hbarc = 197.327

# We enumerate the Eigenvectors to be plotted
Nmax = 2
data_array =[]
Eig = []
for i in range(0,Nmax):
	data_i = 'EigenVec_00'+str(i)+'.txt'
	data_array.append(data_i)
	
# Read in the Eigenvalues of the Eigenvectors
data = open('EigenVectors.txt', 'r')
for line in data:
		line = line.strip()
		columns = line.split()
		Ei = float(columns[1])
		Eig.append(Ei)	

# First we will read in the u(r) and w(r) data for the desired potential.

for i in range(0,Nmax):	
	data = open(data_array[i], 'r')
	x=[]
	psi=[]
	
	for line in data:
		line = line.strip()
		columns = line.split()
		xi = float(columns[0])
		psi_xi = float(columns[1])
		psi_xi = psi_xi**2
		x.append(xi)
		psi.append(psi_xi)
	
	if i==0:
		c = 'b'
	else:
		c = "#%06x" % random.randint(0, 0xFFFFFF)
	
	plt.plot(x,psi,'r-',color=c)

@np.vectorize
def scattering(E,x):
	#0 = c1*Sin[E*L]+c2*Cos[E*L]
	#0 = -c1*Sin[E*L]+c2*Cos[E*L]
	# => F[E,L]  = A1*Sin[E*L] 
	# E*L = 2*Pi
	return math.Sin((E*x)/hbarc)


plt.show()	
