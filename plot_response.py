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
from scipy.integrate import quad
from scipy import integrate
from scipy.optimize import curve_fit
from scipy import interpolate
from scipy import stats
from scipy.optimize import curve_fit
import sys
import random

w=[]
Res =[]

data = open('response_1.txt', 'r')
for line in data:
		line = line.strip()
		columns = line.split()
		wi = float(columns[0])
		Ri = float(columns[1])
		w.append(wi)
		Res.append(Ri)
plt.plot(w,Res)

w=[]
Res=[]
data = open('response_2.txt', 'r')
for line in data:
		line = line.strip()
		columns = line.split()
		wi = float(columns[0])
		Ri = float(columns[1])
		w.append(wi)
		Res.append(Ri)
plt.plot(w,Res)	


w=[]
Res=[]
data = open('response_1_gauss.txt', 'r')
for line in data:
		line = line.strip()
		columns = line.split()
		wi = float(columns[0])
		Ri = float(columns[1])
		w.append(wi)
		Res.append(Ri)
plt.plot(w,Res)	

def normal(E,Nmax):
	mass =1.0
	hbarc = 197.327
	L = hbarc*math.sqrt(Nmax+2.0)
	p = math.sqrt(2.0*E*mass)
	Norm = L + hbarc*(1.0/(2.0*p))*math.cos((2.0*p/(hbarc))*L)
	#Norm = Norm*0.26827449478198306/hbarc
	return Norm		
		
print normal(  3.7962914575317711E-002, 30)/2.7552664629319596
		
plt.xlim(0.0,w[len(w)-1])
plt.show()
