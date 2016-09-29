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
		
		
plt.xlim(1.0,w[len(w)-1])
plt.show()
