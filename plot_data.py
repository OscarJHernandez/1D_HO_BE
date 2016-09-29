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

# constants
hbarc = 197.327
xmax = 3000.0

# We enumerate the Eigenvectors to be plotted
Nmax = 9
data_array =[]
Eig = []

mass =1.0
hw = 2.0
N = 50
L = hbarc*math.sqrt(2.0*(N+2)/(mass*hw))


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
# We plot only the last one
#for i in range(0,Nmax):	
	#data = open(data_array[i], 'r')
	#x=[]
	#psi=[]
	
	#for line in data:
		#line = line.strip()
		#columns = line.split()
		#xi = float(columns[0])
		#psi_xi = float(columns[1])
		#psi_xi = psi_xi**2
		#x.append(xi)
		#psi.append(psi_xi)
	
	#if i==0:
		#c = 'b'
	#else:
		#c = "#%06x" % random.randint(0, 0xFFFFFF)
	
	#plt.plot(x,psi,'r-',color=c)

@np.vectorize
def scatteringSin(E,x):
    p = math.sqrt(2.0*E*mass)
    s = math.sin((p*x)/hbarc)
    s = s*s
    Norm = L - hbarc*(1.0/(2.0*p))*math.sin((2.0*p/(hbarc))*L)
    s = s/Norm
    return s
    
@np.vectorize
def scatteringCos(E,x):
    p = math.sqrt(2.0*E*mass)
    c = math.cos((p*x)/hbarc) 
    c=c*c   
    Norm = L + hbarc*(1.0/(2.0*p))*math.sin((2.0*p/(hbarc))*L)
    c = c/Norm
    return c 
    
@np.vectorize
def scatt_w_phase(E,delta,R,x):
    p = math.sqrt(2.0*E*mass)
    
    if x>0:
        c = math.cos(((p*x)/hbarc)+delta) # The scattering solution with phase shift
    else:
        c = math.cos(((p*x)/hbarc)-delta) # phase shift must change sign
    c=c*c   
    Norm = ((L-R)+(0.5*hbarc/p)*math.sin(2.0*(((p*L)/hbarc)+delta))-(0.5*hbarc/p)*math.sin(2.0*(((p*R)/hbarc)+delta))) 
    
    c = c/Norm
    return c 

# We want to calculate the phase shift!
def phaseshift(i,R):
    h0=0.0001
    p = math.sqrt(2.0*abs(Eig[i])*mass)
    
    data = open(data_array[i], 'r')
    
    x=[]
    psi=[]
	
    for line in data:
		line = line.strip()
		columns = line.split()
		xi = float(columns[0])
		psi_xi = float(columns[1])
		x.append(xi)
		psi.append(psi_xi)
    
    # Now do a spline fit of the function
    f = interpolate.interp1d(x, psi)
    fp = (f(R+h0)-f(R))/h0
    
    delta = -((R*p)/(hbarc))-math.atan((fp/f(R))*(hbarc/p))
    
    #print delta
    
    print 1.0*(p/(hbarc))-math.pi*0.5
    
    
    return delta,f,x
    
# \int_{-L}^{L} dx |psi(x)|^2 * x**p
def operPsi0_Norm(i,p):
    
    data = open(data_array[i], 'r')
    x=[]
    psi=[]
	
    for line in data:
		line = line.strip()
		columns = line.split()
		xi = float(columns[0])
		psi_xi = float(columns[1])
		x.append(xi)
		psi.append(psi_xi)
    
    # Now do a spline fit of the function
    f = interpolate.interp1d(x, psi)
    
    # The function to be integrated
    def func(x):
		return f(x)*f(x)*(x**p)
    
    
    s=quad(func, -L, L)[0]
    
    print s
    
    return s   
    


print operPsi0_Norm(0,4)
    
    
R = 100.0
n=0
delta,psi_spline,x = phaseshift(n,R)    
print delta
xvar = np.arange(x[0], x[len(x)-1], 0.1)
#plt.plot(xvar,scatteringCos(Eig[0],xvar),'-g')
#plt.plot(xvar,scatteringSin(Eig[1],xvar),'-g')

#plt.plot(xvar,scatteringSin(Eig[n],xvar),'--r')
plt.plot(xvar,psi_spline(xvar)**2,'-b')
#plt.plot(xvar,scatt_w_phase(Eig[n],delta,R,xvar),'--g')


#plt.plot(xvar,scatteringSin(Eig[2],xvar),'--r')
#plt.plot(xvar,scatteringCos(Eig[3],xvar),'--r')
plt.xlim(-xmax,xmax)
plt.show()	
