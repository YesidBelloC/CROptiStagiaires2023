# Car race along a track
# ----------------------
# An optimal control problem (OCP),
# solved with direct multiple-shooting.
#
# For more information see: http://labs.casadi.org/OCP
from casadi import *
import numpy as np
import math

N = 100 # number of control intervals

opti = Opti() # Optimization problem

# Input parameters
Rint  = 80*1e-3  #[Ohm]
Ccell = 3 * 3600  #[A.s]
Text  = 25+273.15 #[°K]
Pcell = 40

a1 = -0.5863
a2 = 21.9
a3 = 3.414
a4 = 0.1102
a5 = -0.1718
a6 = 0.008

b1 = 0.1369
b2 = -0.2518
b3 = 0.1609
b4 = -0.041
b5 = 0.0821

k1 = 5.896*(10**-10)
k2 = -18.75
k3 = 0.01388

d1 = -10260
d2 = 17230
d3 = -10130
d4 = 2340
d5 = 684.9

g1 = 8.913*(10**-15)
g2 = -32.23
g3 = 0.031
g4 = 0.007473

h1 = -154100
h2 = 204200
h3 = -4009
h4 = -81240
h5 = 22830
h6 = 7144

#Voc = 3.3
R0 = 0.0796
Rs = 0.0164
Cs = 710.5678
Rl = 0.0395
Cl = 5183.6054

m  = 0.05 #50gr
Cb = 3.56 #$Litium Specific Heat
h  = 0.026 #Air Convention Param
A  = 0.0652 # Hauteur de baterie
B  = 0.009  # Rayon batterie
S  = 2*np.pi*A*B

Vbmax = 3.3

# ---- decision variables ---------
X = opti.variable(4,N+1) # state trajectory
soc  = X[0,:]
vs   = X[1,:]
vl   = X[2,:]
# el   = X[3,:]
Tb   = X[3,:]
U = opti.variable(1,N)   # control trajectory (throttle)
T = opti.variable()      # final time
# T = 1000

# IntegrationLagrange = 0
# for l in range(N): # Lagrange Integration
#    Lagrange = -X[0,l]
#    IntegrationLagrange += Lagrange

# ---- objective          ---------
opti.minimize(Tb[-1])


# ---- dynamic constraints --------
# f = lambda x,u: vertcat(u/(Ccell),u/Cs-x[1]/(Rs*Cs),u/Cl-x[2]/(Rl*Cl),u*R0) # dx/dt = f(x,u)
f = lambda x,u: vertcat(u/(Ccell),u/Cs-x[1]/(Rs*Cs),u/Cl-x[2]/(Rl*Cl),1/(m*Cb)*(R0*u**2+h*S*Text-h*S*x[3])) # dx/dt = f(x,u)

dt = T/N # length of a control interval
for k in range(N): # loop over control intervals
   # Runge-Kutta 4 integration
   k1 = f(X[:,k],         U[:,k])
   k2 = f(X[:,k]+dt/2*k1, U[:,k])
   k3 = f(X[:,k]+dt/2*k2, U[:,k])
   k4 = f(X[:,k]+dt*k3,   U[:,k])
   x_next = X[:,k] + dt/6*(k1+2*k2+2*k3+k4)
   opti.subject_to(X[:,k+1]==x_next) # close the gaps

# ---- states constraints -----------
opti.subject_to(soc<=1)   # track soc limit
opti.subject_to(soc>=0)   # track soc limit
opti.subject_to(vl>=0)   # track vl limit
opti.subject_to(vs>=0)   # track vs limit
opti.subject_to(Tb<=6000+273.15)   # track vs limit
opti.subject_to(opti.bounded(-5,U,10)) # control is limited

# ---- boundary conditions --------
opti.subject_to(soc[0]==0)   # start at soc 0 ...
opti.subject_to(vs[0]==0) # ... from vs
opti.subject_to(vl[0]==0) # ... from vl
opti.subject_to(Tb[0]==25+273.150) # ... from Tb
opti.subject_to(soc[-1]==1)  # finish line at position

Voc = a1*2.7182**(-a2*soc) + a3 + a4*soc + a5*2.7182**(-a6/(1.01-soc))
opti.subject_to(Voc>=0)
# opti.subject_to(Voc<=3.5)

# ---- misc. constraints  ----------
opti.subject_to(T>=0) # Time must be positive

# ---- initial values for solver ---
# opti.set_initial(soc, 0)
opti.set_initial(T, 1000)

# ---- solve NLP              ------
opti.solver("ipopt") # set numerical backend
sol = opti.solve()   # actual solve

# ---- post-processing        ------
from pylab import step, figure, show, spy
import matplotlib.pyplot as plt

figure()
plt.plot(np.arange(0,sol.value(T),sol.value(T)/(N+1)),sol.value(Voc))
plt.xlabel('Temps [s]')
plt.ylabel('Voc')
plt.title('VoC state [-]')
plt.grid()

figure()
plt.plot(np.arange(0,sol.value(T),sol.value(T)/(N+1)),sol.value(soc))
plt.xlabel('Temps [s]')
plt.ylabel('soc')
plt.title('SoC state [-]')
plt.grid()

figure()
plt.plot(np.arange(0,sol.value(T),sol.value(T)/(N+1)),sol.value(vl))
plt.xlabel('Temps [s]')
plt.ylabel('vl')
plt.title('Vl state [v]')
plt.grid()

figure()
plt.plot(np.arange(0,sol.value(T),sol.value(T)/(N+1)),sol.value(vs))
plt.xlabel('Temps [s]')
plt.ylabel('vs')
plt.title('Vs state [v]')
plt.grid()

figure()
plt.plot(np.arange(0,sol.value(T),sol.value(T)/(N+1)),sol.value(Tb-273.15))
plt.xlabel('Temps [s]')
plt.ylabel('Tb')
plt.title('Temperature [C]')
plt.grid()

figure()
plt.plot(np.arange(0,sol.value(T),sol.value(T)/N),sol.value(U),'k')
plt.xlabel('Temps [s]')
plt.ylabel('Current')
plt.title('Input Signal [A]')
plt.grid()

# figure()
# spy(sol.value(jacobian(opti.g,opti.x)))
# figure()
# spy(sol.value(hessian(opti.f+dot(opti.lam_g,opti.g),opti.x)[0]))

show()
