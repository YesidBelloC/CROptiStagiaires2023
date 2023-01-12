# Car race along a track
# ----------------------
# An optimal control problem (OCP),
# solved with direct multiple-shooting.
#
# For more information see: http://labs.casadi.org/OCP
from casadi import *
import numpy as np
import math

N = 20 # number of control intervals

opti = Opti() # Optimization problem

# Input parameters

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


Qcell = 2.6 * 3600  #[A.s]
Vbmax = 4.2
Vbnom = 3.7
Ifast = 1.3
#Voc = 3.3

R0 = 0.0796
Rs = 0.0173
Cs = 745.2904
Rl = 0.0521
Cl = 5786.49

m  = 0.05 #50gr
Cb = 1000 #Lithium Specific Heat
h  = 0.026 #Air Convection Param

A  = 0.0652 # Hauteur de batterie
B  = 0.009  # Rayon batterie
S  = 2*np.pi*A*B

Text  = 25+273.15 #[°K]


# ---- decision variables ---------
X = opti.variable(4,N+1) # state trajectory
soc  = X[0,:]
vs   = X[1,:]
vl   = X[2,:]
Tb   = X[3,:]

U = opti.variable(1,N)   # control trajectory (throttle)

Vcell = opti.variable()

T = opti.variable()      # final time
# T = 180*60
# t = np.linspace(0,T,N+1)

#Lagrange Integration
IntegrationLagrange = 0
for l in range(N): 
    #IntegrationLagrange += (X[0,l]-1)**2
    #IntegrationLagrange += -X[0,l]
    IntegrationLagrange += (X[3,l])
    
    
# ---- objective          ---------


J = IntegrationLagrange
# J = -soc[-1]
# J = T

opti.minimize(J)
#opti.minimize(Tb[-1])
#opti.minimize((T-180*60)**2)
#opti.minimize(T)
#opti.minimize(-soc[-1])
#opti.minimize(IntegrationLagrange)


# ---- dynamic constraints --------
#f = lambda x,u: vertcat(u/(Qcell),u/Cs-x[1]/(Rs*Cs),u/Cl-x[2]/(Rl*Cl)) #dx/dt = f(x,u)
f = lambda x,u: vertcat(u/(Qcell),\
                        u/Cs-x[1]/(Rs*Cs),\
                        u/Cl-x[2]/(Rl*Cl),\
                        (R0*u**2 + (x[1]**2)/Rs + (x[2]**2)/Rl - h*S*(x[3]-Text) ) / (m*Cb)) # dx/dt = f(x,u)


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
opti.subject_to(Tb<=60+273.15)   # track Tb limit
opti.subject_to(opti.bounded(0,U,Ifast)) # control is limited
#opti.subject_to(U==Ifast)

# ---- boundary conditions --------
opti.subject_to(soc[0]==0)   # start at soc 0 ...
opti.subject_to(vs[0]==0) # ... from vs
opti.subject_to(vl[0]==0) # ... from vl
opti.subject_to(Tb[0]==25+273.15) # ... from Tb
#opti.subject_to(U[-1]==0.052)
opti.subject_to(Vcell[-1]==4.2)
opti.subject_to(soc[-1]==1)  # finish line at position


Voc = 0.55 + a1*np.exp(-a2*soc) + a3 + a4*soc + a5*np.exp(-a6/(1-soc))
#Voc = 4.0
#opti.subject_to(Voc>=0)
#opti.subject_to(Voc<=100)


#Vcell = Voc[0:N] - vs[0:N] - vl[0:N] - U*R0
Vcell = Voc[0:N] + vs[0:N] + vl[0:N] + U*R0
#opti.subject_to(Vcell>=0)
#opti.subject_to(Vcell<=4.2)

#opti.subject_to(opti.bounded(0,Vcell,4.2))

# ---- misc. constraints  ----------
opti.subject_to(T>=0) # Time must be positive
#opti.subject_to(T==180*60)

# ---- initial values for solver ---
# opti.set_initial(soc, 0)
opti.set_initial(T, 1000)
#opti.set_initial(T, 180*60)

# ---- solve NLP              ------
opti.solver("ipopt") # set numerical backend
sol = opti.solve()   # actual solve

# ---- post-processing        ------
from pylab import step, figure, show, spy
import matplotlib.pyplot as plt

t = np.linspace(0,sol.value(T),N+1)

# figure()
# plt.plot(t/60,sol.value(Voc))
# plt.xlabel('Temps [min]')
# plt.ylabel('Voc')
# plt.title('VoC state [-]')
# plt.grid()

figure()
plt.plot(t/60,sol.value(soc))
plt.xlabel('Temps [min]')
plt.ylabel('soc')
plt.title('SoC state [-]')
plt.grid()

# figure()
# plt.plot(t/60,sol.value(vl))
# plt.xlabel('Temps [min]')
# plt.ylabel('vl')
# plt.title('Vl state [v]')
# plt.grid()

# figure()
# plt.plot(t/60,sol.value(vs))
# plt.xlabel('Temps [min]')
# plt.ylabel('vs')
# plt.title('Vs state [v]')
# plt.grid()

figure()
plt.plot(t/60,sol.value(Tb-273.15))
plt.xlabel('Temps [min]')
plt.ylabel('Tb')
plt.title('Temperature [C]')
plt.grid()

# figure()
# plt.plot(np.linspace(0,sol.value(T),N)/60,sol.value(U),'k')
# plt.xlabel('Temps [min]')
# plt.ylabel('Current')
# plt.title('Input Signal [A]')
# plt.grid()

# figure()
# plt.plot(np.arange(0,sol.value(T),sol.value(T)/(N)),sol.value(Vcell),'k')
# plt.xlabel('Temps [s]')
# plt.ylabel('Cell voltage')
# plt.title('Input Signal [A]')
# plt.grid()

# figure()
# plt.plot(np.arange(0,sol.value(T),sol.value(T)/N)/60,sol.value(U)*sol.value(Vcell),'k')
# plt.xlabel('Temps [min]')
# plt.ylabel('Power supplied')
# plt.grid()

fig, ax1 = plt.subplots() 
ax1.set_xlabel('Time [min]') 
ax1.set_ylabel('Courant', color = 'red') 
ax1.plot(np.arange(0,sol.value(T),sol.value(T)/N)/60,sol.value(U), color = 'red') 
ax1.tick_params(axis ='y', labelcolor = 'red') 
ax2 = ax1.twinx()   
ax2.set_ylabel('Tension', color = 'blue') 
#ax2.set_ylim([-500,500])
ax2.plot(np.arange(0,sol.value(T),sol.value(T)/N)/60,sol.value(Vcell), color = 'blue') 
ax2.tick_params(axis ='y', labelcolor = 'blue') 
plt.show()



# figure()
# spy(sol.value(jacobian(opti.g,opti.x)))
# figure()
# spy(sol.value(hessian(opti.f+dot(opti.lam_g,opti.g),opti.x)[0]))

show()
