# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 14:38:15 2022

@author: bvenot
"""

import numpy as np
import math
import matplotlib.pyplot as plt 
from pylab import step, figure, show, spy

# Input parameters
Rint  = 0  #[Ohm]
Ccell = 3 * 3600  #[A.s]
Text  = 25 #[°K]
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

# Time
Nech = 1800
tmax = 1800
ts = np.linspace(0, tmax, Nech) #[s]
dt = ts[1]

# Control variable
# I = np.zeros(100)
# Imax = 10
# stage1 = 4*np.ones(25)
# stage2 = 3*np.ones(25)
# stage3 = 1*np.ones(25)
# stage4 = 3*np.ones(25)
# stage12 = np.concatenate((stage1,stage2))
# stage23 = np.concatenate((stage12,stage3))
# I = -np.concatenate((stage23,stage4))

Imax = 10
I = -Imax*np.ones(Nech)
Imin = 0.5
Ncte1 = 600
Ndec = 850
Icte1 = Imax*np.ones(Ncte1)
Idec  = -((Imax-Imin)/Ndec)*(ts[0:Ndec]) + Imax
Icte2 = 0.5*np.ones(Nech-Ncte1-Ndec)
I = np.concatenate((Icte1, Idec))
I = -np.concatenate((I, Icte2))

I=-I
# fig
figure()
plt.plot(I)
plt.ylabel("Courrant profile")
plt.grid()
plt.show()

# State variables
SoC0   = SoC   = 0 #[s.u]
Vcell0 = Vcell = 0 #[V]

Voc = 0
Vs = 0
Vl = 0
Vcell = 0
R0    = 0  #[Ohm]
Rs    = 0  #[Ohm]
Cs    = 0  #[F]
Rl    = 0  #[Ohm]
Cl    = 0  #[F]

SoCs = np.zeros(Nech)
Voc_vec = np.zeros(Nech)
Vcell_vec = np.zeros(Nech)
R0_vec = np.zeros(Nech)
Rs_vec = np.zeros(Nech)
Rl_vec = np.zeros(Nech)
Cs_vec = np.zeros(Nech)
Cl_vec = np.zeros(Nech)
tau_s_vec = np.zeros(Nech)
tau_l_vec = np.zeros(Nech)

for t in range(len(ts)):  
    Voc = a1*math.exp(-a2*SoC) + a3 + a4*SoC + a5*math.exp(-a6/(1-SoC))   
    # R0  = b1*(SoC**4) + b2*(SoC**3) + b3*(SoC**2) +b4*SoC + b5
    # Rs  = k1*math.exp(-k2*SoC) + k3
    # Cs  = d1*(SoC**4) + d2*(SoC**3) + d3*(SoC**2) + d4*SoC + d5
    # Rl  = g1*math.exp(-g2*SoC) + g3 + g4*SoC
    # Cl  = h1*(SoC**5) + h2*(SoC**4) + h3*(SoC**3) + h4*(SoC**2) + h5*SoC + h6
    
    #Voc = 3.3
    R0 = 0.0796
    Rs = 0.0164
    Cs = 710.5678
    Rl = 0.0395
    Cl = 5183.6054
    
    dSoC_dt = I[t]/Ccell
    dVs_dt = I[t]/Cs - Vs/(Rs*Cs)
    dVl_dt = I[t]/Cl - Vl/(Rl*Cl)
    
    SoC += dSoC_dt * dt
    Vs  += dVs_dt * dt
    Vl  += dVl_dt * dt
    
    tau_s = Rs*Cs
    tau_l = Rl*Cl

    Vcell = Voc - R0*I[t] - Vs - Vl
    
    SoCs[t]      = SoC
    Voc_vec[t]   = Voc
    R0_vec[t]    = R0
    Rs_vec[t]    = Rs
    Rl_vec[t]    = Rl
    Cs_vec[t]    = Cs
    Cl_vec[t]    = Cl
    Vcell_vec[t] = Vcell
    tau_s_vec[t] = tau_s
    tau_l_vec[t] = tau_l

print("mean R0 = "+ str(np.mean(R0_vec))+" Ohm")
print("mean Rs = "+ str(np.mean(Rs_vec))+" Ohm")
print("mean Cs = "+ str(np.mean(Cs_vec))+" F")
print("mean Rl = "+ str(np.mean(Rl_vec))+" Ohm")
print("mean Cl = "+ str(np.mean(Cl_vec))+" F")
print("mean tau_s = "+ str(np.mean(tau_s_vec))+" s")
print("mean tau_l = "+ str(np.mean(tau_l_vec))+" s")

# fig
figure()
plt.plot(ts, SoCs)
plt.ylabel("SoC")
plt.xlabel("Time (s)")
plt.grid()
plt.show()

# fig
figure()
plt.plot(ts, Vcell_vec*I)
plt.ylabel("Puissance absorbée (W)")
plt.xlabel("Time (s)")
plt.grid()
plt.show()

# fig
figure()
plt.plot(ts, Vcell_vec)
plt.ylabel("Tension Batterie (v)")
plt.xlabel("Time (s)")
plt.grid()
plt.show()

# fig
figure()
plt.plot(ts, R0_vec)
plt.ylabel("Résistance interne (Ohm)")
plt.xlabel("Time (s)")
plt.grid()
plt.show()

# fig
figure()
plt.plot(ts, Voc_vec)
plt.ylabel("Open Circuit Voltage (V)")
plt.xlabel("Time (s)")
plt.grid()
plt.show()

# fig
figure()
plt.plot(SoCs, Voc_vec)
plt.ylabel("VoC")
plt.xlabel(" (s)")
plt.grid()
plt.show()

# fig
figure()
ax1 = plt.subplots() 
ax1.set_xlabel('Time (s)') 
ax1.set_ylabel('Courant', color = 'red') 
ax1.plot(ts, np.abs(I), color = 'red') 
ax1.tick_params(axis ='y', labelcolor = 'red') 
ax2 = ax1.twinx()   
ax2.set_ylabel('Tension', color = 'blue') 
#ax2.set_ylim([-500,500])
ax2.plot(ts, Vcell_vec, color = 'blue') 
ax2.tick_params(axis ='y', labelcolor = 'blue') 
ax2.grid()
plt.show()

