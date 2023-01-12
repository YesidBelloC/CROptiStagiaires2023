import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.signal

# Input parameters
Rint  = 0  #[Ohm]
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

I0 = I = 0

SoC0 = SoC = 0
vs0 = vs = 0
vl0 = vl = 0

Voc = a1*2.7182**(-a2*SoC) + a3 + a4*SoC + a5*2.7182**(-a6/(1.01-SoC))
Vcell = Voc - R0*I0 - vs - vl

ts = np.linspace(0, 1000, 1000)
dt = ts[1]

Kc = 10
tau_i = 500

Gc = scipy.signal.lti([Kc*tau_i, Kc], [tau_i, 0]).to_ss()

vobj = Vbmax

SoCs = []
vss = []
vls = []

vcells = []
vocs = []
Is = []

Tc = np.zeros([Gc.A.shape[0], 1])
for t in ts:

    e =  - Vcell
    dTcdt = Gc.A.dot(Tc) + Gc.B.dot(e)
    yc = Gc.C.dot(Tc) + Gc.D.dot(e)

    I = I0 + yc[0,0]  # I0 is the controller bias

    dSoCdt = I/(Ccell)
    dvsdt  = I/Cs-vs/(Rs*Cs)
    dvldt  = I/Cl-vl/(Rl*Cl)


    # f = lambda x,u: vertcat(u/(Ccell),u/Cs-x[1]/(Rs*Cs),u/Cl-x[2]/(Rl*Cl),1/(m*Cb)*(R0*u**2+h*S*Text-h*S*x[3])) # dx/dt = f(x,u)

    SoC += dSoCdt*dt
    vs += dvsdt*dt
    vl += dvldt*dt

    Tc += dTcdt*dt

    Voc = a1*2.7182**(-a2*SoC) + a3 + a4*SoC + a5*2.7182**(-a6/(1.01-SoC))
    Vcell = Voc - R0*I - vs - vl

    SoCs.append(SoC)
    vss.append(vs)
    vls.append(vl)

    vcells.append(Vcell)
    vocs.append(Voc)

    Is.append(I)


plt.figure(figsize=(10,4))
plt.plot(ts, SoCs)
#plt.fill_between(DistCum,theta_vals_int,alpha=0.1)
plt.ylabel("SoC (-)")
plt.xlabel("Temps (s)")
plt.grid()
plt.show()

plt.figure(figsize=(10,4))
plt.plot(ts, vss)
#plt.fill_between(DistCum,theta_vals_int,alpha=0.1)
plt.ylabel("Tension RC s (V)")
plt.xlabel("Temps (s)")
plt.grid()
plt.show()

plt.figure(figsize=(10,4))
plt.plot(ts, vls)
#plt.fill_between(DistCum,theta_vals_int,alpha=0.1)
plt.ylabel("Tension RC l (V)")
plt.xlabel("Temps (s)")
plt.grid()
plt.show()

plt.figure(figsize=(10,4))
plt.plot(ts, Is)
#plt.fill_between(DistCum,theta_vals_int,alpha=0.1)
plt.ylabel("Courrent input (A)")
plt.xlabel("Temps (s)")
plt.grid()
plt.show()