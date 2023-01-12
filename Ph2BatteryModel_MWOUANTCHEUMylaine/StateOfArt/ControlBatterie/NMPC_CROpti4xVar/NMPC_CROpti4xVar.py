import numpy as np
import matplotlib.pyplot as plt
from casadi import *
from casadi.tools import *
import pdb
import sys
sys.path.append('C:/Users/crybelloceferin/Documents/CROpti/ControlBatterie/NMPC_CROpti4xVar')
import do_mpc

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import time

from modelCROpti import modelCROpti
from mpcCROpti import mpcCROpti
from simulatorCROpti import simulatorCROpti

#Information:  https://www.do-mpc.com/en/latest/index.html

""" User settings: """
show_animation = False
store_results = False
sample_Time = 1.0

"""
Get configured do-mpc modules:
"""

model = modelCROpti()
mpc = mpcCROpti(model)
simulator = simulatorCROpti(model)
estimator = do_mpc.estimator.StateFeedback(model)



"""
Set initial state
"""
time_list = []
time_start = time.process_time()

SoC_0 = 0.0
Vs_0 = 0.0
Vl_0 = 0.0
Tb_0 = 25+273.15
x0 = np.array([SoC_0, Vs_0, Vl_0, Tb_0])

socs = []
vss = []
vls = []
vss = []
tbs = []
us = [0]

socs.append(float(x0[0]))
vss.append(float(x0[1]))
vls.append(float(x0[2]))
tbs.append(float(x0[3]))

mpc.x0 = x0
simulator.x0 = x0
estimator.x0 = x0

mpc.set_initial_guess()

"""
Run MPC main loop:
"""

for k in range(int(180*60/5)):
    tic = time.process_time()
    u0 = mpc.make_step(x0)
    toc = time.process_time()
    # if float(u0)<0:
    #     p_num = simulator.get_p_template()
    #     p_num['Eff'] = 0.2
    #     def p_fun(t_now):
    #         return p_num
    #     simulator.set_p_fun(p_fun)
    # else:
    #     p_num = simulator.get_p_template()
    #     p_num['Eff'] = 1/0.8
    #     def p_fun(t_now):
    #         return p_num
    #     simulator.set_p_fun(p_fun)
    y_next = simulator.make_step(u0)
    x0 = estimator.make_step(y_next)
    # if float(x0[0])>500:
    #     tvp_template = mpc.get_tvp_template()
    #     def tvp_fun(t_now):
    #         tvp_template['_tvp',:, 'Theta'] = 0.1
    #         return tvp_template
    #     mpc.set_tvp_fun(tvp_fun)
    # else:
    #     tvp_template = mpc.get_tvp_template()
    #     def tvp_fun(t_now):
    #         tvp_template['_tvp',:, 'Theta'] = 0
    #         return tvp_template
    #     mpc.set_tvp_fun(tvp_fun)

    socs.append(float(x0[0]))
    vss.append(float(x0[1]))
    vls.append(float(x0[2]))
    tbs.append(float(x0[3]))
    us.append(float(u0[0]))

    time_list.append(toc-tic)

time_arr = np.array(time_list)
mean = np.round(np.mean(time_arr[1:])*1000)
var = np.round(np.std(time_arr[1:])*1000)
print('mean runtime:{}ms +- {}ms for MPC step'.format(mean, var))

time_elapsed = (time.process_time() - time_start)
print(time_elapsed)

tbsDeg = []
for i in tbs:
    tbsDeg.append(i-273.15)


plt.subplot(4,1,1)
# plt.figure(figsize=(10,4))
plt.plot(socs)
#plt.fill_between(DistCum,theta_vals_int,alpha=0.1)
plt.ylabel("SoC (-)")
plt.xlabel("Temps (s)")
plt.grid()
# plt.show()

plt.subplot(4,1,2)
# plt.figure(figsize=(10,4))
plt.plot(vss)
#plt.fill_between(DistCum,theta_vals_int,alpha=0.1)
plt.ylabel("Tension S (V)")
plt.xlabel("Temps (s)")
plt.grid()
# plt.show()

plt.subplot(4,1,3)
# plt.figure(figsize=(10,4))
plt.plot(vls)
#plt.fill_between(DistCum,theta_vals_int,alpha=0.1)
plt.ylabel("Tension L (V)")
plt.xlabel("Temps (s)")
plt.grid()
# plt.show()

plt.subplot(4,1,4)
# plt.figure(figsize=(10,4))
plt.plot(tbsDeg)
#plt.fill_between(DistCum,theta_vals_int,alpha=0.1)
plt.ylabel("Temperature (degree)")
plt.xlabel("Temps (s)")
plt.grid()
# plt.show()

# plt.subplot(4,1,4)
plt.figure(figsize=(10,4))
plt.plot(us)
#plt.fill_between(DistCum,theta_vals_int,alpha=0.1)
plt.ylabel("Courrant (A)")
plt.xlabel("Temps (s)")
plt.grid()
plt.show()