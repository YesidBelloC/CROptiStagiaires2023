
import numpy as np
from casadi import *
from casadi.tools import *
import pdb
import sys
sys.path.append('../../')
import do_mpc


def modelCROpti():
    model_type = 'continuous' # either 'discrete' or 'continuous'
    model = do_mpc.model.Model(model_type)

    # Certain parameters
    #Parametros ctes
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
    Rs = 0.0164
    Cs = 710.5678
    Rl = 0.0395
    Cl = 5183.6054

    m  = 0.05 #50gr
    Cb = 1000 #Lithium Specific Heat
    h  = 0.026#10 #Air Convection Param 10.45-c+10*sqrt(c) c=Air speed
    A  = 0.0652 # Hauteur de batterie
    B  = 0.009  # Rayon batterie
    S  = 2*np.pi*A*B

    Text  = 25+273.15 #[°K]

    # States struct (optimization variables):
    SoC = model.set_variable('_x',  'SoC')  # distance
    Vs = model.set_variable('_x',  'Vs')  # speed
    Vl = model.set_variable('_x',  'Vl')  # speed
    Tb = model.set_variable('_x',  'Tb')  # speed

    # Input struct (optimization variables):
    I = model.set_variable('_u',  'I')

    # Fixed parameters:
    # Eff = model.set_variable('_p', 'Eff')

    # algebraic equations
    # Faer = 0.5*pair*SCx*V**2; #Pendiente encontrar relacion entre direction del viento y del vehiculo.

    # Differential equations
    model.set_rhs('SoC', I/(Qcell))
    model.set_rhs('Vs', I/Cs-Vs/(Rs*Cs))
    model.set_rhs('Vl', I/Cl-Vl/(Rl*Cl))
    model.set_rhs('Tb', (R0*I**2 + (Vs**2)/Rs + (Vl**2)/Rl - h*S*(Tb-Text) ) / (m*Cb))

    # Build the model
    model.setup()

    return model