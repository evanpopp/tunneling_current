import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sci

from my_libs import Constants as const
from my_libs import Probs as prob

## Constants ##
a = const.a ## Fine structure constant
c = const.c ## Speed of light in a vacuum
hb = const.hb ## reduced plancks constant
h = const.h ## Plancks constant
me = const.me ##electron mass (kg)
q = const.q ## electron charge
kb = const.kb ## Boltzmans constant

## Integral limits ##
int_low = const.E_LOW
int_high = const.E_HIGH

def Constant(m_eff):
    return 4*np.pi*m_eff*me*q/(h**3)

def Fermi(Ex, Ef, T):
    factor = -(Ex - Ef)/(kb*T)
    return np.logaddexp(0, factor)

def N(E, Ef1, Ef2, T):
    return kb*T*(Fermi(E, Ef1, T) - Fermi(E, Ef2, T))

def SF_El(V, E, T):
    Ef1 = V*q
    Ef2 = 0
    return N(E, Ef1, Ef2, T)

def Basic_QW_Current(Voltage, param, T, A): 
    J = lambda E, V: N(E, (V)*q, 0, T)*prob.Basic_QW_Prob(V, E, param)
    I = []
    for val in Voltage:
        result, error = sci.quad(J, (val - int_low)*q, (val + int_high)*q, args = (val,), epsabs=1e-12, epsrel=1e-12)
        I.append(result*Constant(param["m*_lBarrier"])*A) 
    return I

def TMM_QW_Current(Voltage, param, T, A):
    J = lambda E, V: SF_El(V, E, T)*prob.TMM_QW_Prob(V, E, param)
    I = []
    for val in Voltage:
        result, error = sci.quad(J, (val - int_low)*q, (val + int_high)*q, args = (val,), epsabs=1e-12, epsrel=1e-12)
        I.append(result*Constant(param["m*_lBarrier"])*A)  
    return I

def Barrier_Current(Voltage, param, T, A):
    J = lambda E, V: SF_El(V, E, T)*prob.Barrier_Prob(V, E, param)
    I = []
    for val in Voltage:
        result, error = sci.quad(J, (val - int_low)*q, (val + int_high)*q, args = (val,), epsabs=1e-12, epsrel=1e-12)
        I.append(result*Constant(param["m*_Barrier"])*A)  
    return I