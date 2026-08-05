import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

from my_libs import Constants as const

a = const.a ## Fine structure constant
c = const.c ## Speed of light in a vacuum
hb = const.hb ## reduced plancks constant
h = const.h ## Plancks constant
me = const.me ##electron mass (kg)
q = const.q ## electron charge
kb = const.kb ## Boltzmans constant

def e_vel(kz, m_cond):   ## electron velocity
    ''' 
        kz = wavevector of electron
        m_eff = effective mass of electron in conducting medium
    '''
    return (hb*kz)/(me*m_cond)

def k_vector(E, U, m_cond): ## electron k-vector
    '''
        E = energy of electron
        U = barrier height
        m_eff = effective mass of electron in conducting medium
    '''
    return np.sqrt(2*me*m_cond*(E-U))/hb

def decay(U, E, m_bar):  ## decay constant inside of a barrier
    '''
        E = energy of electron
        U = barrier height
        m_bar = effective mass of electron in the barrier
    '''
    return np.sqrt(2*me*m_bar*(U-E))/hb

def Tunnel_1Bar(E, U, width, m_bar):   ## Tunneling thru 1x barrier
    '''
        E = energy of electron
        U = barrier height
        width = barrier width
        m_bar = effective mass of electron in the barrier
    '''
    k = decay(U, E, m_bar)
    prefact = 16*E*(U-E)/(U**2)
    result = prefact*np.exp(-2*k*abs(width))
    return result

def BHeight(height, potential): ## height of a barrier with a bias
    '''
        V = Bias applied
        height = original height of the barrier
        potential = voltage across the barrier
    '''
    delta = potential*q
    average = (height + height + delta)/2
    return average

def Elastic_Prob(E_l, U_l, U_r, width_l, width_r, m_lbar, m_rbar):
    '''
        E_l = left side energy
        U_l = left side barrier height
        U_r = right side barrier height
        width_l = left barrier width
        width_r = right barrier width
        m_bar = effective electron mass in the barrier
    '''
    Tunnel_l = Tunnel_1Bar(E_l, U_l, width_l, m_lbar)
    Tunnel_r = Tunnel_1Bar(E_l, U_r, width_r, m_rbar)
    T = Tunnel_l * Tunnel_r
    
    return T