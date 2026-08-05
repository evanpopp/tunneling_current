import numpy as np
import matplotlib.pyplot as plt

from my_libs import Quantum_Core as qcore
from my_libs import TMM as tmm
from my_libs import Constants as const

q = const.q
hb = const.hb 

def Barrier_Prob(V, E, param):
    EA_Cond = param["EA_Cond"]
    EA_Bar = param["EA_Barrier"]
    m_bar = param["m*_Barrier"]
    width = param["Bar_Thickness"]
    
    l_bar_pot = -V
    bar_height = EA_Cond - EA_Bar
    U = qcore.BHeight(bar_height, l_bar_pot) + V*q
    T = qcore.Tunnel_1Bar(E, U, width, m_bar)
    return T

def Basic_QW_Prob(V, E, param):
    width_l = param["lBar_Thickness"]
    EA_lbar = param["EA_lBarrier"]
    m_lbar = param["m*_lBarrier"]
    width_r = param["rBar_Thickness"]
    EA_rbar = param["EA_rBarrier"]
    m_rbar = param["m*_rBarrier"]
    
    qw_l = param["QW_Length"]
    EA_QW = param["QW_EA"]
    EA_Cond = param["EA_Cond"]
    
    t = param["Dephasing_time"]
    
    x1 = 10e-9
    x2 = x1 + width_l
    x3 = x2 + qw_l
    x4 = x3 + width_r
    
    distance = x4 - x1
    QW_base = EA_Cond - EA_QW
    qw_pos = (width_r + qw_l/2)/distance
    Uqw = QW_base + qw_pos*V*q
    lbar_height = EA_Cond - EA_lbar
    rbar_height = EA_Cond - EA_rbar
    l_bar_pot = -V*width_l/distance
    r_bar_pot = V*width_r/distance

    U_l = qcore.BHeight(lbar_height, l_bar_pot) + V*q
    U_r = qcore.BHeight(rbar_height, r_bar_pot)
    
    Tunnel_l = qcore.Tunnel_1Bar(E, U_l, width_l, m_lbar)
    Tunnel_r = qcore.Tunnel_1Bar(E, U_r, width_r, m_rbar)
    T = Tunnel_l * Tunnel_r
    return T

def TMM_QW_Prob(V, E, param):
    width_l = param["lBar_Thickness"]
    EA_lbar = param["EA_lBarrier"]
    m_lbar = param["m*_lBarrier"]
    width_r = param["rBar_Thickness"]
    EA_rbar = param["EA_rBarrier"]
    m_rbar = param["m*_rBarrier"]
    
    qw_l = param["QW_Length"]
    EA_QW = param["QW_EA"]
    m_qw = param["m*_QW"]
    
    EA_Cond = param["EA_Cond"]
    m_cond = param["m*_Cond"]
    
    t = param["Dephasing_time"]
    
    FWHM = hb/t
    
    x1 = 10e-9
    x2 = x1 + width_l
    x3 = x2 + qw_l
    x4 = x3 + width_r
    
    distance = x4 - x1
    QW_base = EA_Cond - EA_QW
    qw_pos = (width_r + qw_l/2)/distance
    Uqw = QW_base + qw_pos*V*q
    lbar_height = EA_Cond - EA_lbar
    rbar_height = EA_Cond - EA_rbar
    l_bar_pot = -V*width_l/distance
    r_bar_pot = V*width_r/distance
    
    Ul = qcore.BHeight(lbar_height, l_bar_pot) + V*q
    Ur = qcore.BHeight(rbar_height, r_bar_pot)
    
    potentials_meV = [0, 1000*Ul/q, 1000*Uqw/q, 1000*Ur/q, 0]
    boundary_positions_nm = [x1*1e9, x2*1e9, x3*1e9, x4*1e9]
    effective_mass_multipliers = [m_cond, m_lbar, m_qw, m_rbar, m_cond]
    Prob = tmm.solve_transmission_asymmetric(potentials_meV, boundary_positions_nm, effective_mass_multipliers, E, FWHM)
    return Prob