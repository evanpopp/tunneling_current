import numpy as np
import matplotlib.pyplot as plt

from my_libs import Probs as prob
from my_libs import Constants as const
from my_libs import T_E as tecore

q = const.q

def TMM_QW_Prob_Plot(potential, param):
    Prob = []
    for val in potential:
        Prob.append(prob.TMM_QW_Prob(val, (val + 0.00001)*q, param))
    
    plt.plot(potential, np.log10(Prob))
    plt.title("Applied Bias vs Tunneling Probability")
    plt.xlabel("Bias (V)")
    plt.ylabel("Probability (log10)")
    plt.show()
    
def Basic_QW_Prob_Plot(potential, param):
    Prob = []
    for val in potential:
        Prob.append(prob.Basic_QW_Prob(val, (val + 0.00001)*q, param))
    
    plt.plot(potential, np.log10(Prob))
    plt.title("Applied Bias vs Tunneling Probability")
    plt.xlabel("Bias (V)")
    plt.ylabel("Probability (log10)")
    plt.show()
    
def QW_Prob_Plot(potential, param):
    TMM = []
    Basic = []
    for val in potential:
        Basic.append(prob.Basic_QW_Prob(val, (val + 0.00001)*q, param))
        TMM.append(prob.TMM_QW_Prob(val, (val + 0.00001)*q, param))
    
    plt.plot(potential, np.log10(TMM), label = "TMM")
    plt.plot(potential, np.log10(Basic), label = "Basic")
    plt.title("Applied Bias vs Tunneling Probability")
    plt.xlabel("Bias (V)")
    plt.ylabel("Probability (log10)")
    plt.legend()
    plt.show()
    
def Barrier_Prob_Plot(potential, param):
    Prob = []
    for val in potential:
        Prob.append(prob.Barrier_Prob(val, (val + 0.00001)*q, param))
        
    plt.plot(potential, np.log10(Prob))
    plt.title("Applied Bias vs Tunneling Probability")
    plt.xlabel("Bias (V)")
    plt.ylabel("Probability (log10)")
    plt.show()
    
def Barrier_I_Plot(potential, param, T, A):
    I = tecore.Barrier_Current(potential, param, T, A)
        
    plt.semilogy(potential, I)
    plt.title("Applied Bias vs Tunneling Current")
    plt.xlabel("Bias (V)")
    plt.ylabel("Current")
    plt.show()
    
def Basic_QW_I_Plot(potential, param, T, A):
    I = tecore.Basic_QW_Current(potential, param, T, A)
    
    plt.semilogy(potential, I)
    plt.title("Applied Bias vs Tunneling Current")
    plt.xlabel("Bias (V)")
    plt.ylabel("Current")
    plt.show()
    
def TMM_QW_I_Plot(potential, param, T, A):
    I = tecore.TMM_QW_Current(potential, param, T, A)
    
    plt.semilogy(potential, I)
    plt.title("Applied Bias vs Tunneling Current")
    plt.xlabel("Bias (V)")
    plt.ylabel("Current")
    plt.show()
    
def QW_I_Plot(potential, param, T, A):
    Basic = tecore.Basic_QW_Current(potential, param, T, A)
    TMM = tecore.TMM_QW_Current(potential, param, T, A)
    
    plt.semilogy(potential, Basic, label = "Basic")
    plt.semilogy(potential, TMM, label = "TMM")
    plt.title("Applied Bias vs Tunneling Current")
    plt.xlabel("Bias (V)")
    plt.ylabel("Current")
    plt.legend()
    plt.show()
    

    
    
def plot_quantum_well(params):
    """
    Plots the conduction band energy diagram of a quantum well heterostructure.
    Assumes Vacuum Level is at 0 Energy.
    """
    # 1. Define region widths 
    # (Assigning an arbitrary width to the contacts for visualization purposes)
    w_contact = params['lBar_Thickness'] * 1.5 
    w_lBar = params['lBar_Thickness']
    w_qw = params['QW_Length']
    w_rBar = params['rBar_Thickness']
    
    # 2. Calculate x-axis boundaries for each region
    x0 = 0
    x1 = x0 + w_contact
    x2 = x1 + w_lBar
    x3 = x2 + w_qw
    x4 = x3 + w_rBar
    x5 = x4 + w_contact # Assuming symmetric right contact for visual completion
    
    # 3. Calculate Energy levels (Conduction Band = -Electron Affinity)
    e_cond = -params['EA_Cond'] / q
    e_lBar = -params['EA_lBarrier'] / q
    e_qw = -params['QW_EA'] / q
    e_rBar = -params['EA_rBarrier'] / q
    
    # Generate x and y arrays for the step-like potential profile
    # (Assuming the right contact has the same EA as the left contact. Modify if asymmetric)
    x = [x0, x1, x1, x2, x2, x3, x3, x4, x4, x5]
    y = [e_cond, e_cond, e_lBar, e_lBar, e_qw, e_qw, e_rBar, e_rBar, e_cond, e_cond] 

    # 4. Initialize Plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot the Conduction Band Edge
    ax.plot(x, y, color='black', linewidth=2, label='Conduction Band ($E_c$)')
    
    # 5. Shade the physical regions for clarity
    ax.axvspan(x0, x1, color='gray', alpha=0.2, label='Left Contact')
    ax.axvspan(x1, x2, color='blue', alpha=0.15, label='Left Barrier')
    ax.axvspan(x2, x3, color='red', alpha=0.15, label='Quantum Well')
    ax.axvspan(x3, x4, color='blue', alpha=0.15, label='Right Barrier')
    ax.axvspan(x4, x5, color='gray', alpha=0.2, label='Right Contact')
    
    # 6. Plot the Quantum Well Energy State
    # (Assuming params['QW_Energy'] is given relative to the bottom of the well)
    e_state = e_qw + params['QW_Energy'] / q
    ax.hlines(e_state, x2, x3, color='green', linestyle='--', linewidth=2, label='QW Bound State')
    
    # 7. Formatting and Labels
    # Note: If your constants are in Joules, the y-axis will be in Joules. 
    # Divide y-values by elementary charge (1.6e-19) if you prefer to plot in eV.
    ax.set_title('Quantum Well Energy Band Diagram', fontsize=14, fontweight='bold')
    ax.set_xlabel('Position', fontsize=12)
    ax.set_ylabel('Energy (eV)', fontsize=12)
    
    # Adjust y-axis limits slightly for better padding
    ax.set_ylim(min(y) - abs(min(y)*0.1), max(y) + abs(max(y)*0.1))
    
    # Clean up legend and grid
    ax.legend(loc='upper right', bbox_to_anchor=(1.25, 1))
    ax.grid(True, linestyle=':', alpha=0.6)
    
    plt.tight_layout()
    plt.show()
    
def plot_barrier(params):
    """
    Plots the conduction band energy diagram of a barrier heterostructure.
    Assumes Vacuum Level is at 0 Energy.
    """
    # 1. Define region widths 
    # (Assigning an arbitrary width to the contacts for visualization purposes)
    w_contact = params['Bar_Thickness'] * 1.5 
    w_Bar = params['Bar_Thickness']
    
    # 2. Calculate x-axis boundaries for each region
    x0 = 0
    x1 = x0 + w_contact
    x2 = x1 + w_Bar
    x3 = x2 + w_contact # Assuming symmetric right contact for visual completion
    
    # 3. Calculate Energy levels (Conduction Band = -Electron Affinity)
    e_cond = -params['EA_Cond'] / q
    e_Bar = -params['EA_Barrier'] / q
    
    # Generate x and y arrays for the step-like potential profile
    # (Assuming the right contact has the same EA as the left contact. Modify if asymmetric)
    x = [x0, x1, x1, x2, x2, x3]
    y = [e_cond, e_cond, e_Bar, e_Bar, e_cond, e_cond] 

    # 4. Initialize Plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot the Conduction Band Edge
    ax.plot(x, y, color='black', linewidth=2, label='Conduction Band ($E_c$)')
    
    # 5. Shade the physical regions for clarity
    ax.axvspan(x0, x1, color='gray', alpha=0.2, label='Left Contact')
    ax.axvspan(x1, x2, color='blue', alpha=0.15, label='Barrier')
    ax.axvspan(x2, x3, color='red', alpha=0.15, label='Right Contact')
    
    # 6. Formatting and Labels
    # Note: If your constants are in Joules, the y-axis will be in Joules. 
    # Divide y-values by elementary charge (1.6e-19) if you prefer to plot in eV.
    ax.set_title('Quantum Barrier Energy Band Diagram', fontsize=14, fontweight='bold')
    ax.set_xlabel('Position', fontsize=12)
    ax.set_ylabel('Energy (eV)', fontsize=12)
    
    # Adjust y-axis limits slightly for better padding
    ax.set_ylim(min(y) - abs(min(y)*0.1), max(y) + abs(max(y)*0.1))
    
    # Clean up legend and grid
    ax.legend(loc='upper right', bbox_to_anchor=(1.25, 1))
    ax.grid(True, linestyle=':', alpha=0.6)
    
    plt.tight_layout()
    plt.show()