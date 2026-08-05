'''
   !!!!!  MADE BY GEMINI, MODDED BE ME (A HUMAN) !!!!!
'''

import numpy as np
import cmath

from my_libs import Constants as const

# --- Physical Constants (Using meV, nm, and electron mass units) ---
# It's convenient to work in units tailored to quantum semiconductor physics.
# Potential is in meV, distances in nm, and mass in units of the free electron mass (me).
# Using these units, the constant (hbar^2) / (2me) becomes a convenient number.
# In SI: (hbar^2) / (2me) = (1.054e-34 J.s)^2 / (2 * 9.11e-31 kg) ~ 0.609e-38 J.m^2
# Converting to eV.nm^2: 0.609e-38 J.m^2 * (1 eV / 1.6e-19 J) * (1e9 nm / 1 m)^2 ~ 3.8e-1 eV.nm^2
# Let's use a standard value: (hbar^2 / (2me)) = 3.81 eV*nm^2 or 3810 meV*nm^2. Wait, that conversion is hard.

# Let's use standard SI for intermediate calculations and then convert. This is more robust.
hbar = const.hb # J s (Planck's reduced constant)
me = const.me # kg (Electron rest mass)
eV_to_J = const.q # Conversion from eV to J

J_to_meV = 1.0 / (eV_to_J * 1e-3) # Conversion from J to meV
nm_to_m = 1e-9 # Conversion from nm to m

# The user requested asymmetric structure and DIFFERENT MATERIALS for the well
# and each barrier. This overrides the simplification in the image where effective mass
# is treated as constant. The effective mass is defined as a multiplier of the free electron mass (me).

def solve_transmission_asymmetric(potentials_meV, boundary_positions_nm, effective_mass_multipliers, E, FWHM):
    """
    Solves for the transmission probability spectrum using the Transfer Matrix Method (TMM).
    Supports a generalized asymmetric structure with different effective masses in each region.

    Parameters:
    - potentials_meV: A list or array of V levels for each region (e.g., [V1, V2, V3, V4, V5]).
    - boundary_positions_nm: A list or array of the x-coordinates of each interface (e.g., [x1, x2, x3, x4]).
    - effective_mass_multipliers: A list or array of effective mass multipliers for each region (e.g., [m1*, m2*, m3*, m4*, m5*]).
    - E_array_meV: A list or array of incident particle energy values to solve for.

    Returns:
    - transmission_probabilities: A numpy array of T values for each energy in E_array.
    """
    E_J = E + 1j*FWHM
    # --- Preliminary Setup and Pre-calculation ---
    # The image describes 5 regions and 4 interfaces (R1 to R4).
    num_interfaces = len(boundary_positions_nm)
    # The image equation [A1, B1] = R1 R2 R3 R4 [A5, B5] relates the FIRST to the LAST.
    # To relate the FIRST (region 1) to the SECOND (region 2), the total matrix would be R1.
    # The code below implements this logic by iterating interface-by-interface.

    M_total = np.eye(2, dtype=complex) # Initialize the total matrix as an identity matrix

        # --- Transfer Matrix Multiplication ---
        # Matrix to transform from Region i to Region i+1 is called M_{i -> i+1}
        # [A1, B1] = M_{1 -> 2} * [A2, B2]
        # [A2, B2] = M_{2 -> 3} * [A3, B3]
        # ...
        # [A1, B1] = (M_{1 -> 2} * M_{2 -> 3} * ... * M_{N-1 -> N}) * [AN, BN]
        # [A1, B1] = M_total * [AN, BN]

        # The order of multiplication from the equation in the picture ([A1, B1] = R1 R2 R3 R4 [A5, 0])
        # means R1 converts 1->2, R2 converts 2->3, etc. and they multiply from LEFT TO RIGHT.
        # This is equivalent to applying interface matrices in sequence.

        # Calculate initial conditions (k and mass in Region 1)
    k_prev = cmath.sqrt(2.0 * me * effective_mass_multipliers[0] * (E_J - (potentials_meV[0] / J_to_meV))) / hbar
    m_eff_prev = me * effective_mass_multipliers[0]

        # First boundary (x1 between Region 1 and Region 2)
        # Apply matrix R1 in the notation of the user's image.
    x_interface_nm = boundary_positions_nm[0]
    v_next_meV = potentials_meV[1]
    m_eff_multiplier_next = effective_mass_multipliers[1]

    k_next = cmath.sqrt(2.0 * me * m_eff_multiplier_next * (E_J - (v_next_meV / J_to_meV))) / hbar
    m_eff_next = me * m_eff_multiplier_next

    R_interface = get_interface_matrix_general(k_prev, m_eff_prev, k_next, m_eff_next, x_interface_nm)
    M_total = M_total @ R_interface # Multiply to the right

        # Update and continue for the remaining boundaries
    k_prev, m_eff_prev = k_next, m_eff_next

    for i in range(1, num_interfaces):
            x_interface_nm = boundary_positions_nm[i]
            # Next region index is i+1 (e.g., at x2, it's region 3)
            v_next_meV = potentials_meV[i+1]
            m_eff_multiplier_next = effective_mass_multipliers[i+1]

            # Recalculate k for the next region
            k_next = cmath.sqrt(2.0 * me * m_eff_multiplier_next * (E_J - (v_next_meV / J_to_meV))) / hbar
            m_eff_next = me * m_eff_multiplier_next

            R_interface = get_interface_matrix_general(k_prev, m_eff_prev, k_next, m_eff_next, x_interface_nm)
            M_total = M_total @ R_interface # Multiply from left to right, as indicated in the user's diagram equation: [A1, B1] = R1 R2...

            # Update for the next iteration
            k_prev, m_eff_prev = k_next, m_eff_next


        # --- Transmission Coefficient and Flux Factor ---
        # The equation from the image is [A1, B1] = M_total * [A5, 0] = [M_total[0,0]*A5, M_total[1,0]*A5]
        # This gives A1 = M_total[0,0] * A5. The transmission amplitude t = A5 / A1 = 1 / M_total[0,0].
        # In a generalized system (asymmetric, different masses), the transmission probability T is:
        # T = |t|^2 * (Transmitted Flux / Incident Flux)
        # Incident flux is proportional to |A1|^2 * k_in / m_in_eff
        # Transmitted flux is proportional to |A5|^2 * k_out / m_out_eff
        # So, T = |1 / M_total[0,0]|^2 * (k_out / m_out_eff) / (k_in / m_in_eff)

    M11 = M_total[0, 0]
    if abs(M11) == 0:
        T = 1 # Special case, practically shouldn't happen
    else:
        # In and out parameters
        k_in = cmath.sqrt(2.0 * me * effective_mass_multipliers[0] * (E_J - (potentials_meV[0] / J_to_meV))) / hbar
        k_out = cmath.sqrt(2.0 * me * effective_mass_multipliers[-1] * (E_J - (potentials_meV[-1] / J_to_meV))) / hbar
        m_in_eff = me * effective_mass_multipliers[0]
        m_out_eff = me * effective_mass_multipliers[-1]

        t_amp = 1.0 / M11
        T = abs(t_amp)**2 * (k_out.real / m_out_eff) / (k_in.real / m_in_eff) # Only take real flux part to handle complex k (tunneling regions)
        
    return T

def get_interface_matrix_general(k1, m1_eff, k2, m2_eff, x_interface_nm):
    """
    Constructs the 2x2 generalized interface matrix (transfer matrix from region 1 to region 2).
    Matches ψ and (1/m)dψ/dx continuity at the interface.
    This formula is derived to support different masses, in contrast to the constant-mass simplification shown in the user's image.
    When m1_eff = m2_eff, this form reduces to the exact one in the user's image.

    The correct form for differing masses is:
    R = [1/2E1] * [ E2*(1+K)    E2^{-1}*(1-K) ]
        [E1/2]  * [ E2*(1-K)    E2^{-1}*(1+K) ]
    where E1 = e^{ik1x}, E2 = e^{ik2x}, and K = (m1_eff * k2) / (k1 * m2_eff)
    """
    
    # Check for divide-by-zero if k1 is very close to zero
    '''
    if abs(k1) < 1e-18:
        k1 = 1e-18
    
    if abs(k2) < 1e-18:
        k2 = 1e-18
    '''
    x_interface_m = x_interface_nm * nm_to_m
    ##E1 = cmath.exp(1j * k1 * x_interface_m)
    ##E2 = cmath.exp(1j * k2 * x_interface_m)
    
    ##phase_diff = cmath.exp(1j * (k2 - k1) * x_interface_m)
    ##phase_add = cmath.exp(1j * (k2 + k1) * x_interface_m)
    
    diff_phase = 1j * (k2 - k1) * x_interface_m
    add_phase = 1j * (k2 + k1) * x_interface_m

    K = (m1_eff * k2) / (k1 * m2_eff)

    # Building the generalized matrix element-by-element (which is what the previous function's D^-1 * D form would do)
    # The image formula assumes constant mass, leading to K = k2/k1.
    # The general form is:
    # (1,1) = (E2/E1) * (1+K)/2 = e^{i(-k1+k2)x} * (1+K)/2
    # This matches the user's image term: (k1+k2)e^{i(-k1+k2)x} * 1/2k1 = e^{i(-k1+k2)x} * (1+k2/k1)/2

    R11 = cmath.exp(diff_phase) * (1.0 + K) / 2.0 
    R12 = cmath.exp(-add_phase) * (1.0 - K) / 2.0  # (E2^-1 / E1) * (1-K)/2 = e^{i(-k1-k2)x} * (1-K)/2. Correct.
    R21 = cmath.exp(add_phase) * (1.0 - K) / 2.0  # (E1 * E2) * (1-K)/2 = e^{i(k1+k2)x} * (1-K)/2. Correct.
    R22 = cmath.exp(-diff_phase) * (1.0 + K) / 2.0 # (E1 * E2^-1) * (1+K)/2 = e^{i(k1-k2)x} * (1+K)/2. Correct.

    return np.array([[R11, R12], [R21, R22]], dtype=complex)