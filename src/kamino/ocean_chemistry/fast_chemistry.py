import numpy as np
from kamino.utils import smooth_max

# ==========================================
# 1. DATABASE COEFFICIENTS (From basic_v2.txt)
# ==========================================

# Equilibrium Constants
CALCITE_ANALYTIC = [-1.4978e2, -4.8370e-2, 4.8974e3, 6.0458e1, 7.6464e1]
K1_ANALYTIC = [-1.0534e1, 2.1746e-2, 2.5216e3, 7.9125e-1, 3.9351e1]
K1_SIGN = -1.0 
K2_ANALYTIC = [-6.9958e1, -3.3526e-2, -7.0846e1, 2.8224e1, -1.0849]
KW_ANALYTIC = [-6.7506e1, -3.0619e-2, -1.9901e3, 2.8004e1, -3.1033e1]

# Ion Pair Constants
K_CaHCO3_ANALYTIC = [5.5985e1, 3.4639e-2, -3.6972e2, -2.5864e1, -5.7859e0]
K_CaCO3_FILE = [2.3045e2, 5.5350e-2, -8.5056e3, -9.1096e1, -1.3279e2]

# B-dot / Debye-Huckel Parameters (Table from LLNL_AQUEOUS_MODEL_PARAMETERS)
DH_TEMP_C = np.array([0.0, 25.0, 60.0, 100.0, 150.0])
# A parameter (Numerator)
DH_A_VALS = np.array([0.4939, 0.5114, 0.5465, 0.5995, 0.6855])
# B parameter (Denominator, multiplies ion size)
DH_B_VALS = np.array([0.3253, 0.3288, 0.3346, 0.3421, 0.3525])
# B-dot parameter (Linear term)
BDOT_VALS = np.array([0.0374, 0.0410, 0.0438, 0.0460, 0.0470])

# Ion Sizes (a_i) from SOLUTION_SPECIES / standard LLNL values
# Used in: 1 + B * a * sqrt(I)
SIZE_Ca = 6.0
SIZE_HCO3 = 4.0   # Default for monovalent
SIZE_CO3 = 4.5    # Specific for carbonate
SIZE_H = 9.0
SIZE_OH = 3.5
SIZE_CaHCO3 = 4.0 # Estimate for charged pair
SIZE_DEFAULT = 4.0

# Kinetics
R_GAS = 8.314
AA = 5.625   
AC = 62.5    
EA = 16000   
EAC = 48000  
KC = 160     

# ==========================================
# 2. HELPER FUNCTIONS
# ==========================================

def calc_log_k(T_kelvin, coeffs, sign=1.0):
    A1, A2, A3, A4, A5 = coeffs
    log_k = A1 + A2*T_kelvin + A3/T_kelvin + A4*np.log10(T_kelvin) + A5/(T_kelvin**2)
    return sign * log_k

def get_bdot_params(T_kelvin):
    """Interpolates A, B, and bdot parameters."""
    T_C = T_kelvin - 273.15
    T_C = np.clip(T_C, 0, 150)
    
    A = np.interp(T_C, DH_TEMP_C, DH_A_VALS)
    B = np.interp(T_C, DH_TEMP_C, DH_B_VALS)
    bdot = np.interp(T_C, DH_TEMP_C, BDOT_VALS)
    return A, B, bdot

def get_activity_coefficient_bdot(z, I, size_param, A, B, bdot):
    """
    Calculates gamma using the B-dot (Extended Debye-Huckel) equation.
    log gamma = - A z^2 sqrt(I) / (1 + B a sqrt(I)) + bdot * I
    """
    if I < 1e-9: return 1.0
    sqrt_I = np.sqrt(I)
    
    # B-dot Equation
    # Note: bdot term (linear) is usually added for neutral salting out effect
    log_gamma = - (A * (z**2) * sqrt_I) / (1.0 + B * size_param * sqrt_I) + (bdot * I)
    
    return 10**log_gamma

# ==========================================
# 3. MAIN SOLVER
# ==========================================

def get_calcite_data_fast(P, T, Alk, DIC, Ca, Mg=0, S_area=1e-6, scaling_factor=1.0):
    
    T = smooth_max(273.5, T)
    
    # --- 1. Thermodynamics ---
    log_Ksp_calc = calc_log_k(T, CALCITE_ANALYTIC)
    log_K1 = calc_log_k(T, K1_ANALYTIC, sign=K1_SIGN)
    log_K2 = calc_log_k(T, K2_ANALYTIC)
    log_Kw = calc_log_k(T, KW_ANALYTIC)
    log_K_CaHCO3 = calc_log_k(T, K_CaHCO3_ANALYTIC)
    log_K_CaCO3 = calc_log_k(T, K_CaCO3_FILE) - log_K2

    K1 = 10**log_K1
    K2 = 10**log_K2
    Kw = 10**log_Kw
    K_CaHCO3 = 10**log_K_CaHCO3
    K_CaCO3  = 10**log_K_CaCO3
    Ksp = 10**(log_Ksp_calc + log_K2)

    # --- 2. Activity Parameters (B-dot Model) ---
    A_dh, B_dh, bdot = get_bdot_params(T)

    # --- 3. Solver ---
    I = 0.5 * (4*Ca + 4*Mg + 1*Alk + 1*Alk) 
    I = max(I, 1e-7)

    # Calculate individual activity coefficients
    g_Ca    = get_activity_coefficient_bdot(2, I, SIZE_Ca, A_dh, B_dh, bdot)
    g_HCO3  = get_activity_coefficient_bdot(1, I, SIZE_HCO3, A_dh, B_dh, bdot)
    g_CO3   = get_activity_coefficient_bdot(2, I, SIZE_CO3, A_dh, B_dh, bdot)
    g_H     = get_activity_coefficient_bdot(1, I, SIZE_H, A_dh, B_dh, bdot)
    g_OH    = get_activity_coefficient_bdot(1, I, SIZE_OH, A_dh, B_dh, bdot)
    g_CaHCO3= get_activity_coefficient_bdot(1, I, SIZE_CaHCO3, A_dh, B_dh, bdot)
    
    # Apparent Constants (K' = K * g_reactants / g_products)
    K1_app = K1 / (g_H * g_HCO3) # Dissociation
    K2_app = K2 * g_HCO3 / (g_H * g_CO3)
    Kw_app = Kw / (g_H * g_OH)
    
    # Pair Formation K'
    # Ca + HCO3 -> CaHCO3+
    K_CaHCO3_app = K_CaHCO3 * (g_Ca * g_HCO3) / g_CaHCO3
    # Ca + CO3 -> CaCO3(aq) (Assuming g_neutral = 1)
    K_CaCO3_app  = K_CaCO3  * (g_Ca * g_CO3) / 1.0

    h = 10**(-8.0)
    Ca_free = Ca 
    DIC_free = DIC
    
    for _ in range(30):
        # 1. Carbonate Fractions
        denom = h**2 + h*K1_app + K1_app*K2_app
        
        # Safety check: If h exploded, reset or break to avoid NaN
        if np.isinf(denom):
            h = 1e-7 # Reset to neutral if blown up
            denom = h**2 + h*K1_app + K1_app*K2_app

        alpha1 = h*K1_app / denom      # HCO3
        alpha2 = K1_app*K2_app / denom # CO3
        
        # 2. Iterative Mass Balance (Picard)
        for _inner in range(5):
            HCO3 = DIC_free * alpha1
            CO3  = DIC_free * alpha2
            
            # Ca Mass Balance
            denom_Ca = 1.0 + K_CaHCO3_app*HCO3 + K_CaCO3_app*CO3
            Ca_free = Ca / denom_Ca
            
            # DIC Mass Balance
            bound_C = (K_CaHCO3_app*HCO3 + K_CaCO3_app*CO3) * Ca_free
            DIC_free = DIC - bound_C

        HCO3_final = DIC_free * alpha1
        CO3_final  = DIC_free * alpha2
        OH_final   = Kw_app / h
        
        # 3. Alkalinity Balance
        CaHCO3_conc = K_CaHCO3_app * Ca_free * HCO3_final
        CaCO3_conc  = K_CaCO3_app  * Ca_free * CO3_final
        
        Alk_calc = HCO3_final + 2*CO3_final + OH_final - h + CaHCO3_conc + 2*CaCO3_conc
        
        diff = Alk - Alk_calc
        if abs(diff) < 1e-12: break
        
        # --- ROBUST UPDATE STEP ---
        # Calculate step based on relative error
        # Add a small epsilon to denominator to prevent divide-by-zero if Alk is 0
        step = -0.8 * diff / (abs(Alk) + 1e-12)
        
        # CLAMP THE STEP:
        # Prevent h from changing by more than e^2 (~7.4x) or e^-2 (~0.14x) in one step.
        # This prevents the "Infinity" explosion.
        step = np.clip(step, -2.0, 2.0)
        
        h = h * np.exp(step)
        
        # Absolute safety clamp for pH (prevent going outside pH 0-14 range)
        h = np.clip(h, 1e-14, 1.0)

    # --- 4. Saturation Index ---
    act_Ca = Ca_free * g_Ca
    act_CO3 = CO3_final * g_CO3
    
    Omega = (act_Ca * act_CO3) / Ksp
    SI = np.log10(Omega + 1e-16)

    # --- 5. Kinetics (Linear Law) ---
    act_HCO3 = HCO3_final * g_HCO3
    act_H = h * g_H
    act_c = act_HCO3 + act_CO3
    
    carb_term = 1 - (KC * act_c) / (1 + KC * act_c)
    
    rplusa = AA * np.exp(-EA / (R_GAS * T)) * act_H * S_area
    rplusc = AC * np.exp(-EAC / (R_GAS * T)) * carb_term
    
    rplus = (rplusa + rplusc) * scaling_factor
    
    if Omega <= 1.0:
        return 0.0, SI
        
    rate_val = rplus * (Omega - 1.0) ** 2

    return rate_val, SI