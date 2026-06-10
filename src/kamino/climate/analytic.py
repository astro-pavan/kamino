import numpy as np
from scipy.optimize import root, minimize, least_squares, newton, brentq

from kamino.utils import *

I0 = -3.1

# constants for OLR parameterisation for low pco2 < 1bar
Bl = np.array([[87.8373, -311.289, -504.408, -422.929, -134.611],
              [54.9102, -677.741, -1440.63, -1467.04, -543.371],
              [24.7875, 31.3614, -364.617, -747.352, -395.401],
              [75.8917, 816.426, 1565.03, 1453.73, 476.475],
              [43.0076, 339.957, 996.723, 1361.41, 612.967],
              [-31.4994, -261.362, -395.106, -261.600, -36.6589],
              [-28.8846, -174.942, -378.436, -445.878, -178.948]])

# constants for OLR parameterisation for 1 bar < pco2 < 10 bar
Bh = np.array([[87.8373, -52.1056, 35.2800, -1.64935, -3.42858],
              [54.9102, -49.6404, -93.8576, 130.671, -41.1725],
              [24.7875, 94.7348, -252.996, 171.685, -34.7665],
              [75.8917, -180.679, 385.989, -344.020, 101.455],
              [43.0076, -327.589, 523.212, -351.086, 81.0478],
              [-31.4994, 235.321, -462.453, 346.483, -90.0657],
              [-28.8846, 284.233, -469.600, 311.854, -72.4874]])

S0 = 1376

pco2_ref = 280e-6
T_ref = 288

alb_ref = 0.3

def albedo_funtion(pCO2: float, ag: float, cos_zeta: float=0.6666) -> float:

   tau_ray = 0.19513 * pCO2 # pCO2 in bar
   aa = ((0.5 - 0.75*cos_zeta)*(1 - np.exp(-(tau_ray/cos_zeta))) + 0.75*tau_ray)/(1 + 0.75 * tau_ray)
   aa2 = (0.75*tau_ray)/(1 + 0.75*tau_ray)

   return 1 - ((1 - ag) * (1 - aa))/((1 - ag) * aa2 + (1 - aa2))

def average_instellation(albedo: float, instellation: float, tidally_locked) -> float:
   f = 0.66 if tidally_locked else 0.25
   return f * (1 - albedo) * instellation

def get_instellation(T: float, pCO2: float) -> float:
   # S = (4 * OLR(T, pCO2)) / (1 - albedo(pCO2, alb_ref))
   S = (4 * OLR(T, pCO2)) / (1 - alb_ref)
   return S / S0

def OLR(T: float, pCO2: float) -> float:

   pCO2 = 1e-12 if pCO2 <= 0 else pCO2
   log_pCO2 = np.log10(pCO2)
   xi = 0.01 * (T - 250)

   T_mat = np.array([xi**0, xi, xi**2, xi**3, xi**4, xi**5, xi**6])

   v_low  = 0.2 * log_pCO2
   v_high = log_pCO2
   Y_low  = np.array([v_low**0,  v_low,  v_low**2,  v_low**3,  v_low**4])
   Y_high = np.array([v_high**0, v_high, v_high**2, v_high**3, v_high**4])

   low_pCO2_LR  = float(np.einsum('i,ij,j', T_mat, Bl, Y_low))
   high_pCO2_LR = float(np.einsum('i,ij,j', T_mat, Bh, Y_high))

   # Smooth sigmoid blend between the two polynomial fits, centred at 1 bar.
   # Width 0.1 bar keeps the transition region >> the Jacobian perturbation (~0.02 bar),
   # preventing the finite-difference Jacobian from straddling the formula boundary and
   # producing a spurious ~100-year fast mode in the pCO2 ODE.
   blend = 0.5 * (1.0 + np.tanh((pCO2 - 1.0) / 0.1))
   return I0 + (1.0 - blend) * low_pCO2_LR + blend * high_pCO2_LR

def get_T_surface_v0(S, P_CO2, albedo, tidally_locked=False) -> float:
   
   pCO2 = P_CO2 / 1e5
   
   # func = lambda T: average_instellation(albedo(pCO2, alb_ref), S) - OLR(T, pCO2)
   func = lambda T: average_instellation(albedo_funtion(pCO2, albedo), S, tidally_locked) - OLR(T, pCO2)

   sol = newton(func, 273, disp=False)

   T_res = sol

   return float(T_res)

def get_T_surface_analytic(S, P_CO2, albedo, tidally_locked=False):

    pCO2_bar = P_CO2 / 1e5
    A_bond = albedo_funtion(pCO2_bar, albedo)

    f_geo = 0.66 if tidally_locked else 0.25
    F_in = S * (1 - A_bond) * f_geo

    def residual(T):
        return F_in - OLR(np.clip(T, 180, 400), pCO2_bar)

    # The Bl polynomial (low pCO2) has a local OLR maximum near T~360 K and
    # drops sharply toward T=400 K (Komabayashi-Ingersoll saturation plus
    # polynomial boundary artefact). Using T=400 K as the upper bracket end
    # can make residual(400) > 0 even when a stable equilibrium exists at
    # ~300-360 K, hiding the root from brentq. Using T=400 K can also create
    # spurious secondary roots from the OLR dip.
    #
    # Fix: scan [180, 390] K in coarse steps and take the first + -> - sign
    # change (stable equilibrium). This reliably finds the physical root without
    # being confused by the boundary artefact or the unstable branch above the
    # OLR maximum.
    T_nodes = np.linspace(180, 390, 22)   # ~10 K spacing
    r_nodes = [residual(T) for T in T_nodes]

    for i in range(len(T_nodes) - 1):
        if r_nodes[i] >= 0 and r_nodes[i + 1] < 0:
            return float(brentq(residual, T_nodes[i], T_nodes[i + 1]))

    if r_nodes[0] <= 0:
        return 180.0   # Snowball: OLR > F_in even at 180 K
    return 400.0       # Runaway: OLR < F_in everywhere in [180, 390] K
