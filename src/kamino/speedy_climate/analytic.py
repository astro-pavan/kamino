import numpy as np
from scipy.optimize import root, minimize, least_squares, newton

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

def albedo(pCO2: float, ag: float, cos_zeta: float=0.6666) -> float:

   tau_ray = 0.19513 * pCO2 # pCO2 in bar
   aa = ((0.5 - 0.75*cos_zeta)*(1 - np.exp(-(tau_ray/cos_zeta))) + 0.75*tau_ray)/(1 + 0.75 * tau_ray)
   aa2 = (0.75*tau_ray)/(1 + 0.75*tau_ray)

   return 1 - ((1 - ag) * (1 - aa))/((1 - ag) * aa2 + (1 - aa2))

def average_instellation(albedo: float, instellation: float) -> float:
   return 0.25 * (1 - albedo) * instellation

def get_instellation(T: float, pCO2: float) -> float:
   # S = (4 * OLR(T, pCO2)) / (1 - albedo(pCO2, alb_ref))
   S = (4 * OLR(T, pCO2)) / (1 - alb_ref)
   return S / S0

def OLR(T: float, pCO2: float) -> float:

   pCO2 = 1e-12 if pCO2 <= 0 else pCO2
   v = np.where(pCO2 < 1, 0.2 * np.log10(pCO2), np.log10(pCO2))
   xi = 0.01 * (T - 250)

   T_mat = np.array([xi**0, xi, xi**2, xi**3, xi**4, xi**5, xi**6])
   Y_mat = np.array([v**0, v, v**2, v**3, v**4])

   T_mat = T_mat[:, 0] if T_mat.ndim == 2 else T_mat
   Y_mat = Y_mat[:, 0] if Y_mat.ndim == 2 else Y_mat

   high_pCO2_LR = float(np.einsum('i,ij,j', T_mat, Bh, Y_mat))
   low_pCO2_LR = float(np.einsum('i,ij,j', T_mat, Bl, Y_mat))

   if pCO2 < 1:
      return I0 + low_pCO2_LR
   else:
      return I0 + high_pCO2_LR

def get_surface_temp(S: float, pCO2: float) -> float:
   
   pCO2 = pCO2 / 1e5
   
   # func = lambda T: average_instellation(albedo(pCO2, alb_ref), S) - OLR(T, pCO2)
   func = lambda T: np.abs(average_instellation(albedo(pCO2, alb_ref), S) - OLR(T, pCO2))

   sol = least_squares(func, 273)

   T_res = sol.x[0]
   T_res = smooth_max(T_res, 272)

   return T_res
