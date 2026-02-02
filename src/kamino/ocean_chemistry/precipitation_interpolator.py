import numpy as np
import pickle
import os
from scipy.interpolate import RegularGridInterpolator
from kamino.ocean_chemistry.fast_chemistry import get_calcite_data_fast
from kamino.ocean_chemistry.precipitation import get_calcite_precipitation_rate
from tqdm import tqdm

class PrecipitationInterpolator:
    def __init__(self, filename='precip_table_5D.pkl'):
        self.filename = filename
        self.interp_rate = None
        self.interp_SI = None
        
        # --- DEFINITION OF SAFE BOUNDS ---
        self.T_MIN, self.T_MAX = 270.0, 380.0
        self.P_MIN, self.P_MAX = 1e5, 1000e5       # 1 bar to 1000 bar (Deep Ocean support)
        self.DIC_MIN, self.DIC_MAX = 1e-6, 1.0     
        self.CA_MIN, self.CA_MAX   = 1e-6, 1.0     
        self.ALK_RATIO_MIN, self.ALK_RATIO_MAX = 0.1, 2.0
        
        # Grid Definition
        # Pressure Grid: Log-spaced to capture surface (1 bar) and depth (500 bar) accurately
        self.T_grid   = np.linspace(self.T_MIN, self.T_MAX, 20)
        self.P_grid   = np.geomspace(self.P_MIN, self.P_MAX, 5) 
        self.DIC_grid = np.geomspace(self.DIC_MIN, self.DIC_MAX, 15)
        self.Ca_grid  = np.geomspace(self.CA_MIN, self.CA_MAX, 15)
        self.Alk_ratio_grid = np.geomspace(self.ALK_RATIO_MIN, self.ALK_RATIO_MAX, 15)

    def generate_table(self):
        print("Generating 5D Precipitation Table (T, P, DIC, Ca, Ratio)...")
        
        # Pre-allocate arrays (5 Dimensions)
        shape = (len(self.T_grid), len(self.P_grid), len(self.DIC_grid), len(self.Ca_grid), len(self.Alk_ratio_grid))
        rate_data = np.zeros(shape)
        SI_data   = np.zeros(shape)
        
        total_points = rate_data.size
        print(f"  Total Grid Points: {total_points:,}")
        
        # Loop through parameter space
        for i, T in enumerate(self.T_grid):
            print(f"  Processing T={T:.1f}K... ({i+1}/{len(self.T_grid)})")
            
            for j, P in tqdm(enumerate(self.P_grid)):
                # Inner loops
                for k, DIC in enumerate(self.DIC_grid):
                    for l, Ca in enumerate(self.Ca_grid):
                        for m, ratio in enumerate(self.Alk_ratio_grid):
                            
                            Alk = DIC * ratio
                            
                            try:
                                # 1. Run the Physics Solver (With Pressure P)
                                # P passed in Pascals
                                # rate, SI = get_calcite_data_fast(P, T, Alk, DIC, Ca, rate_constant=1.0)
                                rate, SI = get_calcite_precipitation_rate(P, T, Alk, DIC, Ca)
                                
                                if not np.isfinite(rate) or not np.isfinite(SI):
                                    raise ValueError("NaN")
                                    
                                rate_data[i,j,k,l,m] = max(0.0, rate)
                                SI_data[i,j,k,l,m]   = SI
                                
                            except Exception:
                                rate_data[i,j,k,l,m] = 0.0
                                SI_data[i,j,k,l,m]   = -99.0

        print("  Table Generation Complete.")

        # Create 5D Interpolators
        # Dimensions: T, log(P), log(DIC), log(Ca), log(Ratio)
        self.interp_rate = RegularGridInterpolator(
            (self.T_grid, np.log10(self.P_grid), np.log10(self.DIC_grid), np.log10(self.Ca_grid), np.log10(self.Alk_ratio_grid)), 
            rate_data, bounds_error=False, fill_value=None
        )
        
        self.interp_SI = RegularGridInterpolator(
            (self.T_grid, np.log10(self.P_grid), np.log10(self.DIC_grid), np.log10(self.Ca_grid), np.log10(self.Alk_ratio_grid)), 
            SI_data, bounds_error=False, fill_value=None
        )
        
        self.save()

    def save(self):
        with open(self.filename, 'wb') as f:
            pickle.dump((self.interp_rate, self.interp_SI), f)

    def load(self):
        if os.path.exists(self.filename):
            try:
                with open(self.filename, 'rb') as f:
                    self.interp_rate, self.interp_SI = pickle.load(f)
                return True
            except:
                return False
        return False

    def get_rate(self, P, T, Alk, DIC, Ca):
        """
        Robust query function with Pressure.
        """
        # 1. HARD CLAMPING
        T_safe   = np.clip(T,   self.T_MIN,   self.T_MAX)
        P_safe   = np.clip(P,   self.P_MIN,   self.P_MAX)
        DIC_safe = np.clip(DIC, self.DIC_MIN, self.DIC_MAX)
        Ca_safe  = np.clip(Ca,  self.CA_MIN,  self.CA_MAX)
        
        if DIC_safe < 1e-9:
            ratio_safe = self.ALK_RATIO_MIN
        else:
            ratio_safe = np.clip(Alk / DIC_safe, self.ALK_RATIO_MIN, self.ALK_RATIO_MAX)

        # 2. Query in Log-Space
        point = [T_safe, np.log10(P_safe), np.log10(DIC_safe), np.log10(Ca_safe), np.log10(ratio_safe)]
        
        # 3. Fetch Data (CRITICAL FIX HERE)
        # RegularGridInterpolator returns an array. We use .item() to extract the scalar.
        rate_array = self.interp_rate(point)
        SI_array   = self.interp_SI(point)
        
        rate = rate_array.item()
        SI   = SI_array.item()
        
        if np.isnan(rate): rate = 0.0
        if np.isnan(SI):   SI = -99.0
            
        return max(0.0, rate), SI

# --- SINGLETON ACCESSOR ---
_interpolator = None

def get_calcite_data_interpolated(P, T, Alk, DIC, Ca):
    global _interpolator
    if _interpolator is None:
        _interpolator = PrecipitationInterpolator()
        if not _interpolator.load():
            _interpolator.generate_table()
    
    return _interpolator.get_rate(P, T, Alk, DIC, Ca)