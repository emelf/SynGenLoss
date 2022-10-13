import numpy as np 
from numpy import cos, sin, sqrt, arctan

from scipy.optimize import brenth, root
from ..DataClasses import GenDataClass
from ..components.GenSaturationModel_v1 import SaturationModel

import cmath as cm
from ..CapDiag import CapabilityDiagram 
from copy import deepcopy

class TrafoCapabilityDiagram(CapabilityDiagram): 
    def __init__(self, gen_data: GenDataClass, sat_data: SaturationModel, X_T_pu: float): 
        self.md = deepcopy(gen_data)
        self.sat = deepcopy(sat_data)
        self.X_T_pu = X_T_pu
        super().__init__(self.md, self.sat)
        
    def get_Q_lims_plant(self, V_hv: float, P_pu: float): 
        """Return -> (P_g, Q_hv_min, Q_hv_max, V_g_min, V_g_max)"""
        def objective(X, min_max, P_g_in): 
            V_g, Q_hv_lim = X 
            P_g, Q_g_min, Q_g_max = self.get_Q_lims(V_g, P_g_in)
            if min_max == 0: 
                Q_g_lim = Q_g_min 
            else: 
                Q_g_lim = Q_g_max 
            delta_g = self._calc_delta_g(V_hv, V_g, P_g)
            Q_g = self._calc_Q_g(V_hv, V_g, delta_g)
            Q_hv = self._calc_Q_hv(V_g, P_g, Q_g_lim)
            f1 = Q_g_lim - Q_g
            f2 = Q_hv_lim - Q_hv 
            return np.array([f1, f2])
        
        P_g, Q_min_0, Q_max_0 = self.get_Q_lims(V_hv, P_pu)
        X0_min = np.array([V_hv, Q_min_0])
        X0_max = np.array([V_hv, Q_max_0]) 
        sol_min = root(objective, X0_min, args=(0, P_pu)) 
        sol_max = root(objective, X0_max, args=(1, P_pu)) 
        V_g_min, Q_hv_min = sol_min.x 
        V_g_max, Q_hv_max = sol_max.x 
        
        P_g_min, Q_g_min, _ = self.get_Q_lims(V_g_min, P_pu)
        P_g_max, _, Q_g_max = self.get_Q_lims(V_g_max, P_pu)
        return P_g_min, P_g_max, Q_g_min, Q_g_max, Q_hv_min, Q_hv_max, V_g_min, V_g_max
        
    def calc_Q_at_Vg_lim(self, V_hv, V_g_lim, P_g): 
        delta = self._calc_delta_g(V_hv, V_g_lim, P_g)
        Q_g = self._calc_Q_g(V_hv, V_g_lim, delta)
        Q_g_hv = self._calc_Q_hv(V_g_lim, P_g, Q_g)
        return Q_g_hv
    
    def _calc_Q_g(self, V_hv, V_g, delta): 
        return V_g**2/self.X_T_pu - V_hv*V_g*cos(delta)/self.X_T_pu 
    
    def _calc_delta_g(self, V_hv, V_g, P_g): 
        return np.arcsin(P_g * self.X_T_pu / (V_g * V_hv)) 
    
    def _calc_Q_hv(self, V_g, P_g, Q_g): 
        return (Q_g*V_g - self.X_T_pu*(P_g**2 + Q_g**2)) / (V_g**2) 

