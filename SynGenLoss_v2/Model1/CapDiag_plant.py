import numpy as np 
from numpy import cos, sin, sqrt, arctan

from scipy.optimize import brenth, root
from .GenDataClass import Model1DataClass
from .components.GenSaturationModel_v1 import SaturationModel

import cmath as cm
from .CapDiag import CapabilityDiagram 
from .CapDiag_trafo import TrafoCapabilityDiagram
from copy import deepcopy
from typing import Sequence

class PlantCapabilityDiagram: 
    def __init__(self, gen_data: Sequence[Model1DataClass], sat_data: Sequence[SaturationModel], X_Ts_pu: Sequence[float]): 
        self.mds = [deepcopy(gen_d) for gen_d in gen_data] 
        self.sat = [deepcopy(sat_d) for sat_d in sat_data]
        self.X_Ts_pu = X_Ts_pu
        self.CDs = [CapabilityDiagram(gen_d, sat_d) for gen_d, sat_d in zip(self.mds, self.sat)] 
        self.CDs_trafo = [TrafoCapabilityDiagram(gen_d, sat_d, X_T) for gen_d, sat_d, X_T in zip(self.mds, self.sat, self.X_Ts_pu)]
        
    def get_Q_lims_plant(self, V_hv: float, P_gs: Sequence[float]): 
        """Return -> (P_g, Q_hv_min, Q_hv_max, V_g_min, V_g_max) \n 
        Note: Assuming equal pu distribution. """
        Q_gs_min =  np.zeros(len(self.mds))
        Q_gs_max =  np.zeros(len(self.mds))
        Q_hvs_min = np.zeros(len(self.mds))
        Q_hvs_max = np.zeros(len(self.mds))
        V_gs_min =  np.zeros(len(self.mds))
        V_gs_max =  np.zeros(len(self.mds))
        P_gs_min =  np.zeros(len(self.mds))
        P_gs_max =  np.zeros(len(self.mds))
        
        def assign_to(idx, P_g_min, P_g_max, Q_g_min, Q_g_max, Q_hv_min, Q_hv_max, V_g_min, V_g_max): 
            Q_gs_min[idx] = Q_g_min
            Q_gs_max[idx] = Q_g_max
            Q_hvs_min[idx] = Q_hv_min
            Q_hvs_max[idx] = Q_hv_max
            V_gs_min[idx] = V_g_min
            V_gs_max[idx] = V_g_max
            P_gs_min[idx] = P_g_min
            P_gs_max[idx] = P_g_max
        
        for i, (CD_trafo, P_pu) in enumerate(zip(self.CDs_trafo, P_gs)):   
            P_g_min, P_g_max, Q_g_min, Q_g_max, Q_hv_min, Q_hv_max, V_g_min, V_g_max = CD_trafo.get_Q_lims_plant(V_hv, P_pu)
            assign_to(i, P_g_min, P_g_max, Q_g_min, Q_g_max, Q_hv_min, Q_hv_max, V_g_min, V_g_max)
            
        return P_gs_min, P_gs_max, Q_gs_min, Q_gs_max, Q_hvs_min, Q_hvs_max, V_gs_min, V_gs_max
    
    def get_Q_plant(self, V_hv: float, V_gs: Sequence[float], P_gs: Sequence[float]): 
        """Input: V_hv, V_gs -> Return (Q_plant, Q_gs, Q_gs_hv) \n 
        V_hv -> The voltage at the HV busbar \n 
        V_gs -> A sequence of the generator terminal voltages. """
        deltas = np.array([self._calc_delta_g(V_hv, V_g, P_g, X_T) for V_g, P_g, X_T in zip(V_gs, P_gs, self.X_Ts_pu)] )
        Q_gs = np.array([self._calc_Q_g(V_hv, V_g, delta, X_T) for V_g, delta, X_T in zip(V_gs, deltas, self.X_Ts_pu)] )
        Q_hvs = np.array([self._calc_Q_hv(V_g, P_g, Q_g, X_T) for V_g, P_g, Q_g, X_T in zip(V_gs, P_gs, Q_gs, self.X_Ts_pu)] )
        Q_plant = Q_hvs.sum() 
        return Q_plant, Q_gs, Q_hvs
        
    def calc_Q_at_Vg_lim(self, V_hv, V_g_lims: Sequence[float], P_gs: Sequence[float]) -> Sequence[float]: 
        deltas = np.array([self._calc_delta_g(V_hv, V_g, P_g, X_T) for V_g, P_g, X_T in zip(V_g_lims, P_gs, self.X_Ts_pu)] )
        Q_gs = np.array([self._calc_Q_g(V_hv, V_g, delta, X_T) for V_g, delta, X_T in zip(V_g_lims, deltas, self.X_Ts_pu)] )
        Q_hvs = np.array([self._calc_Q_hv(V_g, P_g, Q_g, X_T) for V_g, P_g, Q_g, X_T in zip(V_g_lims, P_gs, Q_gs, self.X_Ts_pu)] )
        return Q_hvs
    
    def _calc_delta_g(self, V_hv, V_g, P_g, X_T): 
        return np.arcsin(P_g * X_T / (V_g * V_hv)) 
        
    def _calc_Q_g(self, V_hv, V_g, delta, X_T): 
        return V_g**2/X_T - V_hv*V_g*cos(delta)/X_T

    def _calc_Q_hv(self, V_g, P_g, Q_g, X_T): 
        return (Q_g*V_g - X_T*(P_g**2 + Q_g**2)) / (V_g**2) 

