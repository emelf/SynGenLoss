import numpy as np 
from numpy import cos, sin, sqrt, arctan
from scipy.interpolate import interp1d 
from scipy.optimize import root
from copy import deepcopy
from numba import njit 
import cmath as cm
from typing import Sequence

from .ModelClasses import GeneratorModel, TrafoModel

class CDFuncs: 
    @staticmethod
    @njit
    def calc_Q1(V, xd_pu, xq_pu, E_min, P_g): 
        """Minimum internal voltage limit / UEL"""     
        c = (0.5*V**2 * (xd_pu+xq_pu))/(xd_pu*xq_pu) - E_min
        R = (0.5*V**2 * (xd_pu-xq_pu))/(xd_pu*xq_pu) + E_min 
        th = np.arcsin(P_g / R)
        Q_g = R*np.cos(th) - c
        return Q_g
    
    @staticmethod
    def calc_Q2(V, xd_pu, xq_pu, P_g, delta_max): 
        """Rotor Stability Limit"""
        def objective(X, delta_max): 
            Q_g, E_q = X
            a = E_q*xq_pu / (V*(xd_pu-xq_pu))
            b = np.sqrt(a**2 + 8) 
            cos_delta = 0.25 * (b - a) #stability factor ?
            delta = np.arccos(cos_delta) 
            c = E_q*V / xd_pu * np.sin(delta) 
            d = 0.5*V**2 * (xd_pu - xq_pu)/xd_pu * xq_pu * np.sin(2*delta)
            P = c + d 
            Q = V*E_q/xd_pu*np.cos(delta_max) - V**2*(np.sin(delta_max)**2/xq_pu + np.cos(delta_max)**2/xd_pu)
            I_a = np.sqrt(P_g**2 + Q_g**2)/V 
            phi = np.arctan(Q_g / P_g)    
            E = V + 1j*cm.rect(xd_pu*I_a, -phi) # + cm.rect(ra_pu*I_a, -phi)
            E = abs(E)
            f2 = Q - Q_g 
            f3 = E - E_q
            return np.array([f2, f3])

        sol = root(objective, x0=np.array([0.5, V]), args=(delta_max))
        Q_g, E_q = sol.x
        stab_m = 0.75 #Stability marigin (comes from??)
        return Q_g*stab_m
    
    @staticmethod
    @njit
    def calc_Q3(V, P_g):
        """Stator current limit \n 
        returns (P_new, Q_min, Q_max) \n 
        equations stems from I^2 = (P_g^2 + Q_g^2)/V_g^2 -> I = 1.0 to reach stator limit. \n 
        Therefore, P_g^2 + Q_g^2 = V_g^2 """
        if P_g >= V: 
            return (V, 0.0, 0.0)
        Q_max = np.sqrt(V**2 - P_g**2)
        Q_min = -Q_max
        return (P_g, Q_min, Q_max)
    
    @staticmethod   
    def calc_Q4(V_g, I_f_max, x_d, x_q, x_p, r_a, b_v, k, Cm, m, P_g): 
        def objective(X, x_d, x_q, x_p, r_a, b_v, k, Cm, m, V_g, P, I_f): 
            e_g, e_p, delta, theta, Q = X 
            I_a = sqrt(P**2 + Q**2)/V_g
            phi = arctan(Q/P)
            f1 = (e_g - e_p)/b_v + k*(e_p + Cm*e_p**m) - I_f
            f2 = V_g*cos(delta) + r_a*I_a*cos(delta+phi) + x_d*I_a*sin(delta+phi) - e_g
            f3 = V_g*cos(theta) + r_a*I_a*cos(theta+phi) + x_p*I_a*sin(theta+phi) - e_p 
            f4 = arctan((I_a*(x_q*cos(phi) - r_a*sin(phi))) / (V_g + I_a*(r_a*cos(phi) + x_q*sin(phi)))) - delta 
            f5 = arctan((I_a*(x_p*cos(phi) - r_a*sin(phi))) / (V_g + I_a*(r_a*cos(phi) + x_p*sin(phi)))) - theta 
            return np.array([f1,f2,f3,f4,f5])
        
        X0 = np.array([V_g, 0.1*V_g, 0, 0, 0.5])
        sol = root(objective, X0, args=(x_d, x_q, x_p, r_a, b_v, k, Cm, m, V_g, P_g, I_f_max))
        Q_g = sol.x[-1]
        return Q_g

class CapabilityDiagram: 
    def __init__(self, gen: GeneratorModel): 
        self.md = deepcopy(gen.md)
        self.sat = deepcopy(gen.satmodel)
        self.E_min = 0.1
        self.delta_max = 50*np.pi/180
        self.P_min = 1e-3
        self.P_max = 1.0
        self.I_f_max = 2.0
    
    def get_Q_lims(self, V: float, P_pu: float): 
        """Finds the Q_pu limits for a given voltage and active power \n 
        The active power value is limited to be within P_min and P_max, considering \n 
        both set limitations and the capability diagram upper limit. \n 
        Returns P_new, Q_min, Q_max """
        Q1 = CDFuncs.calc_Q1(V, self.md.Xd, self.md.Xq, self.E_min, P_pu) 
        Q2 = CDFuncs.calc_Q2(V, self.md.Xd, self.md.Xq, P_pu, self.delta_max) 
        P_g, Q3_min, Q3_max = CDFuncs.calc_Q3(V, P_pu) 
        Q4 = CDFuncs.calc_Q4(V, self.I_f_max, self.md.Xd, self.md.Xq, self.md.Xp, self.md.Ra, self.sat.bv, self.sat.k, self.sat.Cm, self.sat.m, P_pu)         
        Q_max = np.min((Q3_max, Q4))
        if np.isnan(Q1): 
            Q_min = np.max((Q2, Q3_min))
        else: 
            Q_min = np.max((Q1, Q2, Q3_min))
        return P_g, Q_min, Q_max


class TrafoCapabilityDiagram(CapabilityDiagram): 
    def __init__(self, gen: GeneratorModel, trafo: TrafoModel): 
        super().__init__(gen)
        self.trafo_md = deepcopy(trafo.md)
        self.X_T_pu = self.trafo_md.X_T
        
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
    
  
class PlantCapabilityDiagram: 
    def __init__(self, gens: Sequence[GeneratorModel], trafos: Sequence[TrafoModel]): 
        self.gen_mds = [deepcopy(gen.md) for gen in gens] 
        self.sat_mds = [deepcopy(gen.satmodel) for gen in gens]
        self.trafo_mds = [deepcopy(trafo.md) for trafo in trafos]
        self.X_Ts_pu = np.array([trafo_md.X_T for trafo_md in self.trafo_mds])
        self.CDs = [CapabilityDiagram(gen) for gen in gens] 
        self.CDs_trafo = [TrafoCapabilityDiagram(gen, trafo) for gen, trafo in zip(gens, trafos)]
        self.Sn_base = np.array([md.Sn_mva for md in self.gen_mds])
        
    def get_Q_lims_plant(self, V_hv: float, P_gs_mva: Sequence[float]): 
        """Return -> (P_g, Q_hv_min, Q_hv_max, V_g_min, V_g_max) \n 
        Note: Assuming equal pu distribution. """
        Q_gs_min =  np.zeros(len(self.gen_mds))
        Q_gs_max =  np.zeros(len(self.gen_mds))
        Q_hvs_min = np.zeros(len(self.gen_mds))
        Q_hvs_max = np.zeros(len(self.gen_mds))
        V_gs_min =  np.zeros(len(self.gen_mds))
        V_gs_max =  np.zeros(len(self.gen_mds))
        P_gs_min =  np.zeros(len(self.gen_mds))
        P_gs_max =  np.zeros(len(self.gen_mds))
        
        def assign_to(idx, P_g_min, P_g_max, Q_g_min, Q_g_max, Q_hv_min, Q_hv_max, V_g_min, V_g_max): 
            Q_gs_min[idx] = Q_g_min
            Q_gs_max[idx] = Q_g_max
            Q_hvs_min[idx] = Q_hv_min
            Q_hvs_max[idx] = Q_hv_max
            V_gs_min[idx] = V_g_min
            V_gs_max[idx] = V_g_max
            P_gs_min[idx] = P_g_min
            P_gs_max[idx] = P_g_max
        
        for i, (CD_trafo, P_hv_mva) in enumerate(zip(self.CDs_trafo, P_gs_mva)): 
            S_base = CD_trafo.md.Sn_mva
            P_pu = P_hv_mva / S_base 
            P_g_min, P_g_max, Q_g_min, Q_g_max, Q_hv_min, Q_hv_max, V_g_min, V_g_max = CD_trafo.get_Q_lims_plant(V_hv, P_pu)
            assign_to(i, P_g_min*S_base, P_g_max*S_base, Q_g_min*S_base, Q_g_max*S_base, 
                      Q_hv_min*S_base, Q_hv_max*S_base, V_g_min, V_g_max)
            
        return P_gs_min, P_gs_max, Q_gs_min, Q_gs_max, Q_hvs_min, Q_hvs_max, V_gs_min, V_gs_max
    
    def get_Q_plant(self, V_hv: float, V_gs: Sequence[float], P_gs: Sequence[float]): 
        """Input: V_hv, V_gs -> Return (Q_plant, Q_gs, Q_gs_hv) \n 
        V_hv -> The voltage at the HV busbar \n 
        V_gs -> A sequence of the generator terminal voltages. """
        P_gs = P_gs / self.Sn_base
        deltas = np.array([self._calc_delta_g(V_hv, V_g, P_g, X_T) for V_g, P_g, X_T in zip(V_gs, P_gs, self.X_Ts_pu)] )
        Q_gs = np.array([self._calc_Q_g(V_hv, V_g, delta, X_T) for V_g, delta, X_T in zip(V_gs, deltas, self.X_Ts_pu)] ) * self.Sn_base
        Q_hvs = np.array([self._calc_Q_hv(V_g, P_g, Q_g, X_T) for V_g, P_g, Q_g, X_T in zip(V_gs, P_gs, Q_gs, self.X_Ts_pu)] ) * self.Sn_base
        Q_plant = Q_hvs.sum() 
        return Q_plant, Q_gs, Q_hvs
        
    def calc_Q_at_Vg_lim(self, V_hv, V_g_lims: Sequence[float], P_gs: Sequence[float]) -> Sequence[float]: 
        deltas = np.array([self._calc_delta_g(V_hv, V_g, P_g, X_T) for V_g, P_g, X_T in zip(V_g_lims, P_gs, self.X_Ts_pu)] )
        Q_gs = np.array([self._calc_Q_g(V_hv, V_g, delta, X_T) for V_g, delta, X_T in zip(V_g_lims, deltas, self.X_Ts_pu)] ) * self.Sn_base
        Q_hvs = np.array([self._calc_Q_hv(V_g, P_g, Q_g, X_T) for V_g, P_g, Q_g, X_T in zip(V_g_lims, P_gs, Q_gs, self.X_Ts_pu)] ) * self.Sn_base
        return Q_hvs
    
    def _calc_delta_g(self, V_hv, V_g, P_g, X_T): 
        return np.arcsin(P_g * X_T / (V_g * V_hv)) 
        
    def _calc_Q_g(self, V_hv, V_g, delta, X_T): 
        return V_g**2/X_T - V_hv*V_g*cos(delta)/X_T

    def _calc_Q_hv(self, V_g, P_g, Q_g, X_T): 
        return (Q_g*V_g - X_T*(P_g**2 + Q_g**2)) / (V_g**2) 
