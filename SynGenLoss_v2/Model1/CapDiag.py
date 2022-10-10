import numpy as np 
from numpy import cos, sin, sqrt, arctan
from scipy.interpolate import interp1d 
from scipy.optimize import brenth, root
from .GenDataClass import Model1DataClass
from .components.GenSaturationModel_v1 import SaturationModel
from numba import njit 
import cmath as cm

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
    def __init__(self, gen_data: Model1DataClass, sat_data: SaturationModel): 
        self.md = gen_data
        self.sat = sat_data
        self.xq_pu = gen_data.Xq 
        self.xd_pu = gen_data.Xd 
        self.E_min = 0.1
        self.delta_max = 50*np.pi/180
        self.P_min = 1e-3
        self.P_max = 1.0
        self.I_f_max = 2.0
    
    def get_Q_lims(self, V: float, P_pu: float): 
        """Finds the Q_pu limits for a given voltage and active power- \n 
        The active power value is limited to be within P_min and P_max, considering \n 
        both set limitations and the capability diagram upper limit. \n 
        Returns P_new, Q_min, Q_max """
        Q1 = CDFuncs.calc_Q1(V, self.xd_pu, self.xq_pu, self.E_min, P_pu) 
        Q2 = CDFuncs.calc_Q2(V, self.xd_pu, self.xq_pu, P_pu, self.delta_max) 
        P_g, Q3_min, Q3_max = CDFuncs.calc_Q3(V, P_pu) 
        Q4 = CDFuncs.calc_Q4(V, self.I_f_max, self.md.Xd, self.md.Xq, self.md.Xp, self.md.Ra, self.sat.bv, self.sat.k, self.sat.Cm, self.sat.m, P_pu)         
        Q_max = np.min((Q3_max, Q4))
        if np.isnan(Q1): 
            Q_min = np.max((Q2, Q3_min))
        else: 
            Q_min = np.max((Q1, Q2, Q3_min))
        return P_g, Q_min, Q_max
