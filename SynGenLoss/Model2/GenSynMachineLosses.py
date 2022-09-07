import numpy as np 
from numpy import cos, sin, pi 
from scipy.optimize import root 

from .components.GenConstantLossModel_v1 import get_constant_losses 
from .components.GenRotorLossModel_v1 import get_rotor_loss 
from .components.GenStatorLossModel_v1 import get_stator_loss
from .GenDataClass import Model2DataClass

class GenModel_SS: 
    """All values MUST BE UNSATURATED VALUES for the saturation model to work properly! \n
    Sn_va -> Rated power of the machine, in [MVA] \n
    Vn_kv -> Rated terminal voltage [kV] \n 
    f -> Rated frequency [Hz] \n 
    R_a -> Machine armature resistance per phase (at nominal operating point) [pu] \n 
    R_fd -> Field circuit resistance (at nominal operating point) [pu] \n 
    X_d -> Unsaturated synchronous reactance in the d-axis [pu] \n 
    X_q -> Unsaturated synchronous reactance i the q-axis [pu] \n 
    X_l -> Leakage reactance (assuming same for both axes [pu] \n 
    SG10 -> Saturation coefficient from OCC [.] \n 
    SG12 -> Saturation coefficient from OCC [.] \n """ 
    def __init__(self, model_data: Model2DataClass): 
        self.md = model_data
        self.V_base = self.md.Vn_kv
        self.S_base = self.md.Sn_mva 
        self.f = self.md.f 
        self.R_a = self.md.R_a 
        self.R_fd = self.md.R_fd
        self.X_d_un = self.md.X_d 
        self.X_q_un = self.md.X_q 
        self.X_l = self.md.X_l 
        self.X_ad_un = self.md.X_d - self.md.X_l 
        self.X_aq_un = self.md.X_q - self.md.X_l
        self.SG10 = self.md.SG10 
        self.SG12 = self.md.SG12
    
    def define_init_eqs(self, X, I_t, phi, E_t): 
        """The variable X is an array of all unknowns. All 17 equations should equal 0 in order to solve the system of equations. """
        #         0    1   2    3    4   5   6   7   8     9     10     11     12     13    14    15  
        #  X = [delta, Xd, Xq, Xad, Xaq, ed, eq, id, iq, psi_d, psi_q, sat_d, sat_q, c_sat, exp, psi_m]
        delta, X_d, X_q, X_ad, X_aq, e_d, e_q, i_d, i_q, psi_d, psi_q, sat_d, sat_q, c_sat, exp, psi_m = X
        f1 =  delta  -  (np.arctan( (X_q*I_t*cos(phi) - self.R_a*I_t*sin(phi)) / (E_t + self.R_a*I_t*cos(phi) + X_q*I_t*sin(phi)) ) ) # delta
        f2 =  X_d  -    (X_ad + self.X_l) # X_d_sat 
        f3 =  X_q  -    (X_aq + self.X_l) # X_q_sat
        f4 =  X_ad  -   (self.X_ad_un * sat_d) # X_ad_sat
        f5 =  X_aq  -   (self.X_aq_un * sat_q) # X_aq_sat 
        f6 =  e_d  -    (E_t*sin(delta)) # e_d
        f7 =  e_q  -    (E_t*cos(delta)) # e_q
        f8 =  i_d  -    (I_t*sin(phi + delta)) # i_d
        f9 =  i_q  -    (I_t*cos(phi + delta)) # i_q
        f10 = psi_d  -  (e_q + self.R_a*i_q) # psi_d
        f11 = psi_q -   (-e_d - self.R_a*i_d) # psi_q 
        f12 = sat_d -   (1.0 / (1.0 + c_sat)) # sat_d
        f13 = sat_q -   (1.0/(1 + self.X_aq_un/self.X_ad_un * c_sat)) # sat_q
        f14 = c_sat -   (self.SG10 * psi_m**exp / psi_m) # c_sat 
        f15 = exp -     (np.log(1.2 * self.SG12/self.SG10) / np.log(1.2)) # exp 
        f16 = psi_m -   (np.sqrt((psi_d + self.X_l*i_d)**2 + (psi_q + self.X_l*i_q)**2)) # psi_m
        return np.array([f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16]) 
    
    def calc_i_phi(self, P_g, Q_g, E_t): 
        I_t = np.sqrt(P_g**2 + Q_g**2)/E_t 
        phi = np.arctan(Q_g/P_g)
        return I_t, phi 

    def calc_i_e_fd(self, e_q, i_d, i_q, X_d, X_ad): 
        i_fd = (e_q + self.R_a*i_q + X_d*i_d) / X_ad
        e_fd = self.R_fd * i_fd
        return i_fd, e_fd
        
    def init_no_sat(self, P_g, Q_g, E_t):
        """returns I_t, i_d, i_q, phi, e_d, e_q, delta, psi_d, psi_q, i_fd, e_fd""" 
        I_t, phi = self.calc_i_phi(P_g, Q_g, E_t)
        
        delta = np.arctan( (self.X_q_un*I_t*cos(phi) - self.R_a*I_t*sin(phi)) / (E_t + self.R_a*I_t*cos(phi) + self.X_q_un*I_t*sin(phi)) ) 
        e_d = E_t * sin(delta)
        e_q = E_t * cos(delta) 
        i_d = I_t * sin(phi + delta) 
        i_q = I_t * cos(phi + delta) 
        psi_d = e_q + self.R_a * i_q 
        psi_q = -e_d - self.R_a * i_d 
        
        i_fd, e_fd = self.calc_i_e_fd(e_q, i_d, i_q, self.X_d_un, self.X_ad_un)

        return I_t, i_d, i_q, phi, e_d, e_q, delta, psi_d, psi_q, i_fd, e_fd
    
    def calculate(self, P_g, Q_g, E_t): 
        """
        return delta, X_d, X_q, i_d, i_q, i_fd \n
        Inputs: \n 
        P_g -> Active power generation (pu) \n 
        Q_g -> Reactive power generation (pu) \n 
        E_t -> Terminal voltage \n 
        returns: delta, Xd, Xq, id, iq, i_fd \n 
        delta -> Rotor angle \n 
        Xd -> Saturated synchronous reactance in the d axis \n 
        Xq -> Saturated synchronous reactance in the q axis \n
        Xad -> Saturated mutual synchronous reactance between the armature and d axis \n
        Xaq -> Saturated mutual synchronous reactance between the armature and q axis \n
        ed -> Internal direct axis voltage \n 
        eq -> Internal quadrature axis voltage \n 
        id -> Direct axis armature current \n 
        iq -> Quadrature axis armature current \n 
        psi_d -> Direct axis flux linkage \n 
        psi_q -> Quadrature axis flux linkage \n 
        sat_d -> Saturation coefficient in the d-axis \n 
        sat_q -> Saturation coefficient in the d-axis \n 
        c_sat -> Saturation calculation parameter \n 
        exp -> Saturation calculation parameter \n 
        psi_m -> Air-gap flux linkage \n 
        i_fd -> Field current in the rotor circuit. 
        """
        names = ["delta", "Xd", "Xq", "Xad", "Xaq", "ed", "eq", "id", "iq", "psi_d", "psi_q", "sat_d", "sat_q", "c_sat", "exp", "psi_m"]
        I_t, i_d, i_q, phi, e_d, e_q, delta, psi_d, psi_q, i_fd, e_fd = self.init_no_sat(P_g, Q_g, E_t) # Used for initial conditions
        psi_m = np.sqrt((psi_d + self.X_l*i_d)**2 + (psi_q + self.X_l*i_q)**2)
        exp = np.log(1.2*self.SG12/self.SG10)/np.log(1.2) 
        c_sat = self.SG10*psi_m**exp / psi_m 
        sat_d = 1.0 / (1.0 + c_sat) 
        sat_q = 1.0/(1.0 + self.X_aq_un/self.X_ad_un * c_sat)
        X0 = np.array([delta, self.X_d_un, self.X_q_un, self.X_ad_un, self.X_aq_un, e_d, e_q, i_d, i_q, psi_d, psi_q, sat_d, sat_q, c_sat, exp, psi_m])
        sol = root(self.define_init_eqs, X0, args=(I_t, phi, E_t), tol=1e-6) 
        delta, X_d, X_q, X_ad, X_aq, e_d, e_q, i_d, i_q, psi_d, psi_q, sat_d, sat_q, c_sat, exp, psi_m = sol.x
        i_fd, e_fd = self.calc_i_e_fd(e_q, i_d, i_q, X_d, X_ad)
        return delta, X_d, X_q, i_d, i_q, i_fd
    
    def calc_losses_pu(self, P_g, Q_g, E_t): 
        """Calculate generator losses based on P_g, Q_g, and E_t. \n
        returns a tuple of (efficiency, P_loss_stator, P_loss_rotor, P_loss_constant) in [pu]"""
        delta, X_d, X_q, i_d, i_q, i_fd = self.calculate(P_g, Q_g, E_t)
        i_a = np.sqrt(i_d**2 + i_q**2)
        P_loss_stator = get_stator_loss(i_a, self.md.Ia_nom, self.md.P_an, self.md.P_sn)
        P_loss_rotor = get_rotor_loss(i_fd, self.md.If_nom, self.md.P_fn, self.md.P_exn, self.md.P_brn)
        P_loss_constant = get_constant_losses(E_t, self.md.V_nom, self.md.P_cn, self.md.P_wfn, self.md.P_bn)
        P_tot = P_loss_constant + P_loss_stator + P_loss_rotor
        n = P_g/(P_g + P_tot)
        return (n, P_loss_stator, P_loss_rotor, P_loss_constant)
    
if __name__ == "__main__": 
    G1_data = Model2DataClass()
    G1_data.standard_params(103, 0.9, 11, 50, 0.0018, 1e-2, 1.059, 0.676, 0.08, 0.166, 0.436)
    G1_data.nominal_losses(1.0, 1.0, 1.9613396583923945, P_an_kW=186.46, P_sn_kW=89.16, P_fn_kW=173.65, P_brn_kW=2.13, 
                           P_exn_kW=15.88, P_cn_kW=211.92, P_wfn_kW=172.92, P_bn_kW=240.9)
    G1 = GenModel_SS(G1_data) 
    P_vals = [0.9, 0.675, 0.5, 0.25, 1.0, 0.75, 0.5, 0.25]
    PF_vals = [0.9, 0.9, 0.9, 0.9, 1.0, 1.0, 1.0, 1.0]
    Vt = 1.0
    
    #Calc I_fd_nom 
    Q = P_vals[0] * np.tan(np.arccos(PF_vals[0])) 
    _, _, _, _, _, i_fd0 = G1.calculate(P_vals[0], Q, Vt) 
    I_fd_nom = 1065 / i_fd0
    
    for P, PF in zip(P_vals, PF_vals): 
        Q = P * np.tan(np.arccos(PF)) 
        # delta, X_d, X_q, X_ad, X_aq, e_d, e_q, i_d, i_q, psi_d, psi_q, sat_d, sat_q, c_sat, exp, psi_m, i_fd = G1.calculate(P, Q, Vt) 
        _, _, _, _, _, i_fd = G1.calculate(P, Q, Vt) 
        eff, P_stator, P_rotor, P_const = G1.calc_losses_pu(P, Q, Vt) 
        P_loss = sum((P_stator, P_rotor, P_const))*103_000 
        
        print(f"P = {P}, Q = {Q:.2f}, i_fd = {(i_fd*I_fd_nom):.2f} A. Eff = {eff*100:.3f} %, P_loss = {P_loss:.2f} kW")
    