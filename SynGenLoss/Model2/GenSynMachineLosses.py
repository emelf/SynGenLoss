import numpy as np 
from numpy import cos, sin, pi 
from scipy.optimize import root 

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
    def __init__(self, Sn_mva, Vn_kv, f, R_a, R_fd, X_d, X_q, X_l, SG10, SG12): 
        self.V_base = Vn_kv
        self.S_base = Sn_mva 
        self.f = f 
        self.R_a = R_a 
        self.R_fd = R_fd
        self.X_d_un = X_d 
        self.X_q_un = X_q 
        self.X_l = X_l 
        self.X_ad_un = X_d - X_l 
        self.X_aq_un = X_q - X_l
        self.SG10 = SG10 
        self.SG12 = SG12
    
    def define_init_eqs(self, X, I_t, phi, E_t): 
        """The variable X is an array of all unknowns. All 17 equations should equal 0 in order to solve the system of equations. """
        #         0    1   2    3    4   5   6   7   8     9     10     11     12     13    14    15  
        #  X = [delta, Xd, Xq, Xad, Xaq, ed, eq, id, iq, psi_d, psi_q, sat_d, sat_q, c_sat, exp, psi_m]
        delta, X_d, X_q, X_ad, X_aq, e_d, e_q, i_d, i_q, psi_d, psi_q, sat_d, sat_q, c_sat, exp, psi_m = X
        f1 =  delta  - (np.arctan( (X_q*I_t*cos(phi) - self.R_a*I_t*sin(phi)) / (E_t + self.R_a*I_t*cos(phi) + X_q*I_t*sin(phi)) ) ) # delta
        f2 =  X_d  - (X_ad + self.X_l) # X_d_sat 
        f3 =  X_q  - (X_aq + self.X_l) # X_q_sat
        f4 =  X_ad  - (self.X_ad_un * sat_d) # X_ad_sat
        f5 =  X_aq  - (self.X_aq_un * sat_q) # X_aq_sat 
        f6 =  e_d  - (E_t*sin(delta)) # e_d
        f7 =  e_q  - (E_t*cos(delta)) # e_q
        f8 =  i_d  - (I_t*sin(phi + delta)) # i_d
        f9 =  i_q  - (I_t*cos(phi + delta)) # i_q
        f10 = psi_d  - (e_q + self.R_a*i_q) # psi_d
        f11 = psi_q - (-e_d - self.R_a*i_d) # psi_q 
        f12 = sat_d - (1.0 / (1.0 + c_sat)) # sat_d
        f13 = sat_q - (1.0/(1 + self.X_aq_un/self.X_ad_un * c_sat)) # sat_q
        f14 = c_sat - (self.SG10 * psi_m**exp / psi_m) # c_sat 
        f15 = exp - (np.log(1.2 * self.SG12/self.SG10) / np.log(1.2)) # exp 
        f16 = psi_m - (np.sqrt((psi_d + self.X_l*i_d)**2 + (psi_q + self.X_l*i_q)**2)) # psi_m
        return np.array([f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16]) 
    
    def calc_i_phi(self, P_g, Q_g, E_t): 
        I_t = np.sqrt(P_g**2 + Q_g**2)/E_t 
        phi = np.arccos(P_g / (E_t * I_t)) 
        return I_t, phi 

    def calc_i_e_fd(self, e_q, i_d, i_q, X_d, X_ad): 
        i_fd = (e_q + self.R_a*i_q + X_d*i_d) / X_ad
        e_fd = self.R_fd * i_fd
        return i_fd, e_fd
        
    def init_no_sat(self, P_g, Q_g, E_t):
        """returns """ 
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
        """Inputs: \n 
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
        X0 = np.array([delta, self.X_d_un, self.X_q_un, self.X_ad_un, self.X_aq_un, e_d, e_q, i_d, i_q, psi_d, psi_q, 0.5, 0.5, 0.5, 0.5, 0.5])
        sol = root(self.define_init_eqs, X0, args=(I_t, phi, E_t), tol=1e-6) 
        delta, X_d, X_q, X_ad, X_aq, e_d, e_q, i_d, i_q, psi_d, psi_q, sat_d, sat_q, c_sat, exp, psi_m = sol.x
        i_fd, e_fd = self.calc_i_e_fd(e_q, i_d, i_q, X_d, X_ad)
        return delta, X_d, X_q, i_d, i_q, i_fd
    
if __name__ == "__main__": 
    G1 = GenModel_SS(103, 11, 50, 0.0018, 1e-2, 1.059, 0.676, 0.08, 0.1659776903705144, 0.4356246146911903) 
    P_vals = [0.9, 0.675, 0.5, 0.25, 1.0, 0.75, 0.5, 0.25]
    PF_vals = [0.9, 0.9, 0.9, 0.9, 1.0, 1.0, 1.0, 1.0]
    Vt = 1.0
    
    for P, PF in zip(P_vals, PF_vals): 
        Q = P * np.tan(np.arccos(PF)) 
        # delta, X_d, X_q, X_ad, X_aq, e_d, e_q, i_d, i_q, psi_d, psi_q, sat_d, sat_q, c_sat, exp, psi_m, i_fd = G1.calculate(P, Q, Vt) 
        _, _, _, _, _, i_fd = G1.calculate(P, Q, Vt) 
        print(f"P = {P}, Q = {Q}, V = {Vt}, i_fd = {i_fd*525.15}")
    