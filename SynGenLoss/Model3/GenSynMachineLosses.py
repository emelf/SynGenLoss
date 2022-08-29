# import numpy as np 
# from numpy import cos, sin, pi 
# from scipy.optimize import root 

# class GenModel_dyn: 
#     """All values MUST BE UNSATURATED VALUES for the saturation model to work properly! \n
#     Sn_va -> Rated power of the machine, in [MVA] \n
#     Vn_kv -> Rated terminal voltage [kV] \n 
#     f -> Rated frequency [Hz] \n 
#     R_a -> Machine armature resistance per phase (at nominal operating point) [pu] \n 
#     R_fd -> Field circuit resistance (at nominal operating point) [pu] \n 
#     R_1d -> Resistance in amortisseur winding nr. 1 in the d-axis [pu] \n 
#     R_1q -> Resistance in amortisseur winding nr. 1 in the q-axis [pu] \n 
#     R_2q -> Resistance in amortisseur winding nr. 2 in the q-axis [pu] \n 
#     X_d -> Unsaturated synchronous reactance in the d-axis [pu] \n 
#     X_q -> Unsaturated synchronous reactance i the q-axis [pu] \n 
#     X_l -> Leakage reactance (assuming same for both axes [pu] \n 
#     SG10 -> Saturation coefficient from OCC [.] \n 
#     SG12 -> Saturation coefficient from OCC [.] \n 
#     H -> Machine inertia related to the apparent power \n 
#     K_D -> Damper coefficient \n 
#     The machine has 7 states, all represeted by flux linkages. The voltages and currents are calculated from initial conditions, where P, Q, and V are specified. 
#     Flux linkages are represented by differential equations, and voltages and currents by algebraic equations. """ 
#     def __init__(self, Sn_mva, Vn_kv, f, R_a, R_fd, R_1d, R_1q, R_2q, X_d, X_q, X_l, , SG10, SG12, H, K_D): 
#         self.V_base = Vn_kv
#         self.S_base = Sn_mva 
#         self.w_base = 2*np.pi*f
#         self.f = f 
#         self.R_a = R_a 
#         self.R_fd = R_fd
#         self.R_1d = R_1d 
#         self.R_1q = R_1q 
#         self.R_2q = R_2q 
#         self.X_d_un = X_d 
#         self.X_q_un = X_q 
#         self.X_l = X_l 
#         self.X_ad_un = X_d - X_l 
#         self.X_aq_un = X_q - X_l
#         self.SG10 = SG10 
#         self.SG12 = SG12
#         self.H = H 
#         self.K_D = K_D
        
#         self.omega_0 = 1.0
    
#     def f(self, X, y, u): 
#         """ Returns the differential equations for the machine state. \n 
#         X -> State variables = [psi_d, psi_q, psi_0, psi_fd, psi_1d, psi_1q, psi_2q, omega, delta] -> 9 states \n 
#         y -> Algebraic variables = [e_d, e_q, e_0, e_fd, i_d, i_q, i_0, i_1d, i_1q, i_2q, X_d, X_q, X_ad, X_aq, T_e, P_g, Q_g, T_t] -> 18 variables \n 
#         u -> Input variables = [i_fd, T_m] -> 2 inputs"""
#         [psi_d, psi_q, psi_0, psi_fd, psi_1d, psi_1q, psi_2q, omega, delta] = X
#         [e_d, e_q, e_0, e_fd, i_d, i_q, i_0, i_1d, i_1q, i_2q, X_d, X_q, X_ad, X_aq, T_e, P_g, Q_g, T_t] = y
#         [i_fd, T_m] = u
        
#         #Stator voltage equations: 
#         f1 = e_d + psi_q*omega + self.R_a*i_d # diff psi_d 
#         f2 = e_q - psi_d*omega + self.R_a*i_q # diff psi_q 
#         f3 = e_0 + self.R_a*i_0 # diff psi_0 
        
#         #Rotor voltgate equations: 
#         f4 =  e_fd - self.R_fd*i_fd # diff psi_fd
#         f5 = -self.R_1d*i_1d # diff psi_1d
#         f6 = -self.R_1q*i_1q # diff psi_1q 
#         f7 = -self.R_2q*i_2q # diff psi_2q 
        
#         #Equations of the rotor motion 
#         f8 = (T_m - T_e - self.K_D*(omega-self.omega_0)/self.omega_0) / (2*self.H/self.omega_0) # diff omega
#         f9 = omega - self.omega_0 # diff delta
        
#         f = np.array([f1, f2, f3, f4, f5, f6, f7, f8, f9]) 
#         return f
    
#     def define_init_eqs(self, X, I_t, phi, E_t): 
#         """The variable X is an array of all unknowns. All 16 equations should equal 0 in order to solve the system of equations. """
#         #         0    1   2    3    4   5   6   7   8     9     10     11     12     13    14    15
#         #  X = [delta, Xd, Xq, Xad, Xaq, ed, eq, id, iq, psi_d, psi_q, sat_d, sat_q, c_sat, exp, psi_m]
#         f1 =  X[0]  - (np.arctan( (X[2]*I_t*cos(phi) - self.R_a*I_t*sin(phi)) / (E_t + self.R_a*I_t*cos(phi) + X[2]*I_t*sin(phi)) ) ) # delta
#         f2 =  X[1]  - (X[3] + self.X_l) # X_d_sat 
#         f3 =  X[2]  - (X[4] + self.X_l) # X_q_sat
#         f4 =  X[3]  - (self.X_ad_un * X[11]) # X_ad_sat
#         f5 =  X[4]  - (self.X_aq_un * X[12]) # X_aq_sat 
#         f6 =  X[5]  - (E_t * sin(X[0])) # e_d
#         f7 =  X[6]  - (E_t * cos(X[0])) # e_q
#         f8 =  X[7]  - (I_t*sin(phi + X[0])) # i_d
#         f9 =  X[8]  - (I_t*cos(phi + X[0])) # i_q
#         f10 = X[9]  - (X[6] + self.R_a*X[8]) # psi_d
#         f11 = X[10] - (-X[5] - self.R_a*X[7]) # psi_q 
#         f12 = X[11] - (1.0 / (1.0 + X[13])) # sat_d
#         f13 = X[12] - (1.0/(1 + self.X_aq_un/self.X_ad_un * X[13])) # sat_q
#         f14 = X[13] - (self.SG10 * X[15]**X[14] / X[15]) # c_sat 
#         f15 = X[14] - (np.log(1.2 * self.SG12/self.SG10) / np.log(1.2)) # exp 
#         f16 = X[15] - (np.sqrt((X[9] + self.X_l*X[7])**2 + (X[10] + self.X_l*X[8])**2)) # psi_m
#         return np.array([f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16]) 
    
#     def calc_i_phi(self, P_g, Q_g, E_t): 
#         I_t = np.sqrt(P_g**2 + Q_g**2)/E_t 
#         phi = np.arccos(P_g / (E_t * I_t)) 
#         return I_t, phi 

#     def calc_i_e_fd(self, e_q, i_d, i_q, X_d, X_ad): 
#         i_fd = (e_q + self.R_a*i_q + X_d*i_d) / X_ad
#         e_fd = self.R_fd * i_fd
#         return i_fd, e_fd
        
#     def init_no_sat(self, P_g, Q_g, E_t):
#         """returns """ 
#         I_t, phi = self.calc_i_phi(P_g, Q_g, E_t)
        
#         delta = np.arctan( (self.X_q_un*I_t*cos(phi) - self.R_a*I_t*sin(phi)) / (E_t + self.R_a*I_t*cos(phi) + self.X_q_un*I_t*sin(phi)) ) 
#         e_d = E_t * sin(delta)
#         e_q = E_t * cos(delta) 
#         i_d = I_t * sin(phi + delta) 
#         i_q = I_t * cos(phi + delta) 
#         psi_d = e_q + self.R_a * i_q 
#         psi_q = -e_d - self.R_a * i_d 
        
#         i_fd, e_fd = self.calc_i_e_fd(e_q, i_d, i_q, self.X_d_un, self.X_ad_un)

#         return I_t, i_d, i_q, phi, e_d, e_q, delta, psi_d, psi_q, i_fd, e_fd
    
#     def init(self, P_g, Q_g, E_t): 
#         #         0    1   2    3    4   5   6   7   8     9     10     11     12     13    14    15
#         #  X = [delta, Xd, Xq, Xad, Xaq, ed, eq, id, iq, psi_d, psi_q, sat_d, sat_q, c_sat, exp, psi_m]
#         names = ["delta", "Xd", "Xq", "Xad", "Xaq", "ed", "eq", "id", "iq", "psi_d", "psi_q", "sat_d", "sat_q", "c_sat", "exp", "psi_m"]
#         I_t, i_d, i_q, phi, e_d, e_q, delta, psi_d, psi_q, i_fd, e_fd = self.init_no_sat(P_g, Q_g, E_t) # Used for initial conditions
#         X0 = np.array([delta, self.X_d_un, self.X_q_un, self.X_ad_un, self.X_aq_un, e_d, e_q, i_d, i_q, psi_d, psi_q, 0.5, 0.5, 0.5, 0.5, 0.5])
#         sol = root(self.define_init_eqs, X0, args=(I_t, phi, E_t), tol=1e-6) 
#         delta, X_d, X_q, X_ad, X_aq, e_d, e_q, i_d, i_q, psi_d, psi_q, sat_d, sat_q, c_sat, exp, psi_m = sol.x
        
#         #Initialize the rest of the machine states: 
#         i_fd, e_fd = self.calc_i_e_fd(e_q, i_d, i_q, X_d, X_ad)
#         psi_fd = (X_ad + self.X_)

#         return sol.x, i_fd, e_fd