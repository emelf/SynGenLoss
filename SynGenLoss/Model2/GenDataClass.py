from math import sqrt
from numpy import arccos, tan
class Model2DataClass: 
    """A dataclass for storing model parameters. \n 
    Usage: Call the following functions \n 
    standard_params(*params) \n 
    nominal_losses(*params) \n """
    
    def standard_params(self, Sn_mva, cos_phi, Vn_kv, f, R_a, R_fd, X_d, X_q, X_l, SG10, SG12): 
        """Sn_mva: Rated apparent power of the generator. [MVA] \n 
        cos_phi: Rated power factor of the machine [.] \n
        V_n_kV: Rated nominal voltage of the generator. [kV] \n
        R_a: Armature resistance of the generator. [pu] \n 
        R_fd: Field winding resistance in the rotor. [pu] \n
        X_d: Direct axis synchronous reactance of the generator. [pu] \n 
        X_q: Quadrature axis synchronous reactance of the generator. [pu] \n
        X_l: Leakage reactance of the generator. [pu] \n 
        SG10: Saturation coefficient [.] \n 
        SG12: Saturation coefficient [.] """
        self.Vn_kv = Vn_kv
        self.Sn_mva = Sn_mva 
        self.cos_phi = cos_phi
        self.P_nom = cos_phi 
        self.Q_nom = self.P_nom * tan(arccos(cos_phi))
        self.f = f 
        self.R_a = R_a 
        self.R_fd = R_fd
        self.X_d = X_d 
        self.X_q = X_q 
        self.X_l = X_l 
        self.X_ad_un = X_d - X_l 
        self.X_aq_un = X_q - X_l
        self.SG10 = SG10 
        self.SG12 = SG12
        
    def nominal_losses(self, V_nom, Ia_nom, If_nom, P_an_kW, P_sn_kW, P_fn_kW, P_brn_kW, 
                       P_exn_kW, P_cn_kW, P_wfn_kW, P_bn_kW):
        """V_nom: Voltage on generator terminal for given power loss values [pu] \n 
        Ia_nom: Armature current for given power loss values [pu] \n 
        If_nom: Field current for given power loss values [pu] \n 
        P_an_kW: Copper power loss from the armature [kW] \n 
        P_sn_kW: Stray power loss from the machine [kW] \n 
        P_fn_kW: Field power loss from the machine [kW] \n
        P_brn_kW: Brush power loss from the machine [kW] \n 
        P_exn_kW: Exciter power loss from the machine [kW] \n 
        P_cn_kW: Iron core power loss from the machine [kW] \n 
        P_wfn_kW: Winding and friction power loss from the machine [kW] \n 
        P_bn_kW: Bearing power loss from the machine [kW]"""
        self.V_nom = V_nom
        self.Ia_nom = Ia_nom 
        self.If_nom = If_nom 
        self.P_an = P_an_kW/(self.Sn_mva*1000)
        self.P_sn = P_sn_kW/(self.Sn_mva*1000)
        self.P_fn = P_fn_kW/(self.Sn_mva*1000)
        self.P_brn = P_brn_kW/(self.Sn_mva*1000)
        self.P_exn = P_exn_kW/(self.Sn_mva*1000)
        self.P_cn = P_cn_kW/(self.Sn_mva*1000)
        self.P_wfn = P_wfn_kW/(self.Sn_mva*1000)
        self.P_bn = P_bn_kW/(self.Sn_mva*1000)
         
    # def get_standard_params(self): 
    #     return self.Sn_mva, self.V_nom_pu, self.Ia_nom, self.If_nom, self.Ra, self.Xd, self.Xq, self.Xp
    
    # def get_nominal_losses(self): 
    #     return self.P_an, self.P_sn, self.P_fn, self.P_brn, self.P_exn, self.P_cn, self.P_wfn, self.P_bn