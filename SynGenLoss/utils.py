from math import sqrt
class Model1DataClass: 
    """A dataclass for storing model parameters. \n 
    Usage: Call the following functions \n 
    standard_params(*params) \n 
    nominal_losses(*params) \n """
    
    def standard_params(self, Sn_mva, V_nom_kV, cos_phi, If_nom_A, Ra, Xd, Xq, Xp): 
        """Sn_mva: Rated apparent power of the generator. [MVA] \n 
        V_nom_kV: Rated nominal voltage of the generator. [kV] \n
        cos_phi: Power factor at nominal operating condition. [.] \n 
        Ia_nom_A: Nominal armature current of the generator. [A] \n
        If_nom_A: Nominal field current at no-load V = 1.0 pu. [A] \n 
        Ra: Armature resistance of the generator. [pu] \n 
        Xd: Direct axis synchronous reactance of the generator. [pu] \n 
        Xq: Quadrature axis synchronous reactance of the generator. [pu] \n
        Xp: Potier reactance of the generator. [pu] \n"""
        self.Sn_mva = Sn_mva
        self.V_nom_kV = V_nom_kV
        self.cos_phi = cos_phi
        self.Ia_nom_A = Sn_mva*1e6/(sqrt(3)*V_nom_kV*1e3)
        self.If_nom_A = If_nom_A
        self.Ra, self.Xd, self.Xq, self.Xp = Ra, Xd, Xq, Xp
        
    def nominal_losses(self, V_nom, Ia_nom, If_nom, P_an_kW, P_sn_kW, P_fn_kW, P_brn_kW, P_exn_kW, P_cn_kW, P_wfn_kW, P_bn_kW):
        """V_nom: Voltage on generator terminal for given power loss values [pu] \n 
        Ia_nom: Armature current for given power loss values [A] \n 
        If_nom: Field current for given power loss values [A] \n 
        P_an_kW: Copper power loss from the armature [kW] \n 
        P_sn_kW: Stray power loss from the machine [kW] \n 
        P_fn_kW: Field power loss from the machine [kW] \n
        P_brn_kW: Brush power loss from the machine [kW] \n 
        P_exn_kW: Exciter power loss from the machine [kW] \n 
        P_cn_kW: Iron core power loss from the machine [kW] \n 
        P_wfn_kW: Winding and friction power loss from the machine [kW] \n 
        P_bn_kW: Bearing power loss from the machine [kW]"""
        self.V_nom = V_nom
        self.Ia_nom = Ia_nom / self.Ia_nom_A
        self.If_nom = If_nom / self.If_nom_A
        self.P_an = P_an_kW/(self.Sn_mva*1000)
        self.P_sn = P_sn_kW/(self.Sn_mva*1000)
        self.P_fn = P_fn_kW/(self.Sn_mva*1000)
        self.P_brn = P_brn_kW/(self.Sn_mva*1000)
        self.P_exn = P_exn_kW/(self.Sn_mva*1000)
        self.P_cn = P_cn_kW/(self.Sn_mva*1000)
        self.P_wfn = P_wfn_kW/(self.Sn_mva*1000)
        self.P_bn = P_bn_kW/(self.Sn_mva*1000)
        
        
    def get_standard_params(self): 
        return self.Sn_mva, self.V_nom_pu, self.cos_phi, self.Ia_nom, self.If_nom, self.Ra, self.Xd, self.Xq, self.Xp
    
    def get_nominal_losses(self): 
        return self.P_an, self.P_sn, self.P_fn, self.P_brn, self.P_exn, self.P_cn, self.P_wfn, self.P_bn
    
