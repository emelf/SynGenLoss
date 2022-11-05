from math import sqrt, atan, pi
import cmath as cm
from typing import Sequence

class GenDataClass: 
    """A dataclass for storing generator model parameters. \n 
    Usage: Call the following functions \n 
    standard_params(*params) \n 
    nominal_losses(*params) \n """
    
    def standard_params(self, Sn_mva:float, V_nom_kV: float, cos_phi: float, If_nom_A: float, 
                        Ra: float, Xd: float, Xq: float, Xp: float) -> None: 
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
        
    def nominal_losses(self, V_nom:float, Ia_nom:float, If_nom:float, P_an_kW:float, P_sn_kW:float, 
                       P_fn_kW:float, P_brn_kW:float, P_exn_kW:float, P_cn_kW:float, P_wfn_kW:float, 
                       P_bn_kW:float) -> None:
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
        
    def get_standard_params(self) -> Sequence[float]: 
        return self.Sn_mva, self.V_nom_pu, self.cos_phi, self.Ia_nom, self.If_nom, self.Ra, self.Xd, self.Xq, self.Xp
    
    def get_nominal_losses(self) -> Sequence[float]: 
        return self.P_an, self.P_sn, self.P_fn, self.P_brn, self.P_exn, self.P_cn, self.P_wfn, self.P_bn


class TrafoDataClass: 
    """A dataclass for storing transformer model parameters. \n 
    Usage: Call the following functions \n 
    define_params(*params) \n """
    def __init__(self, Sn_mva: float, V_hv_kv: float, V_lv_kv: float, V_SHC_pu: float, P_Cu_pu: float, I_E: float, P_Fe_pu: float) -> None:
        self.define_params(Sn_mva, V_hv_kv, V_lv_kv, V_SHC_pu, P_Cu_pu, I_E, P_Fe_pu)
    
    
    def define_params(self, Sn_mva: float, V_hv_kv: float, V_lv_kv: float, V_SHC_pu: float, P_Cu_pu: float, I_E:float, P_Fe_pu: float) -> None: 
        """Sn_mva: Rated apparent power [MVA] \n 
        V_hv_kv: Rated voltage at the HV side [kV] \n
        V_lv_kv: Rated voltage at the LV side [kV] \n
        V_SHC_pu: Short circuit voltage [pu] \n
        P_Cu_pu: Nominal copper losses of the machine (equivalent to R_T) [pu] \n 
        P_Fe_pu: No-load losses [pu] \n 
        Source: Ch. 3.2.1 in "Power System Dynamics: Stability and Control", Jan Machowski
        """
        self.Sn_mva = Sn_mva
        self.V_hv_kv = V_hv_kv
        self.V_lv_kv = V_lv_kv
        self.V_SHC_pu = V_SHC_pu
        self.P_Cu_pu = P_Cu_pu
        self.I_E = I_E
        self.P_Fe_pu = P_Fe_pu
        
        # Calculation parameters: 
        self.G_Fe = self.P_Fe_pu 
        self.B_mu = sqrt(self.I_E**2 - self.G_Fe**2)
        self.Y_E = self.G_Fe - 1j*self.B_mu
        
        self.Z_T = V_SHC_pu
        self.R_T = P_Cu_pu 
        self.X_T = sqrt(self.Z_T**2 - self.R_T**2)
        self.Z_angle = atan(self.X_T / self.R_T) if self.R_T != 0.0 else pi/2
        self.Z_T = cm.rect(self.Z_T, self.Z_angle)
        
    def get_input_params(self) -> Sequence[float]:
        """Return Sn_mva, V_hv_kv, V_lv_kv, V_SHC_pu, P_Cu_pu, P_Fe_pu"""
        return self.Sn_mva, self.V_hv_kv, self.V_lv_kv, self.V_SHC_pu, self.P_Cu_pu, self.P_Fe_pu
    
    def get_calc_params(self) -> Sequence[float]: 
        """Return G_Fe, Y_E, B_mu, Z_T, R_T, X_T \n 
        All quantities in pu. """
        return self.G_Fe, self.B_mu, self.Y_E, self.B_mu, self.Z_T, self.R_T, self.X_T


class LineDataClass: 
    """A dataclass for storing line/cable model parameters. Assuming the pi model.\n 
    Usage: Call the following functions \n 
    define_params(*params) \n """
    
    def define_params(self, Sn_mva:float, V_kv: float, r: float, x: float, b: float, length: float) -> None: 
        """Sn_mva: Rated apparent power [MVA] \n 
        V_kv: Rated voltage \n
        r: resistance [Ohm/km] \n
        x: reactance [Ohm/km] \n 
        b: susceptance [uS/km] \n 
        length: Line length [km] \n
        Source: Ch. 3.2.1 in "Power System Dynamics: Stability and Control", Jan Machowski
        """
        self.Sn_mva = Sn_mva
        self.V_kv = V_kv
        self.r = r 
        self.x = x 
        self.b = b 
        self.length = length 
        
        # Calculation parameters: 
        self.Z_base = self.V_kv**2/self.Sn_mva
        
    def get_params(self) -> Sequence[float]:
        """Return V_kv, r, x, b, length"""
        return self.V_kv, self.r, self.x, self.b, self.length
    
