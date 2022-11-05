import numpy as np
import cmath as cm
from typing import Sequence, Tuple
from .components.GenConstantLossModel_v1 import get_constant_losses 
from .components.GenRotorLossModel_v1 import get_rotor_loss 
from .components.GenSaturationModel_v1 import SaturationModel 
from .components.GenStatorLossModel_v1 import get_stator_loss
from .DataClasses import GenDataClass, TrafoDataClass, LineDataClass

class GeneratorModel: 
    """ Main class for the generator loss model. Requires model data and saturation model to be defined before use. """
    def __init__(self, model_data: GenDataClass, saturationmodel: SaturationModel) -> None: 
        self.md = model_data
        self.satmodel = saturationmodel
        
    def _calc_phi(self, P_el:float, Q_el:float) -> float:
        if P_el == 0 and Q_el == 0: 
            return 0
        elif P_el == 0 and not Q_el == 0: 
            return np.pi/2 * np.sign(Q_el)
        else: 
            return np.arctan(Q_el/P_el) 
    
    def calc_currents(self, P_pu: float, Q_pu: float, V_t: float) -> Tuple[float, float, float]: 
        """Calculates the stator and rotor currents (and load angle) based on given inputs. \n
        returns ia [pu], ifd [pu], delta [rad]"""
        ia = np.sqrt(P_pu**2 + Q_pu**2)/V_t
        if hasattr(ia, "__len__"): #Checks if ia is a list/array or a scalar 
            phi = np.array([self._calc_phi(P_el, Q_el) for P_el, Q_el in zip(P_pu, Q_pu)])
        else: 
            phi = self._calc_phi(P_pu, Q_pu)
        ifd, delta, _ = self.satmodel.calc_ifd(V_t, ia, phi, self.md)
        return ia, ifd, delta

    def calc_losses_pu(self, P_pu: float, Q_pu: float, V_t: float) -> Tuple[float, float, float, float]: 
        """Calculate generator losses based on P, Q, and Vt. \n
        returns a tuple of (efficiency, P_loss_stator, P_loss_rotor, P_loss_constant) in [pu]"""
        ia, ifd, _ = self.calc_currents(P_pu, Q_pu, V_t) 
        P_loss_stator = get_stator_loss(ia, self.md.Ia_nom, self.md.P_an, self.md.P_sn)
        P_loss_rotor = get_rotor_loss(ifd, self.md.If_nom, self.md.P_fn, self.md.P_exn, self.md.P_brn)
        P_loss_constant = get_constant_losses(V_t, self.md.V_nom, self.md.P_cn, self.md.P_wfn, self.md.P_bn)
        P_tot = P_loss_constant + P_loss_stator + P_loss_rotor
        n = P_pu/(P_pu + P_tot)
        return (n, P_loss_stator, P_loss_rotor, P_loss_constant)
    
    
class BranchModel: 
    def __init__(self, model_data, Z, Y) -> None: 
        self.md = model_data
        self._A = self._D = 1 + Z*Y/2
        self._B = Z 
        self._C = Y*(1 + Z*Y/4)
        self._mat = np.array([[self._A, self._B], [self._C, self._D]])
        
    def calc_PQV_sending(self, P_r_mva: float, Q_r_mva: float, V_r_pu: float) -> Tuple[float, float, float, float]: 
        """Returns (P_s_mva, Q_s_mva, V_s_pu, delta_sr)"""
        P = P_r_mva / self.md.Sn_mva # Convert to pu 
        Q = Q_r_mva / self.md.Sn_mva 
        I_r = (P - 1j*Q)/V_r_pu 
        V_s, I_s = self._mat @ np.array([V_r_pu, I_r])
        S_send = V_s * I_s.conjugate() * self.md.Sn_mva
        return S_send.real, S_send.imag, abs(V_s), np.angle(V_s)
    
    def calc_PQV_receiving(self, P_s_mva: float, Q_s_mva: float, V_s_pu: float) -> Tuple[float, float, float, float]: 
        """Returns (P_r_mva, Q_r_mva, V_r_pu, delta_sr)"""
        P = P_s_mva / self.md.Sn_mva # Convert to pu 
        Q = Q_s_mva / self.md.Sn_mva 
        I_s = (P - 1j*Q)/V_s_pu 
        V_r, I_r = np.linalg.inv(self._mat) @ np.array([V_s_pu, I_s])
        S_r = V_r * I_s.conjugate() * self.md.Sn_mva
        return S_r.real, S_r.imag, abs(V_r), np.angle(V_r) 
    
    
class TrafoModel(BranchModel): 
    def __init__(self, model_data: TrafoDataClass) -> None: 
        md = model_data 
        super().__init__(model_data, model_data.Z_T, model_data.Y_E)
        
        
class LineModel(BranchModel): 
    def __init__(self, model_data: LineDataClass) -> None: 
        self.md = model_data 
        self.R = self.md.r * self.md.length
        self.X = self.md.x * self.md.length
        self.Z = self.R + 1j*self.X
        self.Y = 1e-6j*self.md.b * self.md.length
        super().__init__(model_data, self.Z, self.Y)
        
        
class PowerPlantModel: 
    """Assuming each generator is connected through its own step-up transformer. """
    def __init__(self, gen_models: Sequence[GeneratorModel], trafo_models: Sequence[TrafoModel]) -> None: 
        self.gen_models = gen_models 
        self.trafo_models = trafo_models  
        
    def calc_gen_terminals(self, P_gs_hv_mva: Sequence[float], Q_gs_hv_mva: Sequence[float], V_hv_pu: float
                           ) -> Tuple[Sequence[float], Sequence[float], Sequence[float], Sequence[float]]:
        """Return P_gs, Q_gs, V_gs, delta_gs"""
        P_gs = np.zeros(len(P_gs_hv_mva))
        Q_gs = np.zeros(len(P_gs_hv_mva))
        V_gs = np.zeros(len(P_gs_hv_mva))
        delta_gs = np.zeros(len(P_gs_hv_mva))
        for i, (trafo, P, Q) in enumerate(zip(self.trafo_models, P_gs_hv_mva, Q_gs_hv_mva)): 
            P_gs[i], Q_gs[i], V_gs[i], delta_gs[i] = trafo.calc_PQV_sending(P, Q, V_hv_pu) 
        return P_gs, Q_gs, V_gs, delta_gs
    
    def calc_plant_losses(self, P_gs_hv_mva: Sequence[float], Q_gs_hv_mva: Sequence[float], V_hv_pu: float
                          ) -> Tuple[Sequence[float], Sequence[float], Sequence[float], Sequence[float]]:
        """Return P_gs_loss, P_ts_loss"""
        P_gs, Q_gs, V_gs, delta_gs = self.calc_gen_terminals(P_gs_hv_mva, Q_gs_hv_mva, V_hv_pu)
        P_ts_loss = np.abs(P_gs_hv_mva - P_gs) 
        P_gs_loss = np.zeros(len(P_gs_hv_mva))
        for i, (gen, P, Q, V) in enumerate(zip(self.gen_models, P_gs, Q_gs, V_gs)): 
            S_n = gen.md.Sn_mva
            P = P / S_n 
            Q = Q / S_n
            n, P_loss_stator, P_loss_rotor, P_loss_constant = gen.calc_losses_pu(P, Q, V)
            P_gs_loss[i] = (P_loss_stator + P_loss_rotor + P_loss_constant)*S_n
        return P_gs_loss, P_ts_loss 