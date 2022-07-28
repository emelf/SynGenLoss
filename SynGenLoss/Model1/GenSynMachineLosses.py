import numpy as np
from SynGenLoss.components.GenSaturationModel_v1 import SaturationModel
from SynGenLoss.components.GenStatorLossModel_v1 import get_stator_loss
from SynGenLoss.components.GenRotorLossModel_v1 import get_rotor_loss
from SynGenLoss.components.GenConstantLossModel_v1 import get_constant_losses
from SynGenLoss.utils import Model1DataClass

class GeneratorLossModel: 
    """ Main class for the generator loss model. Requires model data and saturation model to be defined before use. """
    def __init__(self, model_data: Model1DataClass, saturationmodel: SaturationModel): 
        self.md = model_data
        self.satmodel = saturationmodel
        
    def calc_phi(self, P_el, Q_el):
        if P_el == 0 and Q_el == 0: 
            return 0
        elif P_el == 0 and not Q_el == 0: 
            return np.pi/2 * np.sign(Q_el)
        else: 
            return np.arctan(Q_el/P_el) 
    
    def calc_currents(self, P, Q, Vt): 
        """Calculates the stator and rotor currents (and load angle) based on given inputs. \n
        returns ia [pu], ifd [pu], delta [rad]"""
        ia = np.sqrt(P**2 + Q**2)/Vt
        if hasattr(ia, "__len__"): #Checks if ia is a list/array or a scalar 
            phi = np.array([self.calc_phi(P_el, Q_el) for P_el, Q_el in zip(P, Q)])
        else: 
            phi = self.calc_phi(P, Q)
        ifd, delta, _ = self.satmodel.calc_ifd(Vt, ia, phi, self.md)
        return ia, ifd, delta

    def calc_losses_pu(self, P, Q, Vt): 
        """Calculate generator losses based on P, Q, and Vt. \n
        returns a tuple of (efficiency, P_loss_stator, P_loss_rotor, P_loss_constant) in [pu]"""
        ia, ifd, _ = self.calc_currents(P, Q, Vt) 
        P_loss_stator = get_stator_loss(ia, self.md.Ia_nom, self.md.P_an, self.md.P_sn)
        P_loss_rotor = get_rotor_loss(ifd, self.md.If_nom, self.md.P_fn, self.md.P_exn, self.md.P_brn)
        P_loss_constant = get_constant_losses(Vt, self.md.V_nom, self.md.P_cn, self.md.P_wfn, self.md.P_bn)
        P_tot = P_loss_constant + P_loss_stator + P_loss_rotor
        n = P/(P + P_tot)
        return (n, P_loss_stator, P_loss_rotor, P_loss_constant)
    