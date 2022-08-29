import numpy as np
from SynGenLoss.Model1 import SaturationModel
from SynGenLoss.Model1 import Model1DataClass
from SynGenLoss.Model1 import GeneratorLossModel 

model_data = Model1DataClass() 
model_data.standard_params(Sn_mva=103, V_nom_kV=11.0, cos_phi=0.9, If_nom_A=525.15, Ra=0.00182, Xd=1.059, Xq=0.676, Xp=0.141)
model_data.nominal_losses(V_nom=1.0, Ia_nom=5406.1, If_nom=1065, P_an_kW=186.46, P_sn_kW=89.16, P_fn_kW=173.65, P_brn_kW=2.13, P_exn_kW=15.88, 
                          P_cn_kW=211.92, P_wfn_kW=172.92, P_bn_kW=240.90)

sat_model = SaturationModel(bv=1.0, k=1.0308, Cm=0.160, m=7)

class AabjoraModel(GeneratorLossModel): 
    def __init__(self): 
        super().__init__(model_data, sat_model)

