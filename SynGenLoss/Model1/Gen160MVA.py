import numpy as np
from .components.GenSaturationModel_v1 import SaturationModel
from .GenDataClass import Model1DataClass
from .GenSynMachineLosses import GeneratorLossModel 

model_data = Model1DataClass() 
model_data.standard_params(Sn_mva=160, V_nom_kV=15.0, cos_phi=0.95, If_nom_A=594.0, Ra=0.002322, Xd=0.8, Xq=0.6, Xp=0.18)
model_data.nominal_losses(V_nom=1.0, Ia_nom=5406, If_nom=1047.0, P_an_kW=327.05, P_sn_kW=237.07, P_fn_kW=477.81, P_brn_kW=5.3, P_exn_kW=33.96, 
                          P_cn_kW=539.87, P_wfn_kW=710.47, P_bn_kW=156.17)

sat_model = SaturationModel(bv=1.0, k=1.0, Cm=0.160, m=7)

class Gen160MVA(GeneratorLossModel): 
    def __init__(self): 
        super().__init__(model_data, sat_model)

