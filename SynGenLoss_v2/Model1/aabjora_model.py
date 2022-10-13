import numpy as np
from .components.GenSaturationModel_v1 import SaturationModel
from .DataClasses import GenDataClass, TrafoDataClass
from .ModelClasses import GeneratorModel, TrafoModel, PowerPlantModel

model_data = GenDataClass() 
model_data.standard_params(Sn_mva=103, V_nom_kV=11.0, cos_phi=0.9, If_nom_A=525.15, Ra=0.00182, Xd=1.059, Xq=0.676, Xp=0.141)
model_data.nominal_losses(V_nom=1.0, Ia_nom=5406.1, If_nom=1065, P_an_kW=186.46, P_sn_kW=89.16, P_fn_kW=173.65, P_brn_kW=2.13, P_exn_kW=15.88, 
                          P_cn_kW=211.92, P_wfn_kW=172.92, P_bn_kW=240.90)

sat_model = SaturationModel(bv=1.0, k=1.0308, Cm=0.160, m=7)

AabjoraGen = GeneratorModel(model_data, sat_model)

T1_data = TrafoDataClass()
T1_data.define_params(103, 132, 11, 0.129, 0.0031, 0.001)
AabjoraTrafo = TrafoModel(T1_data)

AabjoraPlant = PowerPlantModel([AabjoraGen], [AabjoraTrafo])