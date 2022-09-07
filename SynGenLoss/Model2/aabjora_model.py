from .GenSynMachineLosses import GenModel_SS
from .GenDataClass import Model2DataClass

class AabjoraModel(GenModel_SS): 
    def __init__(self): 
        G1_data = Model2DataClass()
        G1_data.standard_params(103, 0.9, 11, 50, 0.0018, 1e-2, 1.059, 0.676, 0.08, 0.198, 0.5) # This is the Original one
        # G1_data.standard_params(103, 0.9, 11, 50, 0.0018, 1e-2, 1.059, 0.676, 0.08, 0.1659776903705144, 0.4356246146911903)
        G1_data.nominal_losses(1.0, 1.0, 2.028, P_an_kW=186.46, P_sn_kW=89.16, P_fn_kW=173.65, P_brn_kW=2.13, 
                            P_exn_kW=15.88, P_cn_kW=211.92, P_wfn_kW=172.92, P_bn_kW=240.9)
        super().__init__(G1_data) 

