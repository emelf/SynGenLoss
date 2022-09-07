import numpy as np
class SaturationModel: 
    def __init__(self, bv, k, Cm, m): 
        self.bv = bv
        self.k = k
        self.Cm = Cm
        self.m = m
        
    def calc_ifd(self, Vt, ia, phi, model_data_class) -> float: 
        delta = np.arctan(ia*(model_data_class.Xq*np.cos(phi)-(model_data_class.Ra*np.sin(phi)))/(Vt+(model_data_class.Ra*ia*np.cos(phi))+model_data_class.Xq*ia*np.sin(phi))) #Power load angle 
        egu = Vt*np.cos(delta) + (model_data_class.Ra*ia*np.cos(delta+phi)) + model_data_class.Xd*ia*np.sin(delta+phi) 
        th = np.arctan(ia*(model_data_class.Xp*np.cos(phi) - model_data_class.Ra*np.sin(phi)) / (Vt + (model_data_class.Ra*ia*np.cos(phi)) + model_data_class.Xp*ia*np.sin(phi)) ) #Potiervinkel 
        ep = Vt*np.cos(th) + (model_data_class.Ra*ia*np.cos(th+phi)) + model_data_class.Xp*ia*np.sin(th+phi) 
        ifu = egu/self.bv #Field current ved luftgapslinje 
        ep_bool = ep > 0.55 #Lager en array av True eller False (1 eller 0) slik at ifs kan lett regnes ut av i neste linje i mange tilfeller samtidig
        ifs = ((ep + self.Cm*ep**self.m)*self.k - ep/self.bv) * ep_bool #Field inkludert saturation
        ifd = (ifu + ifs) 
        return ifd, delta, egu
