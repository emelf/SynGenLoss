def get_constant_losses(Vt: float, V_nom: float, P_cn: float, P_wfn: float, P_bn: float) -> float: 
    """Extrapolates the constant losses with respect to the terminal voltage. \n 
    Vt: Terminal generator voltage. \n 
    P_cn: Nominal core losses at Vt = 1.0. \n 
    P_wfn: Windage and friction losses. \n 
    P_bn: Nominal bearings losses. \n 
    Output: Constant losses [pu]"""
    P_constant = P_cn*(Vt/V_nom)**2 + P_wfn + P_bn
    return P_constant