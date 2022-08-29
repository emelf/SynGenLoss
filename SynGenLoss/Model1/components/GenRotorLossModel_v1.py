def get_rotor_loss(If: float, If_nom: float, P_fn: float, P_exn: float, P_br: float) -> float: 
    """Extrapolates the rotor copper, brush, and exciter losses from nominal values using the field current If. \n 
    If: Field current of given operating point. \n 
    If_nom: Nominal field current for no-load where Vt = 1.0 \n 
    P_fn: Nominal field losses [pu] \n 
    P_exn: Nominal exciter losses [pu] \n 
    P_br: Nominal brush losses [pu] \n 
    Output: Rotor power losses [pu]"""
    P_rotor = (P_fn + P_br)*(If/If_nom)**2 + P_exn*If/If_nom
    return P_rotor