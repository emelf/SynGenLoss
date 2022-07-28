
def get_stator_loss(Ia: float, Ia_nom: float, P_an: float, P_sn: float) -> float:
    """Extrapolates the stator armature and stray losses based on the armature current. \n 
    Ia: Armature current [pu] \n Ia_nom: Nominal armature current [pu] \n
    P_an: Nominal armature losses [pu]  \n
    P_sn: Nominal stray losses [pu] \n 
    Output: Stator losses [pu]"""
    P_stator = (P_an + P_sn)*(Ia/Ia_nom)**2
    return P_stator
