import numpy as np 
from numpy import cos, sin, pi

# Park transformation: 
T = lambda theta: 2/3 * np.array([[cos(theta), cos(theta - 2*pi/3), cos(theta + 2*pi/3)], 
                                  [-sin(theta), -sin(theta - 2*pi/3), -sin(theta+2*pi/3)], 
                                  [0.5, 0.5, 0.5]]) 
T_inv = lambda theta: np.array([[cos(theta), -sin(theta), 1],
                                [cos(theta - 2*pi/3), -sin(theta - 2*pi/3), 1], 
                                [cos(theta + 2*pi/3), -sin(theta + 2*pi/3), 1]]) 



def calc_L_dq0(Laa0, Laa2, Lab0): 
    """For salient pole machines. \n 
    Returns (Ld, Lq, L0"""
    Ld = Laa0 + Lab0 + 3/2*Laa2 
    Lq = Laa0 + Lab0 - 3/2*Laa2 
    L0 = Laa0 - 2*Lab0 
    return Ld, Lq, L0