import math
import numpy
import inputs

def get_C_N():
    return W/(0.5*rho*V**2*S)

def get_Cm_delta():

    Cm_delta = (-1/delta_e)*C_N*(delta_x_cg/c_bar)
    return Cm_delta

def get_Cm_alpha():
    Cm_delta = get_Cm_delta()
    d_delta_d_alpha = get_elevator_trim_slope()
    Cm_alpha = -d_delta_d_alpha*Cm_delta
    return Cm_alpha

def main():
    Cm_delta = get_Cm_delta()
