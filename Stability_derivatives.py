import math
import numpy
import inputs
import Aircraft_curves

def get_C_N():
    return W/(0.5*rho*V**2*S)

def get_CL_alpha(input_matr):
    Alpha_list, CL_list = Aircraft_curves.lift_curve(input_matr)
    CL_alpha = (CL_list[-1]-CL_list[0])/(Alpha_list[-1]-Alpha_list[0])
    return CL_alpha

def get_Cm_delta(input_matr):
    CN_list = Aircraft_curves.calc_CL(input_matr) # C_N approximated as C_L
    delta_e_list =
    delta_xcg_list =
    c_bar =
    Cm_delta = (-1/delta_e)*C_N*(delta_x_cg/c_bar)
    return Cm_delta

def get_Cm_alpha():
    Cm_delta = get_Cm_delta()
    d_delta_d_alpha = get_elevator_trim_slope()
    Cm_alpha = -d_delta_d_alpha*Cm_delta
    return Cm_alpha
