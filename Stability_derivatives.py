import math
import numpy as np
import inputs
import Aircraft_curves
import cgLocation

def get_C_N():
    return W/(0.5*rho*V**2*S)

def get_CL_alpha(input_matr):
    Alpha_list, CL_list = Aircraft_curves.lift_curve(input_matr)
    CL_alpha = (CL_list[-1]-CL_list[0])/(Alpha_list[-1]-Alpha_list[0])
    return CL_alpha

def get_Cm_delta(trim_matr, general_matr, ref_input):
    CN_list = Aircraft_curves.calc_CL(general_matr) # C_N approximated as C_L
    delta_e_list = []
    delta_xcg_list = []
    Cm_delta_list = []
    c_bar = inputs.c_bar
    for row in trim_matr:
        delta_e_list.append(row[6])
        x_cg_temp = -cgLocation.deltaCg(row[11], 0, ref = ref_input) #should this be a minus #TODO
        delta_xcg_list.append(x_cg_temp)
    for i in range(len(CN_list)):
        if delta_e_list[i] != 0:
            Cm_delta_temp = (-1/delta_e_list[i])*CN_list[i]*(delta_xcg_list[i]/c_bar)
            Cm_delta_list.append(Cm_delta_temp)
    return np.sum(Cm_delta_list)/(len(Cm_delta_list)-1)

def get_Cm_alpha():
    Cm_delta = get_Cm_delta()
    d_delta_d_alpha = get_elevator_trim_slope()
    Cm_alpha = -d_delta_d_alpha*Cm_delta
    return Cm_alpha

print(get_Cm_delta(inputs.trim_matrix, inputs.measurement_matrix, True))