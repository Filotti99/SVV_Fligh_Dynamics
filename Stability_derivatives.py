import math
import numpy as np
import inputs
import Aircraft_curves
import cgLocation
import matplotlib.pyplot as plt

def get_CL_alpha(input_matr, ref):
    Alpha_list, CL_list = Aircraft_curves.lift_curve(input_matr, ref)
    CL_alpha = (CL_list[-1]-CL_list[0])/(Alpha_list[-1]-Alpha_list[0])
    return CL_alpha

def get_Cm_delta(trim_matr, general_matr, ref_input):
    CN_list = Aircraft_curves.calc_CL(general_matr, ref_input) # C_N approximated as C_L
    delta_e_list = []
    delta_xcg_list = []
    Cm_delta_list = []
    c_bar = inputs.c_bar
    missing = 0
    for row in trim_matr:
        delta_e_list.append(row[6])
        x_cg_temp = -cgLocation.deltaCg(row[11], 0, ref = ref_input) #should this be a minus #TODO
        delta_xcg_list.append(x_cg_temp)
    for i in range(len(CN_list)):
        if delta_e_list[i] != 0:
            Cm_delta_temp = (-1/delta_e_list[i])*CN_list[i]*(delta_xcg_list[i]/c_bar)
            Cm_delta_list.append(Cm_delta_temp)
        else:
            missing += 1
    return np.sum(Cm_delta_list)/(len(Cm_delta_list)-missing)

def get_Cm_delta_new(delta_matr, ref_input):
    CN_list = Aircraft_curves.calc_CL(delta_matr, ref_input) # C_N approximated as C_L
    c_bar = inputs.c_bar
    delta_x_cg = cgLocation.deltaCg(delta_matr[1][11], delta_matr[0][11], ref=ref_input)
    delta_e = delta_matr[1][6] - delta_matr[0][6]
    Cm_delta = (-1 / delta_e) * np.average(CN_list) * (delta_x_cg / c_bar)
    #print('delta_e =',delta_e,"CN =",np.average(CN_list),"delta_xcg",delta_x_cg)
    return Cm_delta

def get_Cm_alpha(trim_matr, delta_matr, ref_input):
    Cm_delta = get_Cm_delta_new(delta_matr, ref_input)
    alpha, delta = Aircraft_curves.elevator_curve(trim_matr)
    d_delta_d_alpha = (delta[-1]-delta[0])/(alpha[-1]-alpha[0])
    print(d_delta_d_alpha)
    Cm_alpha = -d_delta_d_alpha*Cm_delta
    return Cm_alpha

# Based on reference data
#print("Cm_delta =",get_Cm_delta(inputs.trim_matrix, True)*180/math.pi)
print("Cm_delta =",get_Cm_delta_new(inputs.delta_matrix, True)*180/math.pi)
print("Cm_alpha =",get_Cm_alpha(inputs.trim_matrix, inputs.delta_matrix, True)*180/math.pi)
#print("CL_alpha =",get_CL_alpha(inputs.measurement_matrix, True)*180/math.pi)

# Based on real flight test data
#print("Cm_delta =",get_Cm_delta(inputs.trim_matrix_real, False)*180/math.pi)
print("Cm_delta =",get_Cm_delta_new(inputs.delta_matrix_real, False)*180/math.pi)
print("Cm_alpha =",get_Cm_alpha(inputs.trim_matrix_real, inputs.delta_matrix_real, False)*180/math.pi)
#print("CL_alpha =",get_CL_alpha(inputs.measurement_matrix_real, False)*180/math.pi)