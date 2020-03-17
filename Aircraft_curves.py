import math
import numpy as np
import matplotlib.pyplot as plt
import inputs

def calc_e():
    Clalpha = 2*math.pi*inputs.AR/(2+math.sqrt(4+inputs.AR**2))
    CLalpha = Clalpha*(inputs.AR/(inputs.AR+2))
    e = CLalpha/Clalpha
    return e

def calc_M(measurement_matrix):
    M_array = []
    for row in measurement_matrix:
        M = row[4]/math.sqrt(inputs.gamma*inputs.R*row[9])
        M_array.append(M)
    return M_array

def calc_deltaT(measurement_matrix):
    deltaT_array = []
    for row in measurement_matrix:
        T_ISA = inputs.T_0 + (row[3]*inputs.a_layer)
        T_delta = T_ISA-row[9]
        deltaT_array.append(T_delta)
    return deltaT_array

def calc_CL(measurement_matrix):
    C_L_array = []
    for row in measurement_matrix:
        # nr, time, ET, altitude, IAS, alpha, FFl, FFr, Fused, TAT
        rho = (inputs.p_0*(1+(inputs.a_layer*row[3]/inputs.T_0))**(-inputs.g_0/(inputs.a_layer*inputs.R)))/(inputs.R*row[9]) # change to ISA equation
        W = 60500 # change to varying function
        C_L = W/(0.5*rho*row[4]**2*inputs.S)
        C_L_array.append(C_L)
    return C_L_array

def calc_CD(measurement_matrix):
    C_D_array = []
    C_L_usage = calc_CL(measurement_matrix)
    counter = 0
    for row in measurement_matrix:
        # nr, time, ET, altitude, IAS, alpha, FFl, FFr, Fused, TAT
        C_D = 0.04 + (C_L_usage[counter]**2)/(math.pi*inputs.AR*calc_e())
        C_D_array.append(C_D)
        counter += 1
    return C_D_array

def drag_polar(measurement_matrix):
    C_L_array = calc_CL(measurement_matrix)
    C_D_array = calc_CD(measurement_matrix)
    plt.plot(C_L_array, C_D_array)
    plt.show()

def lift_curve(measurement_matrix):
    Alpha_array = [row[5] for row in measurement_matrix]
    C_L_array = calc_CL(measurement_matrix)
    plt.plot(Alpha_array, C_L_array)
    plt.show()

def drag_curve(measurement_matrix):
    Alpha_array = [row[5] for row in measurement_matrix]
    C_D_array = calc_CD(measurement_matrix)
    plt.plot(Alpha_array, C_D_array)
    plt.show()


print(drag_polar(inputs.measurement_matrix))
print(lift_curve(inputs.measurement_matrix))
print(drag_curve(inputs.measurement_matrix))
print(calc_M(inputs.measurement_matrix))
print(calc_deltaT(inputs.measurement_matrix))