import math
import numpy as np
import matplotlib.pyplot as plt
import inputs


def get_Thrust():
    Thrust_matrix = []
    for row in np.genfromtxt("Thrust_reference.dat"):
        T = sum(row)
        Thrust_matrix.append(T)
    return Thrust_matrix

def calc_e():
    Clalpha = 2*math.pi*inputs.AR/(2+math.sqrt(4+inputs.AR**2))
    CLalpha = Clalpha*(inputs.AR/(inputs.AR+2))
    e = CLalpha/Clalpha
    return e

def calc_CD_curve(measurement_matrix):
    D_array = get_Thrust()
    CL_array = calc_CL(measurement_matrix)
    CD_array = []
    for i in range(len(measurement_matrix)):
        rho = (inputs.p_0*(1+(inputs.a_layer*measurement_matrix[i][3]/inputs.T_0))**(-inputs.g_0/(inputs.a_layer*inputs.R)))/(inputs.R*measurement_matrix[i][9])
        CD_array.append(D_array[i]/(0.5*rho*measurement_matrix[i][4]**2*inputs.S))
    
    e_list = []
    for i in range(len(D_array)-1):
        slope = (CD_array[i+1] -CD_array[i]) / ((CL_array[i+1]**2) -(CL_array[i]**2))
        e_list.append((slope*math.pi*inputs.AR)**-1)
    e = np.average(e_list)
    
    CD0_list = []
    for i in range(len(CD_array)):
        CD0_list.append(CD_array[i] -(CL_array[i]**2/(math.pi*inputs.AR*e)))
    CD0 = np.average(CD0_list)
    
    return e,e_list,CD0,CD0_list,CD_array

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
        # nr, time, ET, altitude, IAS, alpha, FFl, FFr, Fused, TAT, W
        rho = (inputs.p_0*(1+(inputs.a_layer*row[3]/inputs.T_0))**(-inputs.g_0/(inputs.a_layer*inputs.R)))/(inputs.R*row[9]) # change to ISA equation
        C_L = row[10]/(0.5*rho*row[4]**2*inputs.S)
        C_L_array.append(C_L)
    return C_L_array

def calc_CD(measurement_matrix):
    C_D_array = []
    C_L_usage = calc_CL(measurement_matrix)
    counter = 0
    for row in measurement_matrix:
        # nr, time, ET, altitude, IAS, alpha, FFl, FFr, Fused, TAT, W
        C_D = 0.04 + (C_L_usage[counter]**2)/(math.pi*inputs.AR*calc_e())
        C_D_array.append(C_D)
        counter += 1
    return C_D_array

def drag_polar(measurement_matrix):
    C_L_array = calc_CL(measurement_matrix)
    C_D_array = calc_CD(measurement_matrix)
    plt.plot(C_L_array, C_D_array)
    plt.show()
    return C_L_array, C_D_array

def lift_curve(measurement_matrix):
    Alpha_array = [row[5] for row in measurement_matrix]
    C_L_array = calc_CL(measurement_matrix)
    plt.plot(Alpha_array, C_L_array)
    plt.show()
    return Alpha_array, C_L_array

def drag_curve(measurement_matrix):
    Alpha_array = [row[5] for row in measurement_matrix]
    C_D_array = calc_CD(measurement_matrix)
    plt.plot(Alpha_array, C_D_array)
    plt.show()
    return Alpha_array, C_D_array


#print(drag_polar(inputs.measurement_matrix_real))
#print(lift_curve(inputs.measurement_matrix_real))
#print(drag_curve(inputs.measurement_matrix_real))
#print(calc_CL(inputs.measurement_matrix_real))
#print(calc_M(inputs.measurement_matrix_real))
#print(calc_deltaT(inputs.measurement_matrix_real))
