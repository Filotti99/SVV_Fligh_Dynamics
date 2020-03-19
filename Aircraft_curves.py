import math
import numpy as np
import matplotlib.pyplot as plt
import inputs
from tools import interpolate


def get_Thrust(reality): #reality must be boolean True if real data, False if reference
    if reality:
        fname = "Thrust_reference.dat"
    else:
        fname = "Thrust_real.dat"
    Thrust_matrix = []
    for row in np.genfromtxt(fname):
        Thrust_matrix.append(sum(row))
    return Thrust_matrix

def calc_W(w_f0: float,meas_mat: np.ndarray, ref = True) -> np.ndarray:

    path = "cg_data/pass_w_ref.dat" if ref else "cg_data/pass_w.dat"

    w_pass = np.genfromtxt(path)*inputs.g_0
    w_f = w_f0 - meas_mat[:,-2]

    return np.sum(w_pass)+ w_f + inputs.w_oew

#def calc_e(): #old version, use calc_CD_curve
#    Clalpha = 2*math.pi*inputs.AR/(2+math.sqrt(4+inputs.AR**2))
#    CLalpha = Clalpha*(inputs.AR/(inputs.AR+2))
#    e = CLalpha/Clalpha
#    return e

def calc_M(measurement_matrix):
    M_array = []
    for row in measurement_matrix:
        M = row[4]/math.sqrt(inputs.gamma*inputs.R*row[9])
        M_array.append(M)
    return M_array

def V_e_red(meas_matrix: np.ndarray, ref: bool, tilda = True,):
    p   = inputs.p_0*(1+inputs.a_layer*meas_matrix[:,3]/inputs.T_0)**(-inputs.g_0/(inputs.R*inputs.gamma))
    M   = np.sqrt((2/(inputs.gamma-1))*((1+inputs.p_0/p*((1+(inputs.gamma-1)/(2*inputs.gamma)*inputs.rho_0/inputs.p_0*meas_matrix[:,4]**2)**(inputs.gamma/(inputs.gamma-1))-1))**((inputs.gamma-1)/inputs.gamma)))
    T   = meas_matrix[:,-2]/(1+(inputs.gamma-1)/2*M**2)
    V   = M*np.sqrt(inputs.gamma*p/inputs.rho_0)

    w_f0 = 4050 if ref else 2640
    w_f0 *= inputs.lbs*inputs.g_0

    return V*np.sqrt(inputs.W_s/meas_matrix[:,-1]) if tilda else V

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
        #rho = (inputs.p_0*(1+(inputs.a_layer*row[3]/inputs.T_0))**(-inputs.g_0/(inputs.a_layer*inputs.R)))/(inputs.R*row[9]) # change to ISA equation
        rho = inputs.rho_0
        C_L = row[10]/(0.5*rho*row[4]**2*inputs.S)
        C_L_array.append(C_L)
    return C_L_array

#def calc_CD(measurement_matrix): #Old method, use calc_CD_curve
#    C_D_array = []
#    C_L_usage = calc_CL(measurement_matrix)
#    counter = 0
#    for row in measurement_matrix:
#        # nr, time, ET, altitude, IAS, alpha, FFl, FFr, Fused, TAT, W
#        C_D = 0.04 + (C_L_usage[counter]**2)/(math.pi*inputs.AR*calc_e())
#        C_D_array.append(C_D)
#        counter += 1
#    return C_D_array


def calc_CD_curve(measurement_matrix,reality):
    D_array = get_Thrust(reality)
    CL_array = calc_CL(measurement_matrix)
    CD_array = []
    for i in range(len(measurement_matrix)):
        #rho = (inputs.p_0*(1+(inputs.a_layer*measurement_matrix[i][3]/inputs.T_0))**(-inputs.g_0/(inputs.a_layer*inputs.R)))/(inputs.R*measurement_matrix[i][9])
        rho = inputs.rho_0
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

    return e,CD0,CD_array

def drag_polar(measurement_matrix,reality):
    C_L_array = calc_CL(measurement_matrix)
    e, CD0, C_D_array = calc_CD_curve(measurement_matrix,reality)
    e = 0.8
    CD0 = 0.04
    C_D_calculated = []
    error=[]
    for i in range(len(C_L_array)):
        C_D_calculated.append(CD0 + C_L_array[i]**2 / (math.pi*inputs.AR*e))
        error.append((abs(C_D_calculated[i]-C_D_array[i])/C_D_calculated[i])*100)
    plt.plot(C_L_array, C_D_array, label='measured')
    plt.plot(C_L_array, C_D_calculated, label='calculated')
    plt.title('CL-CD polar')
    plt.xlabel('CL')
    plt.ylabel('CD')
    plt.legend()
    plt.show()

#    fig, ax1 = plt.subplots()
#    ax1.set_xlabel('CL')
#    ax1.set_ylabel('CD')
#    ax1.plot(C_L_array, C_D_array, color='blue')
#    ax1.plot(C_L_array, C_D_calculated, color='green')
#    ax1.tick_params(axis='y', labelcolor='blue')
#
#    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
#    ax2.set_ylabel('% error', color='red')  # we already handled the x-label with ax1
#    ax2.plot(C_L_array, error, color='red')
#    ax2.tick_params(axis='y', labelcolor='red')
#
#    fig.tight_layout()  # otherwise the right y-label is slightly clipped
#    plt.show()

    return C_L_array, C_D_array

def lift_curve(measurement_matrix):
    Alpha_array = [row[5] for row in measurement_matrix]
    C_L_array = calc_CL(measurement_matrix)
    plt.plot(Alpha_array, C_L_array)
    plt.title('Lift coefficient curve as a function of the angle of attack')
    plt.xlabel('Angle of attack [degrees]')
    plt.ylabel('Lift coefficient [-]')
    plt.show()
    return Alpha_array, C_L_array

def drag_curve(measurement_matrix,reality):
    Alpha_array = [row[5] for row in measurement_matrix]
    e, CD0, C_D_array = calc_CD_curve(measurement_matrix,reality)
    plt.plot(Alpha_array, C_D_array)
    plt.title('Lift coefficient curve as a function of the angle of attack')
    plt.xlabel('Angle of attack [deg]')
    plt.ylabel('Drag coefficient [-]')
    plt.show()
    return Alpha_array, C_D_array

def elevator_curve(measurement_matrix):
    Combined_array = [[row[5],row[6]] for row in measurement_matrix]
    ordered_array = []
    for i in range(len(Combined_array)):
        counter = 0
        first = 100
        best = 0
        for item in Combined_array:
            if item[0] < first:
                first = item[0]
                best = counter
            counter += 1
        ordered_array.append(Combined_array[best])
        del Combined_array[best]
    Alpha_array = [row[0] for row in ordered_array]
    De_array = [row[1] for row in ordered_array]
    plt.plot(Alpha_array, De_array)
    plt.title('Elevator trim curve')
    plt.show()
    return Alpha_array, De_array

#elevator_curve(inputs.trim_matrix)
#print(drag_polar(inputs.measurement_matrix_real))
print(lift_curve(inputs.measurement_matrix_real))
#print(drag_curve(inputs.measurement_matrix_real))
#print(calc_CL(inputs.measurement_matrix))
#print(calc_M(inputs.measurement_matrix_real))
#print(calc_deltaT(inputs.measurement_matrix_real))
