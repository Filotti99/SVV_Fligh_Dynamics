import math
import numpy as np
import matplotlib.pyplot as plt
import inputs
from tools import interpolate
from scipy import stats

def almost_equal_perc(a, b, percentage=1, return_value=False):
    """
    Tests if two floating point numbers are almost equal by comparing their percentual difference against a
    set value. Default significance is 1%.

    :param a, b: Floating point numbers to compare
    :param percentage = Maximum allowed difference in percentage
    :raises: error when difference falls outside bounds
    """
    assert 100*abs(a-b)/(0.5*(a+b)) <= percentage
    if return_value:
        return 100*abs(a-b)/(0.5*(a+b))

def almost_equal_abs(a, b, difference=0.01, return_value=False):
    """
    Tests if two floating point numbers are almost equal by comparing their absolute difference

    :param a, b: Floating point numbers to compare
    :param difference = Maximum allowed absolute difference
    :raises: error when difference falls outside bounds
    """
    assert abs(a-b) <= difference
    if return_value:
        return abs(a-b)


def get_Thrust(reality:bool, nominal:bool, trim:bool):
    fname = "Thrust"
    if trim:
        fname += "_trim"
    if reality:
        fname += "_real"
    if nominal:
        fname += "_nominal"
    fname += ".dat"

    Thrust_matrix = []
    for row in np.genfromtxt(fname):
        Thrust_matrix.append(sum(row))
    return Thrust_matrix

def calc_Tc(measurement_matrix, reality:bool, nominal:bool, trim:bool):
    Thrust_matrix = get_Thrust(reality,nominal,trim)
    Tc_array = []
    for i in range(len(measurement_matrix)):
        Tc_array.append(Thrust_matrix[i] / (0.5*inputs.rho_0*measurement_matrix[i][4]**2*inputs.d**2))
    return Tc_array


def calc_W(w_f0: float,meas_mat: np.ndarray, ref = True) -> np.ndarray:

    path = "cg_data/pass_w_ref.dat" if ref else "cg_data/pass_w.dat"

    w_pass = np.genfromtxt(path)*inputs.g_0
    w_f = w_f0 - meas_mat[:,-2]

    return np.sum(w_pass)+ w_f + inputs.w_oew


def calc_M(measurement_matrix):
    '''
    Inputs:
     - measurement_matrix = a matrix of the inputs file

     Outputs:
      - An array with the Mach number at each measurement point (= row) of the matrix
    '''
    M_array = []
    for row in measurement_matrix:
        M = row[4]/math.sqrt(inputs.gamma*inputs.R*row[9])
        M_array.append(M)
    return M_array


def V_e_red(meas_matrix: np.ndarray, ref: bool, tilda = True, vtas = False):
    p   = inputs.p_0*(1+inputs.a_layer*meas_matrix[:,3]/inputs.T_0)**(-inputs.g_0/(inputs.R*inputs.a_layer))
    M   = np.sqrt((2/(inputs.gamma-1))*((1+inputs.p_0/p*((1+(inputs.gamma-1)/(2*inputs.gamma)*inputs.rho_0/inputs.p_0*meas_matrix[:,4]**2)**(inputs.gamma/(inputs.gamma-1))-1))**((inputs.gamma-1)/inputs.gamma)-1))
    T   = meas_matrix[:,-2]/(1+(inputs.gamma-1)/2*M**2)
    if vtas:
        print(M*np.sqrt(inputs.gamma*inputs.R*T))
        return M*np.sqrt(inputs.gamma*inputs.R*T)
    V   = M*np.sqrt(inputs.gamma*p/inputs.rho_0)

    w_f0 = 4050 if ref else 2640
    w_f0 *= inputs.lbs*inputs.g_0

    return V*np.sqrt(inputs.W_s/meas_matrix[:,-1]) if tilda else V


def de_red(meas_mat: np.ndarray, c_md: float, Tcs: np.ndarray, Tc: np.ndarray):
    if meas_mat.shape[1] < 13:
        return 0

    c_mtc = -0.0064

    return meas_mat[:, 6] - (c_mtc/c_md)*(Tcs-Tc)


def calc_deltaT(measurement_matrix):
    '''
    Inputs:
     - measurement_matrix = a matrix of the inputs file

     Outputs:
      - An array with the temperature differential with ISA at each measurement point (= row) of the matrix
    '''
    deltaT_array = []
    for row in measurement_matrix:
        T_ISA = inputs.T_0 + (row[3]*inputs.a_layer)
        T_delta = T_ISA-row[9]
        deltaT_array.append(T_delta)
    return deltaT_array


def calc_CL(measurement_matrix, ref):
    '''
    Inputs:
     - measurement_matrix = a matrix of the inputs file
     - ref = a bool that is true if the reference data is used and false if the flight test data is used

     Outputs:
      - An array with the C_L at each measurement point (= row) of the matrix

      Note:
      - An option to use rho at altitude and V_t is left in the code, but is put in comments for now
    '''
    C_L_array = []
    C_L_other = []
    V_e_array = V_e_red(measurement_matrix, ref, False, False)  # array with the equivalent airspeed
    V_t_array = V_e_red(measurement_matrix, ref, False, True)  # array with the true airspeed
    counter = 0
    for row in measurement_matrix:
        # nr, time, ET, altitude, IAS, alpha, FFl, FFr, Fused, TAT, W
        rho = inputs.rho_0
        C_L = row[-1]/(0.5*rho*V_e_array[counter]**2*inputs.S)
        C_L_array.append(C_L)
       # rho = (inputs.p_0*(1+(inputs.a_layer*row[3]/inputs.T_0))**(-inputs.g_0/(inputs.a_layer*inputs.R)))/(inputs.R*row[9]) # change to ISA equation
       # C_L = row[10] / (0.5 * rho * V_t_array[counter] ** 2 * inputs.S)
       # C_L_other.append(C_L)
        counter += 1
    return C_L_array#, C_L_other

def calc_CD_curve(measurement_matrix, reality, ref):
    '''
    Inputs:
     - measurement_matrix = a matrix of the inputs file
     - reality = a bool that is true if the flight test data is used and false if the reference data is used
     - ref = a bool that is true if the reference data is used and false if the flight test data is used

     Outputs:
      - The estimated value of 'e'
      - The estimated value of C_D0
      - An array with the C_D at each measurement point (= row) of the matrix
    '''
    D_array = get_Thrust(reality, False, False)
    CL_array = calc_CL(measurement_matrix, not reality)
    V_e_array = V_e_red(measurement_matrix, not reality, False, False)
    CD_array = []
    CL2_array = []
    for i in range(len(measurement_matrix)):
        rho = inputs.rho_0
        CD_array.append(D_array[i]/(0.5*rho*V_e_array[i]**2*inputs.S))
        CL2_array.append(CL_array[i]**2)

    slope, CD0, r_value, p_value, std_err = stats.linregress(CL2_array,CD_array)
    e = (slope * math.pi * inputs.AR)**-1

    return e,CD0,CD_array

def drag_polar(measurement_matrix,reality):
    '''
    Inputs:
     - measurement_matrix = a matrix of the inputs file
     - reality = a bool that is true if the flight test data is used and false if the reference data is used

     Outputs:
      - A plot of the C_L values on the x-axis an the C_D values on the y-axis (not returned)
      - Returns an array with the C_L at each measurement point (= row) of the matrix
      - Returns an array with the C_D at each measurement point (= row) of the matrix
    '''
    C_L_array = calc_CL(measurement_matrix, not reality)
    e, CD0, C_D_array = calc_CD_curve(measurement_matrix, reality, not reality)
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

def lift_curve(measurement_matrix, ref):
    '''
    Inputs:
     - measurement_matrix = a matrix of the inputs file
     - ref = a bool that is true if the reference data is used and false if the flight test data is used

     Outputs:
      - A plot of the alpha values on the x-axis an the C_L values on the y-axis (not returned)
      - Returns an array with the alpha at each measurement point (= row) of the matrix
      - Returns an array with the C_L at each measurement point (= row) of the matrix
    '''
    Alpha_array = [row[5] for row in measurement_matrix]
    C_L_array = calc_CL(measurement_matrix, ref)
    plt.plot(Alpha_array, C_L_array)
    plt.title('Lift coefficient curve as a function of the angle of attack')
    plt.xlabel('Angle of attack [degrees]')
    plt.ylabel('Lift coefficient [-]')
    plt.show()
    return Alpha_array, C_L_array


def drag_curve(measurement_matrix, reality):
    '''
    Inputs:
     - measurement_matrix = a matrix of the inputs file
     - reality = a bool that is true if the flight test data is used and false if the reference data is used

     Outputs:
      - A plot of the alpha values on the x-axis an the C_D values on the y-axis (not returned)
      - Returns an array with the alpha at each measurement point (= row) of the matrix
      - Returns an array with the C_D at each measurement point (= row) of the matrix
    '''
    Alpha_array = [row[5] for row in measurement_matrix]
    e, CD0, C_D_array = calc_CD_curve(measurement_matrix, reality, not reality)
    plt.plot(Alpha_array, C_D_array)
    plt.title('Drag coefficient curve as a function of the angle of attack')
    plt.xlabel('Angle of attack [deg]')
    plt.ylabel('Drag coefficient [-]')
    plt.show()
    return Alpha_array, C_D_array

def elevator_curve_alpha(measurement_matrix):
    '''
    Inputs:
     - measurement_matrix = a matrix of the inputs file

     Outputs:
      - A plot of the alpha values on the x-axis an the elevator deflection values on the y-axis (not returned)
      - Returns an array with the alpha at each measurement point (= row) of the matrix
      - Returns an array with the elevator deflection at each measurement point (= row) of the matrix
    '''
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
    plt.xlabel("Angle of attack [degrees]")
    plt.ylabel("Elevator deflection [rad]")
    plt.title('Elevator trim curve')
    plt.show()
    return Alpha_array, De_array

def red_elevator_curve(trim_mat:np.ndarray, ref: bool, c_md: float):
    Tcs = np.array(calc_Tc(trim_mat, not ref, True))
    Tc  = np.array(calc_Tc(trim_mat, not ref, False))

    V_e_tilda = V_e_red(trim_mat, ref, tilda=True)
    d_e_star  = de_red(trim_mat, c_md, Tcs, Tc)

    plt.figure("Reduced Elevator Deflection Curve")
    plt.plot(V_e_tilda, d_e_star)
    plt.xlabel("Reduced Airspeed [m/s]")
    plt.ylabel("Reduced Elevator Deflection [$^{\cirlce}$]")
    plt.title("Reduced Elevator Deflection Curve")
    plt.savefig("figures/red_el_curve.png")
    plt.show()

# Test functions (commented to prevent plots from being spammed when running the file)

# Reference Data
# print(calc_M(inputs.measurement_matrix))
# print(calc_deltaT(inputs.measurement_matrix))
# print(calc_CL(inputs.measurement_matrix, ref = True))
# print(calc_CD_curve(inputs.measurement_matrix, reality = False, ref = True))
# print(drag_polar(inputs.measurement_matrix, reality = False))
# print(lift_curve(inputs.measurement_matrix, ref = True))
# print(drag_curve(inputs.measurement_matrix, reality = False))
# print(elevator_curve_alpha(inputs.measurement_matrix))

# Flight Test Data
# print(calc_M(inputs.measurement_matrix_real))
# print(calc_deltaT(inputs.measurement_matrix_real))
# print(calc_CL(inputs.measurement_matrix_real, ref = False))
# print(calc_CD_curve(inputs.measurement_matrix_real, reality = True, ref = False))
# print(drag_polar(inputs.measurement_matrix_real, reality = True))
# print(lift_curve(inputs.measurement_matrix_real, ref = False))
# print(drag_curve(inputs.measurement_matrix_real, reality = True))
# print(elevator_curve_alpha(inputs.measurement_matrix_real))

# Unit tests performed if you run this file


if __name__ == '__main__':
    """
    calcM Test 1
    At ISA 0ft values, a list of M values is given V values
    https://www.engineeringtoolbox.com/specific-heat-ratio-d_602.html
    """
    x = np.array([[0,0,0,0,100,0,0,0,0,288.15], [0,0,0,0,200,0,0,0,0,288.15]])
    y = [100/math.sqrt(1.401*287.057*288.15), 200/math.sqrt(1.401*287.057*288.15)]
    x2 = np.array([[0,0,0,0,0,0,0,0,0,288.15]])
    y2 = [0]
    for i in range(len(y)):
        almost_equal_perc(calc_M(x)[i], y[i], 0.1)
    for i in range(len(x2)):
        almost_equal_abs(calc_M(x2)[i], y2[i], 10**(-2))

    """
    calcdeltaT Test 1
    Difference between T_ISA equation and ISA from online sources should be small
    https://www.engineeringtoolbox.com/specific-heat-ratio-d_602.html
    """
    x = np.array([[0,0,0,0,0,0,0,0,0,288.15], [0,0,0,1000,0,0,0,0,0,281.650], [0,0,0,10000,0,0,0,0,0,223.150]])
    y = [0, 0, 0]
    for i in range(len(x)):
        almost_equal_abs(calc_deltaT(x)[i], y[i], 10**(-2))