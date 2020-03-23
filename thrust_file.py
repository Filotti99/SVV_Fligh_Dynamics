'''
This file generates the matlab.dat files needed to run thrust.exe and compute the thrust values for different measurement points
'''
import math
import numpy as np
import inputs
import Aircraft_curves as ac

m_dot_fs = [0.048] * 6 #nominal fuel flow for sustained flight, used for Tcs
m_dot_fs_trim = [0.048] * 7 


def calc_M(measurement_matrix):
    '''
    Inputs:
     - measurement_matrix = a matrix of the inputs file
     Outputs:
      - An array with the Mach number at each measurement point (= row) of the matrix
    '''
    M_array = []
    for row in measurement_matrix:
        M = row[4]/math.sqrt(inputs.gamma*inputs.R*row[-2])
        M_array.append(M)
    return M_array

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
        T_delta = T_ISA -row[-2]
        deltaT_array.append(T_delta)
    return deltaT_array


#verification of functions:
"""
calcM Test 1
At ISA 0ft values, a list of M values is given V values
https://www.engineeringtoolbox.com/specific-heat-ratio-d_602.html
"""
x = np.array([[0,0,0,0,100,0,0,0,0,288.15], [0,0,0,0,200,0,0,0,0,288.15]])
y = [100/math.sqrt(inputs.gamma*inputs.R*inputs.T_0), 200/math.sqrt(inputs.gamma*inputs.R*inputs.T_0)]
x2 = np.array([[0,0,0,0,0,0,0,0,0,288.15]])
y2 = [0]
for i in range(len(y)):
    ac.almost_equal_perc(calc_M(x)[i], y[i], 0.1, True)
for i in range(len(x2)):
    ac.almost_equal_abs(calc_M(x2)[i], y2[i], 10**(-2), True)

"""
calcdeltaT Test 1
Difference between T_ISA equation and ISA from online sources should be small
https://www.digitaldutch.com/atmoscalc/
"""
x = np.array([[0,0,0,0,0,0,0,0,0,288.15], [0,0,0,1000,0,0,0,0,0,281.650], [0,0,0,10000,0,0,0,0,0,223.150]])
y = [0, 0, 0]
for i in range(len(x)):
    ac.almost_equal_abs(calc_deltaT(x)[i], y[i], 10**(-2), True)


#assembling matrices:
thrust_matrix = np.transpose(np.array([inputs.measurement_matrix[:,3], calc_M(inputs.measurement_matrix), calc_deltaT(inputs.measurement_matrix), inputs.measurement_matrix[:,6], inputs.measurement_matrix[:,7] ]))
thrust_matrix_nominal = np.transpose(np.array([inputs.measurement_matrix[:,3], calc_M(inputs.measurement_matrix), calc_deltaT(inputs.measurement_matrix), m_dot_fs, m_dot_fs ]))

thrust_matrix_real = np.transpose(np.array([inputs.measurement_matrix_real[:,3], calc_M(inputs.measurement_matrix_real), calc_deltaT(inputs.measurement_matrix_real), inputs.measurement_matrix_real[:,6], inputs.measurement_matrix_real[:,7] ]))
thrust_matrix_real_nominal =np.transpose(np.array([inputs.measurement_matrix_real[:,3], calc_M(inputs.measurement_matrix_real), calc_deltaT(inputs.measurement_matrix_real), m_dot_fs, m_dot_fs ]))

thrust_trim_matrix = np.transpose(np.array([inputs.trim_matrix[:,3], calc_M(inputs.trim_matrix), calc_deltaT(inputs.trim_matrix), inputs.trim_matrix[:,9], inputs.trim_matrix[:,10] ]))
thrust_trim_matrix_nominal = np.transpose(np.array([inputs.trim_matrix[:,3], calc_M(inputs.trim_matrix), calc_deltaT(inputs.trim_matrix), m_dot_fs_trim, m_dot_fs_trim ]))

thrust_trim_matrix_real = np.transpose(np.array([inputs.trim_matrix_real[:,3], calc_M(inputs.trim_matrix_real), calc_deltaT(inputs.trim_matrix_real), inputs.trim_matrix_real[:,9], inputs.trim_matrix_real[:,10] ]))
thrust_trim_matrix_real_nominal = np.transpose(np.array([inputs.trim_matrix_real[:,3], calc_M(inputs.trim_matrix_real), calc_deltaT(inputs.trim_matrix_real), m_dot_fs_trim, m_dot_fs_trim ]))


#for verification purposes
testH = [0]*5
testM = calc_M([[0,0,0,0,100,0,0,0,0,288.15]])*5
testT = calc_deltaT([[0,0,0,0,100,0,0,0,0,288.15]])*5
testmdotl = [0,100,200,300,400,500]
testmdotr = [50,150,250,350,450,550]
thrust_test = np.transpose(np.array([testH,testM,testT,testmdotl,testmdotr]))


thrust_input = np.savetxt("matlab.DAT", thrust_test, fmt="%.6f", delimiter=" ") #writes file with whatever matrix you want
