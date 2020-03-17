import math
import numpy as np
import matplotlib.pyplot as plt
import inputs
import Aircraft_curves

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


thrust_matrix = np.transpose(np.array([inputs.measurement_matrix[:,3], calc_M(inputs.measurement_matrix), calc_deltaT(inputs.measurement_matrix), inputs.measurement_matrix[:,6], inputs.measurement_matrix[:,7] ]))

thrust_input = np.savetxt("matlab.DAT", thrust_matrix, delimiter=" ")