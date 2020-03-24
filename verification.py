import numpy as np
import inspect
from math import pi
import matplotlib.pyplot as plt

from State_Space import kts2mps, deg2rad, lbs2Kg, celsius2K, minutes2Seconds
from Aircraft_curves import almost_equal_perc, almost_equal_abs, V_e_red, de_red

def verify(f, inputs:np.ndarray, predicted_output, err: float, perc: bool):
    assert inputs.shape[0] == len(inspect.getfullargspec(f)[0])

    outputs = f(*inputs)
    #print(outputs.shape, predicted_output.shape)
    assert outputs.shape == predicted_output.shape
    actual_err = np.abs(outputs-predicted_output)*100/predicted_output if perc else np.average(np.abs(outputs-predicted_output))
    if actual_err < err:
        print("\nTest passed for input: \n", inputs,"\nError = ", actual_err,"\n")
    else:
        print("\nTest failed for input: \n", inputs,"\nError = ", actual_err,"\n")


if __name__ == '__main__':
    plt.close('all')
    #meas_mat = np.array([[0,0,0,10000,100, 268.338, 0]])
    #V_e_red
#    inputs = np.array([[np.array([[0,0,0,5000,100, 255.65, 0]]), True, False, True],
#                       [np.array([[0,0,0,1000,100, 281.65, 0]]), True, False, True],
#                       [np.array([[0,0,0,0,100, 288.15, 0]]), True, False, True],
#                       [np.array([[0,0,0,5000,100, 255.65, 0]]), True, False, False],
#                       [np.array([[0,0,0,1000,100, 281.65, 0]]), True, False, False],
#                       [np.array([[0,0,0,0,100, 288.15, 0]]), True, False, False]])
#    predicted_output = np.array([[127.843],[104.836],[100],[99.102],[99.868],[100]])
#    for i,item in enumerate(inputs):
#        verify(V_e_red, item, predicted_output[i] , err = 1e-3, perc = False)

    #kts2mps
    print("kts2mps")
    inputs = [1,5,12.3,-1,0]
    predicted = [0.514444, 2.57222, 6.327667, -0.514444, 0]
    for i in range(len(inputs)-1):
        print(inputs[i],kts2mps(inputs[i]),predicted[i],almost_equal_perc(kts2mps(inputs[i]),predicted[i],0.1,True))
    print(inputs[-1],kts2mps(inputs[-1]),predicted[-1],almost_equal_abs(kts2mps(inputs[-1]),predicted[-1],1e-5,True))
    print()
    
    #deg2rad
    print("deg2rad")
    inputs = [180,-180,310,-27.212,360,0]
    predicted = [pi, -pi, np.radians(310), np.radians(-27.212), 2*pi, 0]
    for i in range(len(inputs)-1):
        print(inputs[i],deg2rad(inputs[i]),predicted[i],almost_equal_perc(deg2rad(inputs[i]),predicted[i],0.1,True))
    print(inputs[-1],deg2rad(inputs[-1]),predicted[-1],almost_equal_abs(deg2rad(inputs[-1]),predicted[-1],1e-5,True))
    print()
    
    #lbs2Kg
    print("lbs2Kg")
    f = lbs2Kg
    inputs = [1,1/3,107.3,0]
    predicted = [0.453592,0.151197666666666,48.670461,0]
    for i in range(len(inputs)-1):
        print(inputs[i],f(inputs[i]),predicted[i],almost_equal_perc(f(inputs[i]),predicted[i],0.1,True))
    print(inputs[-1],f(inputs[-1]),predicted[-1],almost_equal_abs(f(inputs[-1]),predicted[-1],1e-5,True))
    print()
    
    #celsius2K
    print("celsius2K")
    f = celsius2K
    inputs = [0,-22.321,15,-273.15]
    predicted = [273.15,-22.321+273.15,288.15,0]
    for i in range(len(inputs)-1):
        print(inputs[i],f(inputs[i]),predicted[i],almost_equal_perc(f(inputs[i]),predicted[i],0.1,True))
    print(inputs[-1],f(inputs[-1]),predicted[-1],almost_equal_abs(f(inputs[-1]),predicted[-1],1e-5,True))
    print()
    
    #minutes2Seconds
    print("minutes2Seconds")
    f = minutes2Seconds
    inputs = [1,0.5,22.9,1/3,1/7,0]
    predicted = [60,30,1374,20,60/7,0]
    for i in range(len(inputs)-1):
        print(inputs[i],f(inputs[i]),predicted[i],almost_equal_perc(f(inputs[i]),predicted[i],0.1,True))
    print(inputs[-1],f(inputs[-1]),predicted[-1],almost_equal_abs(f(inputs[-1]),predicted[-1],1e-5,True))
    print()
    