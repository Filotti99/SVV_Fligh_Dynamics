import numpy as np
import inspect
from math import pi, sin, cos
import matplotlib.pyplot as plt
#import State_Space as SS
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

def calculateSideslip():
		file = open(path+ "Body Yaw Rate[deg_p_s].csv")
		lines = file.readlines()
		file.close()

		yawRates = np.genfromtxt(lines)[indexStart:indexEnd]
		sideslip = [beta_0]
		print(sideslip)
		for i in range(1,len(yawRates)):
			sideslip.append(yawRates[i]*dt + sideslip[-1])

		return np.array(sideslip),yawRates

def calculateSideslip_test():
		file = open("yaw_rate_test.dat")
		lines = file.readlines()
		file.close()

		yawRates = np.genfromtxt(lines)
		sideslip = [beta_0]
		print(sideslip)
		for i in range(1,len(yawRates)):
			sideslip.append(yawRates[i]*dt + sideslip[-1])

		return np.array(sideslip),yawRates

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
#
#    #kts2mps
#    print("kts2mps")
#    f=SS.kts2mps
#    inputs = [1,5,12.3,-1,0]
#    predicted = [0.514444, 2.57222, 6.327667, -0.514444, 0]
#    for i in range(len(inputs)-1):
#        print(inputs[i],f(inputs[i]),predicted[i],almost_equal_perc(f(inputs[i]),predicted[i],0.1,True))
#    print(inputs[-1],f(inputs[-1]),predicted[-1],almost_equal_abs(f(inputs[-1]),predicted[-1],1e-5,True))
#    print()
#
#    #deg2rad
#    print("deg2rad")
#    f=SS.deg2rad
#    inputs = [180,-180,310,-27.212,360,0]
#    predicted = [pi, -pi, np.radians(310), np.radians(-27.212), 2*pi, 0]
#    for i in range(len(inputs)-1):
#        print(inputs[i],f(inputs[i]),predicted[i],almost_equal_perc(f(inputs[i]),predicted[i],0.1,True))
#    print(inputs[-1],f(inputs[-1]),predicted[-1],almost_equal_abs(f(inputs[-1]),predicted[-1],1e-5,True))
#    print()
#
#    #lbs2Kg
#    print("lbs2Kg")
#    f = SS.lbs2Kg
#    inputs = [1,1/3,107.3,0]
#    predicted = [0.453592,0.151197666666666,48.670461,0]
#    for i in range(len(inputs)-1):
#        print(inputs[i],f(inputs[i]),predicted[i],almost_equal_perc(f(inputs[i]),predicted[i],0.1,True))
#    print(inputs[-1],f(inputs[-1]),predicted[-1],almost_equal_abs(f(inputs[-1]),predicted[-1],1e-5,True))
#    print()
#
#    #celsius2K
#    print("celsius2K")
#    f = SS.celsius2K
#    inputs = [0,-22.321,15,-273.15]
#    predicted = [273.15,-22.321+273.15,288.15,0]
#    for i in range(len(inputs)-1):
#        print(inputs[i],f(inputs[i]),predicted[i],almost_equal_perc(f(inputs[i]),predicted[i],0.1,True))
#    print(inputs[-1],f(inputs[-1]),predicted[-1],almost_equal_abs(f(inputs[-1]),predicted[-1],1e-5,True))
#    print()
#
#    #minutes2Seconds
#    print("minutes2Seconds")
#    f = SS.minutes2Seconds
#    inputs = [1,0.5,22.9,1/3,1/7,0]
#    predicted = [60,30,1374,20,60/7,0]
#    for i in range(len(inputs)-1):
#        print(inputs[i],f(inputs[i]),predicted[i],almost_equal_perc(f(inputs[i]),predicted[i],0.1,True))
#    print(inputs[-1],f(inputs[-1]),predicted[-1],almost_equal_abs(f(inputs[-1]),predicted[-1],1e-5,True))
#    print()

    #calculateSideslip
    print('calculateSideslip')
    tStart = 0
    tEnd = 20
    dt = 0.1
    t = np.arange(tStart, tEnd, dt)
    indexStart = int((1/dt)*tStart)
    indexEnd  = int((1/dt)*tEnd)
    path = r"flight_data\matlab_files\\"
    beta_0 = 0
    
#    sideslip, yawrate = calculateSideslip()
#    plt.figure()
#    plt.plot(t,sideslip,color='blue',label='sideslip')
#    plt.plot(t,yawrate,color='red',label='yawrate')
#    plt.xlabel('Time [s]')
#    plt.title('Sideslip and yawrate')
#    plt.legend()
#    plt.show()
    
    yaw_rate_test = []
    sideslip_true = []
    for i in range(len(t)):
        yaw_rate_test.append(cos(t[i]))
        sideslip_true.append(sin(t[i]))
    np.savetxt("yaw_rate_test.dat",yaw_rate_test)
    
    sideslip, yawrate = calculateSideslip_test()
    plt.figure()
    plt.plot(t,sideslip,color='blue',label='sideslip')
    plt.plot(t,yawrate,color='green',label='yawrate')
    plt.xlabel('Time [s]')
    plt.title('Sideslip and yawrate')
    plt.legend()
    plt.show()
    
    error = []
    for i in range(len(sideslip)):
        error.append(almost_equal_abs(sideslip[i],sideslip_true[i],1,True))
    plt.figure()
    plt.plot(t,error,color='red')
    plt.xlabel('Time [s]')
    plt.ylabel('Error')
    plt.title('Error of sideslip calculations')
    
    print(yawrate[np.argmax(sideslip)])