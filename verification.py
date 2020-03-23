import numpy as np
import inspect

from Aircraft_curves import V_e_red, de_red

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
    #meas_mat = np.array([[0,0,0,10000,100, 268.338, 0]])
    inputs = np.array([[np.array([[0,0,0,5000,100, 255.65, 0]]), True, False, True],
                       [np.array([[0,0,0,1000,100, 281.65, 0]]), True, False, True],
                       [np.array([[0,0,0,0,100, 288.15, 0]]), True, False, True],
                       [np.array([[0,0,0,5000,100, 255.65, 0]]), True, False, False],
                       [np.array([[0,0,0,1000,100, 281.65, 0]]), True, False, False],
                       [np.array([[0,0,0,0,100, 288.15, 0]]), True, False, False]])
    predicted_output = np.array([[127.843],[104.836],[100],[99.102],[99.868],[100]])
    for i,item in enumerate(inputs):
        verify(V_e_red, item, predicted_output[i] , err = 1e-3, perc = False)
