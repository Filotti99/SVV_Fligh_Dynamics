import math
import numpy as np
import matplotlib.pyplot as plt
import inputs
import Aircraft_curves

#with open("matlab.DAT", "w") as file:
print(inputs.measurement_matrix[:,3])
np.savetxt("matlab.DAT", inputs.measurement_matrix[:,3], delimiter=" ")