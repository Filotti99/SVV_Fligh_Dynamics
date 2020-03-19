# Flight Dynamics Assignment
from Aircraft_curves import red_elevator_curve
from Stability_derivatives import get_CL_alpha, get_Cm_delta, get_Cm_alpha
import inputs

c_md = get_Cm_delta(inputs.trim_matrix_real, inputs.measurement_matrix_real, False)
red_elevator_curve(inputs.measurement_matrix_real, False, c_md)
