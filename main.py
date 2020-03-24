# Flight Dynamics Assignment
import numpy as np
from Aircraft_curves import red_elevator_curve, red_force_curve, calc_W
from Stability_derivatives import get_CL_alpha, get_Cm_delta, get_Cm_alpha
from cgLocation import deltaCg
import inputs

plot = False
if plot:
    c_md = get_Cm_delta(inputs.delta_matrix_real, False)
    red_elevator_curve(inputs.trim_matrix_real, False, c_md)
    red_force_curve(inputs.trim_matrix, False)


def deltaCg_mat(meas_mat, l_p1):
    cg_mat = np.array([])
    for meas in meas_mat:
        cg_mat = np.append(cg_mat, deltaCg(meas[-3],inputs.w_fuel_real, a = False, ref = False, l_p0 = l_p, l_p1 = l_p1))
    return cg_mat


l_p = np.array([131,131,214,214,251,251,288,288,170])*0.0254
l_p_delt = np.array([131,131,214,214,251,251,134,288,170])*0.0254
x_cg = 8.515204514568698
weight_meas = calc_W(inputs.w_fuel_real, inputs.measurement_matrix_real, ref=False)
weight_trim = calc_W(inputs.w_fuel_real, inputs.trim_matrix_real, ref=False)
weight_delt = calc_W(inputs.w_fuel_real, inputs.delta_matrix_real, ref=False)
weight = np.hstack((weight_meas,weight_trim,weight_delt))
cg_meas = deltaCg_mat(inputs.measurement_matrix_real, l_p) + x_cg
cg_trim = deltaCg_mat(inputs.trim_matrix_real, l_p) + x_cg
cg_delt = np.hstack((deltaCg(inputs.delta_matrix_real[0,-3],inputs.w_fuel_real, a = False, ref = False, l_p0 = l_p, l_p1 = l_p),
                    deltaCg(inputs.delta_matrix_real[0,-3],inputs.w_fuel_real, a = False, ref = False, l_p0 = l_p, l_p1 = l_p_delt))) + x_cg
cg = np.hstack((cg_meas,cg_trim, cg_delt))
weight_cg = np.vstack((weight,cg)).transpose()
np.savetxt("table.dat", weight_cg, fmt='%.4f',delimiter = "&", newline='\\ \n')
