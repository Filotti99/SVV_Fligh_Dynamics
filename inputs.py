import numpy as np

# Aircraft Parameters

S = 30 # m^2
b = 15.911
c_bar = 2.0569
AR = b*b/S

# Atmosphere Parameters

p_0 = 101325 # Pa
a_layer = -0.0065 # K/m
T_0 = 288.15 # K
g_0 = 9.80665 # m/s^2
R = 287.0
gamma = 1.4

# Static measurement matrix
# nr, time, ET, altitude, IAS, alpha, FFi, FFr, Fused, TAT
measurement_matrix = np.array([[1,30, 2000, 5010, 249, 1.7, 798, 813, 360, 12.5],[2,2137, 2000, 5020, 221, 2.4, 633, 682, 412, 10.5],[3,2436, 2000, 5020, 192, 3.6, 561, 579, 447, 8.8],[4,2604, 2000, 5030, 163, 5.4, 463, 484, 478, 7.2],[5,2947, 2000, 5020, 130, 8.7, 443, 467, 532, 6],[6,3200, 2000, 5110, 118, 10.6, 474, 499, 570, 5.2] ])
for row in measurement_matrix:
    row[3] = row[3] * 0.3048
    row[4] = row[4] * 0.514444
    row[6] = row[6] * (0.453592/3600)
    row[7] = row[7] * (0.453592/3600)
    row[8] = row[8] * 0.453592
    row[9] = row[9] + 273.15