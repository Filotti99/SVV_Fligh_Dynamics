import numpy as np
import control.matlab as ml
import control as c
import matplotlib.pyplot as plt

print("start")

#Standard atmosphere
g = 9.81


####INPUTS, Change before use
V = 250
u_0= V
mass = 18200
Weight = mass*g
density  = 1.225

#Using reference values
oswald_factor = 0.8
CD_0 = 0.04
CL_alpha = 5.084
alpha_0 = 0
theta_0 = 0

span = 15.911
cord = 2.0569
S = 30.00
Aspect_ratio = span**2/S

D_b = span/V
D_c = cord/V

mu_b = mass / (density*S*span)
print(mu_b)
mu_c = mass / (density*S*cord)
Kx = np.sqrt(0.019)
Ky = np.sqrt(1.25*1.114)
Kz = np.sqrt(0.042)
Kxz = 0.002

CL = 2*Weight/(density*V**2*S)
print(CL)
CD = CD_0 + ((CL_alpha*alpha_0)**2/(np.pi*Aspect_ratio*oswald_factor))


Cx_0 = Weight * np.sin(theta_0)/(0.5*density*V**2*S)
Cx_q = -0.28170
Cx_u = -0.095
Cx_alpha = 0.47966
Cx_alpha_dot = 0.08330
Cx_delta_e = -0.03728

Cy_beta = -0.75000
Cy_beta_dot = 0
Cy_p = -0.0304
Cy_r = 0.8495
Cy_delta_a = -0.0400
Cy_delta_r = 0.23

Cz_alpha = -5.74340
Cz_alpha_dot = -0.0035
Cz_0 = -Weight*np.cos(theta_0)/(0.5*density*V**2*S)
Cz_u = -0.37616
Cz_q = -5.66290
Cz_delta_e = -0.69612

Cl_beta = -0.10260
Cl_p = -0.71085
Cl_r = 0.23760
Cl_delta_a = -0.23088
Cl_delta_r = 0.03440

Cn_beta = 0.1348
Cn_beta_dot = 0
Cn_p = -0.0602
Cn_r = -0.2061
Cn_delta_a = -0.0120
Cn_delta_r = -0.0939


#Replace values with calculated ones
Cm_q = -8.79415
Cm_u = 0.06990
Cm_alpha = -0.5626
Cm_alpha_dot = 0.1780
Cm_delta_e = -1.1642


#C1sym = np.matrix([[(-2*mu_c*cord)/(V**2),					0,								0,											0],
#				   [0,										(Cz_alpha_dot-2*mu_c)*cord/V,		0,											0],
#				   [0,										0,								-cord/V,									0],
#				   [0,										Cm_alpha_dot*cord/V,				0,										-2*mu_c*Ky**2*cord**2/(V**2)]])

#C2sym = np.matrix([[Cx_u/V,								Cx_alpha,						Cz_0,										Cx_q*cord/V],
#				   [Cz_u/V,									Cz_alpha,						-Cx_0,										(Cz_q+2*mu_c)*cord/V],
#				   [0,										0,								0,											cord/V],
#				   [Cm_u/V,									Cm_alpha,						0,											Cm_q*cord/V]])

#C3sym = np.matrix([[Cx_delta_e],
#				   [Cz_delta_e],
#				   [0],
#				   [Cm_delta_e]])


#C1asym = np.matrix([[(Cy_beta_dot-2*mu_b)*span/V,			0,								0,											0],
#					[0,										-span/(2*V),					0,											0],
#					[0,										0,								(-4*mu_b*Kx**2*span**2)/(2*V**2),			(4*mu_b*Kxz*span**2)/(2*V**2)],
#					[(Cn_beta_dot*span)/(2*V),				0,								(4*mu_b*Kxz*span**2)/(2*V**2),				(-4*mu_b*Kz**2*span**2)/(2*V**2)]])

#C2asym = np.matrix([[Cy_beta,								CL,								Cy_p*span/(2*V),							(Cy_r - 4*mu_b)*span/(2*V)],
#					[0,										0,								span/(2*V),									0],
#					[Cl_beta,								0,								Cl_p*span/(2*V),							Cl_r*span/(2*V)],
#					[Cn_beta,								0,								Cn_p*span/(2*V),							Cn_r*span/(2*V)]])

#C3asym = np.matrix([[Cy_delta_a,							Cy_delta_r],
#					[0,										0],
#					[Cl_delta_a,							Cl_delta_r],
#					[Cn_delta_a,							Cn_delta_r]])

C1sym = D_c*np.matrix([[-2*mu_c			,						0,								0,											0],
						[0,										(Cz_alpha_dot-2*mu_c),			0,											0],
						[0,										0,								-1		,									0],
						[0,										Cm_alpha_dot		,			0,											-2*mu_c*Ky**2]])

C2sym = np.matrix([[Cx_u,									Cx_alpha,						Cz_0,										Cx_q],
				   [Cz_u,									Cz_alpha,						-Cx_0,										Cz_q+2*mu_c],
				   [0,										0,								0,											1],
				   [Cm_u,									Cm_alpha,						0,											Cm_q]])

C3sym = np.matrix([[Cx_delta_e],
				   [Cz_delta_e],
				   [0],
				   [Cm_delta_e]])

C1asym = D_b* np.matrix([[(Cy_beta_dot-2*mu_b),					0,								0,											0],
						[0,										-0.5,							0,											0],
						[0,										0,								-4*mu_b*Kx**2,								4*mu_b*Kxz],
						[Cn_beta_dot,							0,								4*mu_b*Kxz,									-4*mu_b*Kz**2]])

C2asym = np.matrix([[Cy_beta,								CL,								Cy_p,										Cy_r - 4*mu_b],
					[0,										0,								1,											0],
					[Cl_beta,								0,								Cl_p,										Cl_r],
					[Cn_beta,								0,								Cn_p,										Cn_r]])

C3asym = np.matrix([[Cy_delta_a,							Cy_delta_r],
					[0,										0],
					[Cl_delta_a,							Cl_delta_r],
					[Cn_delta_a,							Cn_delta_r]])

Asym = -np.linalg.inv(C1sym)*C2sym
Bsym = -np.linalg.inv(C1sym)*C3sym


Aasym = -np.linalg.inv(C1asym)*C2asym
#Aasym[0] = [0,0,0,2/D_b]




print("symmetric: ")
print("A: ",Asym)
print("B: ",Bsym)
print(V/cord)

print("Asymmetric: ")
print("A: ",Aasym)
Basym = -np.linalg.inv(C1asym)*C3asym
print("B: ",Basym)

print(2*V/span)

print("debug: ")
print("cy_beta_dot = ", Cy_beta_dot)
print("mu_b = ", mu_b)
print("Cn_bet_dot = ", Cn_beta_dot)
print("Kx = ", Kx)
print("Kz = ", Kz)
print("Kxz = ", Kxz)
print("Db = ", D_b)

print("C1-1 = ", np.linalg.inv(C1asym))

print("C2 = ", C2asym)


Csym = np.matrix(np.identity(4))
Casym = np.matrix(np.identity(4))


Dsym = np.matrix([[0],
				  [0],
				  [0],
				  [0]])
Dasym = np.matrix([[0,0],
				   [0,0],
				   [0,0],
				   [0,0]])

systemSym = ml.ss(Asym,Bsym,Csym,Dsym)
systemAsym = ml.ss(Aasym,Basym,Casym,Dasym)

print("Eigenvalues symmetric case: ",np.linalg.eigvals(systemSym.A))
print("Eigenvalues asymmetric case: ", np.linalg.eigvals(systemAsym.A))
Tin = np.arange(0,1000,0.01)
Uin = np.zeros_like(Tin)
Uin[0:100] = np.radians(5)

TSym, ySym, xOut = c.forced_response(systemSym,T = Tin,U=Uin,X0=[0,alpha_0,theta_0,0])
#print(ySym)
#ySym, TSym = ml.step(systemSym,X0=[0,alpha_0,theta_0,0],T = Tin)

ySym[0] = ySym[0]*u_0 + u_0
ySym[3] = ySym[3]*(u_0/cord)

#yAsym, TAsym = ml.impulse(systemAsym,X0=[0,0,0,0],T=Tin,input = 0)
UinAsym = np.zeros((100000,2))
UinAsym[0:50,0] = 0.025#np.radians(5)
print(UinAsym)
TAsym, yAsym, xOutAsym = c.forced_response(systemAsym,T = Tin,U=np.transpose(UinAsym),X0=[0,0,0,0])

yAsym[2] = yAsym[2] * ((2*u_0)/span)
yAsym[3] = yAsym[3] * ((2*u_0)/span)
#print(TSym)


#plt.figure()
#plt.plot(TSym,ySym[0])
#plt.grid(True)
#plt.ylabel("u")
#plt.figure()
#plt.plot(TSym,ySym[1])
#plt.ylabel("alpha")
#plt.grid(True)
#plt.figure()
#plt.plot(TSym,ySym[2])
#plt.ylabel("theta")
#plt.grid(True)
#plt.figure()
#plt.plot(TSym,ySym[3])
#plt.ylabel("q")
#plt.grid(True)

plt.figure()
plt.plot(TSym,yAsym[0])
plt.grid(True)
plt.ylabel("sideslip")
plt.figure()
plt.plot(TSym,yAsym[1])
plt.ylabel("roll angle")
plt.grid(True)
plt.figure()
plt.plot(TSym,yAsym[2])
plt.ylabel("roll rate")
plt.grid(True)
plt.figure()
plt.plot(TSym,yAsym[3])
plt.ylabel("yaw rate")
plt.grid(True)

plt.show()