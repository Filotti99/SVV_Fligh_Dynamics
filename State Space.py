import numpy as np
import control.matlab as ml
import control as c
import matplotlib.pyplot as plt
import cmath
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
chord = 2.0569
S = 30.00
Aspect_ratio = span**2/S

D_b = span/V
D_c = chord / V

mu_b = mass / (density*S*span)
#print(mu_b)
mu_c = mass / (density * S * chord)
Kx = np.sqrt(0.019)
Ky = np.sqrt(1.25*1.114)
Kz = np.sqrt(0.042)
Kxz = 0.002

CL = 2*Weight/(density*V**2*S)
#print(CL)
CD = CD_0 + ((CL_alpha*alpha_0)**2/(np.pi*Aspect_ratio*oswald_factor))

### Stability Derivatives ###
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

### Equation of Motion Matrices - V1 (Outdated) ###
#C1sym = np.matrix([[(-2*mu_c*chord)/(V**2),					0,								0,											0],
#				   [0,										(Cz_alpha_dot-2*mu_c)*chord/V,		0,											0],
#				   [0,										0,								-chord/V,									0],
#				   [0,										Cm_alpha_dot*chord/V,				0,										-2*mu_c*Ky**2*chord**2/(V**2)]])

#C2sym = np.matrix([[Cx_u/V,								Cx_alpha,						Cz_0,										Cx_q*chord/V],
#				   [Cz_u/V,									Cz_alpha,						-Cx_0,										(Cz_q+2*mu_c)*chord/V],
#				   [0,										0,								0,											chord/V],
#				   [Cm_u/V,									Cm_alpha,						0,											Cm_q*chord/V]])

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

### Equation of Motion Matrices - V2 ###
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
# Aasym[0] = [0,0,0,2/D_b]
Basym = -np.linalg.inv(C1asym)*C3asym


def PrintAB(ShouldPrint):
	'''
	Can be used to print the A and B matrices of both the symmetric and asymmetric motions
	:param ShouldPrint: True/False - determines if the matrices will be printed
	:return: Nothing (it does print stuff though :D )
	'''
	if ShouldPrint is True:
		print("symmetric: ")
		print("A: ", Asym)
		print("B: ", Bsym)
		print(V / chord)

		print("Asymmetric: ")
		print("A: ", Aasym)
		print("B: ", Basym)

def PrintStabilityDerivatives(ShouldPrint):
	'''
	Can be used to print a selection of stability derivatives
	:param ShouldPrint: if True, the derivatives will be printed
	:return: None (it does print stuff though :D )
	'''
	if ShouldPrint is True:
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

systemSym = ml.ss(Asym, Bsym, Csym, Dsym)
systemAsym = ml.ss(Aasym, Basym, Casym, Dasym)

Tin = np.arange(0,20,0.01)
Uin = np.zeros_like(Tin)
Uin[0:100] = np.radians(10)

TSym, ySym, xOut = c.forced_response(systemSym, T = Tin, U = Uin, X0 = [0, alpha_0, theta_0, 0])
# ySym, TSym = ml.step(systemSym,X0=[0,alpha_0,theta_0,0],T = Tin)

ySym[0] = ySym[0]*u_0 + u_0
ySym[3] = ySym[3]*(u_0 / chord)

#yAsym, TAsym = ml.impulse(systemAsym,X0=[0,0,0,0],T=Tin,input = 0)
UinAsym = np.zeros((2000,2))
UinAsym[0:100,0] = np.radians(10)

TAsym, yAsym, xOutAsym = c.forced_response(systemAsym,T = Tin,U=np.transpose(UinAsym),X0=[0,0,0,0])

yAsym[2] = yAsym[2] * ((2*u_0)/span)
yAsym[3] = yAsym[3] * ((2*u_0)/span)



def PrintEigvals(ShouldPrint):
	'''
	Can be used to print the eigenvalues of the symmetric and asymmetric systems
	:param ShouldPrint: if True, the eigenvalues will be printed
	:return: None (it does print stuff though :D )
	'''
	if ShouldPrint is True:
		print("Eigenvalues symmetric case: \n", np.linalg.eigvals(systemSym.A))
		print("Eigenvalues asymmetric case: \n", np.linalg.eigvals(systemAsym.A))


def PlotSym(ShouldPlot):
	'''
	Can be used to plot the response of the symmetric system to a disturbance
	:param ShouldPlot: if True, will plot the response
	:return: None (it does plot stuff though :D )
	'''
	if ShouldPlot is True:
		plt.figure()
		plt.subplot(2, 2, 1)
		plt.plot(TSym, ySym[0])
		plt.grid()
		plt.ylabel("u")
		plt.subplot(2, 2, 2)
		plt.plot(TSym, ySym[1])
		plt.grid()
		plt.ylabel("alpha")
		plt.subplot(2, 2, 3)
		plt.plot(TSym, ySym[2])
		plt.grid()
		plt.ylabel("theta")
		plt.subplot(2, 2, 4)
		plt.plot(TSym, ySym[3])
		plt.grid()
		plt.ylabel("q")


def PlotAsym(ShouldPlot):
	'''
	Can be used to plot the response of the symmetric system to a disturbance
	:param ShouldPlot: if True, will plot the response
	:return: None (it does plot stuff though :D )
	'''
	if ShouldPlot is True:
		plt.figure()
		plt.subplot(2, 2, 1)
		plt.plot(TSym, yAsym[0])
		plt.grid()
		plt.ylabel("sideslip")
		plt.subplot(2, 2, 3)
		plt.plot(TSym, yAsym[1])
		plt.grid()
		plt.ylabel("roll angle")
		plt.subplot(2, 2, 2)
		plt.plot(TSym, yAsym[2])
		plt.grid()
		plt.ylabel("roll rate")
		plt.subplot(2, 2, 4)
		plt.plot(TSym, yAsym[3])
		plt.grid()
		plt.ylabel("yaw rate")

PrintAB(False)
PrintStabilityDerivatives(False)
PrintEigvals(True)
PlotSym(True)
PlotAsym(True)
plt.show()

def ShortPeriodOscillation():
	print("Short Period Oscillation")
	sa1 = -2 * mu_c * Ky ** 2 * (Cz_alpha_dot - 2 * mu_c)
	sb1 = -2 * mu_c * Ky ** 2 * Cz_alpha + Cm_q * (Cz_alpha_dot - 2 * mu_c) - Cm_alpha_dot * (Cz_q + 2 * mu_c)
	sc1 = Cz_alpha * Cm_q - Cm_alpha * (Cz_q + 2 * mu_c)
	#Ja echt Ivo wat is dit voor form?
	print("Lambda1 ", (V / chord) * (-sb1 - cmath.sqrt(sb1 ** 2 - 4 * sa1 * sc1)) / (2 * sa1))
	print("Lambda2 ", (V / chord) * (-sb1 + cmath.sqrt(sb1 ** 2 - 4 * sa1 * sc1)) / (2 * sa1))
	print()

def PhugoidMotion():
	print("Phugoid Motion")
	sa2 = 2 * mu_c * ((Cz_alpha * Cm_q) - 2 * mu_c * Cm_alpha)
	sb2 = 2 * mu_c * ((Cx_u * Cm_alpha) - Cm_u * Cx_alpha) + Cm_q * ((Cz_u * Cx_alpha) - (Cx_u * Cz_alpha))
	sc2 = Cz_0 * ((Cm_u * Cz_alpha) - (Cz_u * Cm_alpha))

	lambda21 = (-sb2 - cmath.sqrt(sb2 ** 2 - 4 * sa2 * sc2)) / (2 * sa2)
	lambda22 = (-sb2 + cmath.sqrt(sb2 ** 2 - 4 * sa2 * sc2)) / (2 * sa2)
	print("Lambda1 ", ((V / chord) * lambda21))
	print("Lambda2 ", ((V / chord) * lambda22))
	print()





def HeavilyDampedAperiodicRollingMotion():
	Lambda_b1= V/span * Cl_p/(4*mu_b*Kx**2)
	print("Heavily Damped Aperiodic Rolling Motion:")
	print("Lambda_b1: ", Lambda_b1)
	print()

def DutchRollMotion():
	A = 8*mu_b**2 * Kz**2
	B = -2*mu_b*(Cn_r + 2*Kz**2 *Cy_beta)
	C = 4*mu_b*Cn_beta + Cy_beta*Cn_r
	Lambda1 = V/span * (-B + cmath.sqrt(B ** 2 - 4 * A * C)) / (2 * A)
	Lambda2 = V/span * (-B - cmath.sqrt(B ** 2 - 4 * A * C)) / (2 * A)
	print("Dutch Roll Motion")
	print("Lambda1: ", Lambda1)
	print("Lambda2: ", Lambda2)
	print()

def AperiodicSpiralMotion():
	Lambda_b4 = V/span* (2*CL*(Cl_beta*Cn_r - Cn_beta*Cl_r))/(Cl_p*(Cy_beta*Cn_r + 4*mu_b*Cn_beta) - Cn_p*(Cy_beta*Cl_r + 4*mu_b*Cl_beta))
	print("AperiodicSpiralMotion")
	print("Lambda_b4: ", Lambda_b4)
	print()

def DutchRollMotionAndAperiodicRollingMotion():
	A = 4*mu_b**2 *(Kx**2 * Kz**2 - Kxz**2)
	B = -mu_b*((Cl_r+Cn_p)*Kxz + Cn_r*Kx**2 + Cl_p*Kz**2)
	C = 2*mu_b*(Cl_beta*Kxz+Cn_beta*Kx**2) + 1/4*(Cl_p*Cn_r - Cn_p*Cl_r)
	D = 1/2 * (Cl_beta*Cn_p - Cn_beta*Cl_p)
	Lambda1, Lambda2, Lambda3 = V/span * np.roots([A, B, C, D])

	print("Dutch Roll Motion and Aperiodic Rolling Motion")
	print("Lambda1: ", Lambda1)
	print("Lambda2: ", Lambda2)
	print("Lambda3: ", Lambda3)
	print()

ShortPeriodOscillation()
PhugoidMotion()
HeavilyDampedAperiodicRollingMotion()
DutchRollMotion()
AperiodicSpiralMotion()
DutchRollMotionAndAperiodicRollingMotion()


