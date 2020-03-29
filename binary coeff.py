# author Dominik Capkovic 
# contact: domcapkovic@gmail.com; https://www.linkedin.com/in/dominik-čapkovič-b0ab8575/
# GitHub: https://github.com/kilimetr
# With using virial coeff B the binary coeff of mixture can be evaluated
# Biphenyl & Benzoic acid

import numpy as np

Tc      = np.empty(4) # [K]
pc      = np.empty(4) # [kPa]
Tb      = np.empty(2)
Trb     = np.empty(4)
omegaEd = np.empty(8)
omegaLK = np.empty(8)
alpha   = np.empty(4)
a       = np.empty(4)
b       = np.empty(4)
deltaB  = np.empty(4)
f0      = np.empty(4)
f1      = np.empty(4)
B       = np.empty(8)
Bmix    = np.empty(8)

# B2 = np.empty(4)

Tc[0] = 515.7+273.15 # bi
Tc[1] = 780 # +-20
Tc[2] = 479+273.15 # ba
Tc[3] = 741.85

pc[0] = 4050 # bi
pc[1] = 3500 # +-6
pc[2] = 4559.625 # ba
pc[3] = 5073.01

Tb[0] = 527 # bi
Tb[1] = 522 # ba


T = 150.5 + 273.15 # average temp
R = 8.314


# bi
Trb[0] = Tb[0] / Tc[0]
Trb[1] = Tb[0] / Tc[1]

omegaEd[0] = 3/7 / (Tc[0]/Tb[0] - 1) * np.log10(pc[0]/101.325) - 1 # pc v MPa = ?
# omegaEd[1] = 3/7 / (Tc[1]/Tb[0] - 1) * np.log10(pc[0]/101.325) - 1
# omegaEd[2] = 3/7 / (Tc[0]/Tb[0] - 1) * np.log10(pc[1]/101.325) - 1
omegaEd[3] = 3/7 / (Tc[1]/Tb[0] - 1) * np.log10(pc[1]/101.325) - 1

# omegaLK[0] = np.log(101.325/pc[0]) - 5.92714 + 6.09648/Trb[0] + 1.28862*np.log(Trb[0]) - 0.169347*pow(Trb[0], 6) / \
# 			 (15.2518 - 15.6875/Trb[0] - 13.4721*np.log(Trb[0]) + 0.43577*pow(Trb[0], 6))
# omegaLK[1] = np.log(101.325/pc[1]) - 5.92714 + 6.09648/Trb[0] + 1.28862*np.log(Trb[0]) - 0.169347*pow(Trb[0], 6) / \
# 			 (15.2518 - 15.6875/Trb[0] - 13.4721*np.log(Trb[0]) + 0.43577*pow(Trb[0], 6))
# omegaLK[2] = np.log(101.325/pc[0]) - 5.92714 + 6.09648/Trb[1] + 1.28862*np.log(Trb[1]) - 0.169347*pow(Trb[1], 6) / \
# 			 (15.2518 - 15.6875/Trb[1] - 13.4721*np.log(Trb[1]) + 0.43577*pow(Trb[1], 6))
# omegaLK[3] = np.log(101.325/pc[1]) - 5.92714 + 6.09648/Trb[1] + 1.28862*np.log(Trb[1]) - 0.169347*pow(Trb[1], 6) / \
# 			 (15.2518 - 15.6875/Trb[1] - 13.4721*np.log(Trb[1]) + 0.43577*pow(Trb[1], 6))

alpha[0] = pow(1+(1-pow(T/Tc[0], 0.5)) * (0.37464 + 1.54226*omegaEd[0] - 0.26992*pow(omegaEd[0], 2)), 2)
alpha[1] = pow(1+(1-pow(T/Tc[1], 0.5)) * (0.37464 + 1.54226*omegaEd[3] - 0.26992*pow(omegaEd[3], 2)), 2)

a[0] = 0.4572355 * pow(R, 2) * Tc[0] * alpha[0] / (pc[0]*1e+03)
a[1] = 0.4572355 * pow(R, 2) * Tc[1] * alpha[1] / (pc[1]*1e+03)

b[0] = 0.077796 * R * Tc[0] / (pc[0]*1e+03)
b[1] = 0.077796 * R * Tc[1] / (pc[1]*1e+03)

deltaB[0] = a[0]/pow(T/Tc[0], 6) - b[0]/pow(T/Tc[0], 8) # both too big?
deltaB[1] = a[1]/pow(T/Tc[1], 6) - b[1]/pow(T/Tc[1], 8)

f0[0] = 0.1445 - 0.33/(T/Tc[0]) - 0.1385/pow(T/Tc[0], 2) - 0.0121/pow(T/Tc[0], 3) - 6.07e-03/pow(T/Tc[0], 8)
f0[1] = 0.1445 - 0.33/(T/Tc[1]) - 0.1385/pow(T/Tc[1], 2) - 0.0121/pow(T/Tc[1], 3) - 6.07e-03/pow(T/Tc[1], 8)

f1[0] = 0.0637 + 0.331/pow(T/Tc[0], 2) - 0.423/pow(T/Tc[0], 3) - 0.008/pow(T/Tc[0], 8)
f1[1] = 0.0637 + 0.331/pow(T/Tc[1], 2) - 0.423/pow(T/Tc[1], 3) - 0.008/pow(T/Tc[1], 8)

B[0] = R*Tc[0]/pc[0] * (f0[0] + omegaEd[0]*f1[0])# + deltaB[0]
B[1] = R*Tc[1]/pc[1] * (f0[1] + omegaEd[0]*f1[1])# + deltaB[1]
B[2] = R*Tc[0]/pc[0] * (f0[0] + omegaEd[3]*f1[0])# + deltaB[0]
B[3] = R*Tc[1]/pc[1] * (f0[1] + omegaEd[3]*f1[1])# + deltaB[1]


# ba
Trb[2] = Tb[1] / Tc[2]
Trb[3] = Tb[1] / Tc[3]

omegaEd[4] = 3/7 / (Tc[2]/Tb[1] - 1) * np.log10(pc[2]/101.325) - 1
# omegaEd[5] = 3/7 / (Tc[3]/Tb[1] - 1) * np.log10(pc[2]/101.325) - 1
# omegaEd[6] = 3/7 / (Tc[2]/Tb[1] - 1) * np.log10(pc[3]/101.325) - 1
omegaEd[7] = 3/7 / (Tc[3]/Tb[1] - 1) * np.log10(pc[3]/101.325) - 1

# omegaLK[4] = np.log(101.325/pc[2]) - 5.92714 + 6.09648/Trb[2] + 1.28862*np.log(Trb[2]) - 0.169347*pow(Trb[2], 6) / \
# 			 (15.2518 - 15.6875/Trb[2] - 13.4721*np.log(Trb[2]) + 0.43577*pow(Trb[2], 6))
# omegaLK[5] = np.log(101.325/pc[3]) - 5.92714 + 6.09648/Trb[2] + 1.28862*np.log(Trb[2]) - 0.169347*pow(Trb[2], 6) / \
# 			 (15.2518 - 15.6875/Trb[2] - 13.4721*np.log(Trb[2]) + 0.43577*pow(Trb[2], 6))
# omegaLK[6] = np.log(101.325/pc[2]) - 5.92714 + 6.09648/Trb[3] + 1.28862*np.log(Trb[3]) - 0.169347*pow(Trb[3], 6) / \
# 			 (15.2518 - 15.6875/Trb[3] - 13.4721*np.log(Trb[3]) + 0.43577*pow(Trb[3], 6))
# omegaLK[7] = np.log(101.325/pc[3]) - 5.92714 + 6.09648/Trb[3] + 1.28862*np.log(Trb[3]) - 0.169347*pow(Trb[3], 6) / \
# 			 (15.2518 - 15.6875/Trb[3] - 13.4721*np.log(Trb[3]) + 0.43577*pow(Trb[3], 6))

alpha[2] = pow(1+(1-pow(T/Tc[2], 0.5)) * (0.37464 + 1.54226*omegaEd[4] - 0.26992*pow(omegaEd[4], 2)), 2)
alpha[3] = pow(1+(1-pow(T/Tc[3], 0.5)) * (0.37464 + 1.54226*omegaEd[7] - 0.26992*pow(omegaEd[7], 2)), 2)

a[2] = 0.4572355 * pow(R, 2) * Tc[2] * alpha[2] / (pc[2]*1e+03)
a[3] = 0.4572355 * pow(R, 2) * Tc[3] * alpha[3] / (pc[3]*1e+03)

b[2] = 0.077796 * R * Tc[2] / (pc[2]*1e+03)
b[3] = 0.077796 * R * Tc[3] / (pc[3]*1e+03)

deltaB[2] = a[2]/pow(T/Tc[2], 6) - b[2]/pow(T/Tc[2], 8) # both too big?
deltaB[3] = a[3]/pow(T/Tc[3], 6) - b[3]/pow(T/Tc[3], 8)

f0[2] = 0.1445 - 0.33/(T/Tc[2]) - 0.1385/pow(T/Tc[2], 2) - 0.0121/pow(T/Tc[2], 3) - 6.07e-03/pow(T/Tc[2], 8)
f0[3] = 0.1445 - 0.33/(T/Tc[3]) - 0.1385/pow(T/Tc[3], 2) - 0.0121/pow(T/Tc[3], 3) - 6.07e-03/pow(T/Tc[3], 8)

f1[2] = 0.0637 + 0.331/pow(T/Tc[2], 2) - 0.423/pow(T/Tc[2], 3) - 0.008/pow(T/Tc[2], 8)
f1[3] = 0.0637 + 0.331/pow(T/Tc[3], 2) - 0.423/pow(T/Tc[3], 3) - 0.008/pow(T/Tc[3], 8)

B[4] = R*Tc[2]/pc[2] * (f0[2] + omegaEd[4]*f1[2])# + deltaB[2]
B[5] = R*Tc[3]/pc[3] * (f0[3] + omegaEd[4]*f1[3])# + deltaB[3]
B[6] = R*Tc[2]/pc[2] * (f0[2] + omegaEd[7]*f1[2])# + deltaB[2]
B[7] = R*Tc[3]/pc[3] * (f0[3] + omegaEd[7]*f1[3])# + deltaB[3]

# B2[0] = 0.073 + 0.46/(T/Tc[0]) - 0.5*pow(T/Tc[0], -2) - 0.097*pow(T/Tc[0], -3) - 0.0073*pow(T/Tc[0], -8)


# print(Trb)
print(omegaEd)
# print(omegaLK) # without any sense
# print(alpha)
# print(a)
# print(b)
print(deltaB)
# print(f0)
# print(f1)
print(B)
# print(B2)




