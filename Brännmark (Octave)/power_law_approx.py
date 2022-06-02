# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                         #
#      Power Law Approximation of Hill-Type Kinetics                      #
#                                                                         #
# OUTPUT: Returns the parameter values of the power law approximation of  #
#    Hill-type kinetic reactions (v29 and v33) of Brännmark et al's       #
#    insulin resistance signaling model [1]. We follow the procedure      #
#    detailed in Appendix B of Fortun et al [2] to approximate a non-     #
#    power law kinetic equation with power law kinetics.                  #
#                                                                         #
# References:                                                             #
#    [1] Brännmark C, Nyman E, Fagerholm S, Bergenholm L, Ekstrand E,     #
#           Cedersund G, Stralfors P (2013) Insulin signaling in type 2   #
#           diabetes: experimental and modeling analyses reveal           #
#           mechanisms of insulin resistance in human adipocytes. J Biol  #
#           Chem 288(14):9867--9880.                                      #
#           https://doi.org/10.1074/jbc.m112.432062                       #
#    [2] Fortun N, Mendoza E, Razon L, Lao A (2018) A deficiency-one      #
#           algorithm for power-law kinetic systems with reactant-        #
#           determined interactions. J Math Chem 56:2929--2962.           #
#           https://doi.org/10.1007/s10910-018-0925-2                     #
#                                                                         #
# Created: 28 May 2022                                                    #
# Last Modified: 2 June 2022                                              #
#                                                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Needed library
from sympy import *

# Define symbols
v29, k29, x28, n6, km6, x34 = symbols('v29 k29 x28 n6 km6 x34')
v33, k33, x36, x31, n9, km9 = symbols('v33 k33 x36 x31 n9 km9')

# Original equations
v29 = k29 * ((x28**n6)/(km6**n6 + x28**n6)) * x34
v33 = k33 * x36 * ((x31**n9)/(km9**n9 + x31**n9))

# Partial derivatives
pv29x28 = diff(v29, x28)
pv29x34 = diff(v29, x34)

pv33x36 = diff(v33, x36)
pv33x31 = diff(v33, x31)

# Evaluate partial derivatives
#    - Based on Matlab file from BioModels:
#         https://www.ebi.ac.uk/biomodels/BIOMD0000000448)
pv29x28 = pv29x28.subs([(k29, 36.93), (x28, 37.923), (n6, 2.137), (km6, 30.54), (x34, 29.576)])
pv29x34 = pv29x34.subs([(k29, 36.93), (x28, 37.923), (n6, 2.137), (km6, 30.54), (x34, 29.576)])

pv33x36 = pv33x36.subs([(k33, 0.1298), (x36, 96.047), (x31, 78.791), (n9, 0.9855), (km9, 5873.0)])
pv33x31 = pv33x31.subs([(k33, 0.1298), (x36, 96.047), (x31, 78.791), (n9, 0.9855), (km9, 5873.0)])

# Evaluate original reaction equations
v29 = v29.subs([(k29, 36.93), (x28, 37.923), (n6, 2.137), (km6, 30.54), (x34, 29.576)])
v33 = v33.subs([(k33, 0.1298), (x36, 96.047), (x31, 78.791), (n9, 0.9855), (km9, 5873.0)])

# Power law approximations
#    v29 = k29p * x28**p29 * x34**q29
#    v33 = k33p * x36**p33 * x31**q33

# Solve for power law approximation exponents
#    - Run the Octave file from BioModels for t = linspace(0,1000,1000)
x28 = 37.923
x34 = 29.576
x36 = 96.047
x31 = 78.791

print('p29 =')
p29 = pv29x28 * x28 / v29
print(p29)
print('\n')
print('q29 =')
q29 = pv29x34 * x34 / v29
print(q29)
print('\n')

print('p33 =')
p33 = pv33x36 * x36 / v33
print(p33)
print('\n')
print('q33 =')
q33 = pv33x31 * x31 / v33
print(q33)
print('\n')

# Solve for new rate constants
print('k29p =')
k29p = v29 / (x28**p29 * x34**q29)
print(k29p)
print('\n')

print('k33p =')
k33p = v33 / (x36**p33 * x31**q33)
print(k33p)
print('\n')




