# Mesh Generatio Parameters

import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
import sympy as sp

# Parameters
y_plus = 50

# Fluid Properties
rho = 1.225 #kg/m^3
mu = 1.7894e-5 #kg/ms
nu = mu/rho #m^2/s

# Flow Properties
U = 10 #m/s
L = 1 #m
Re = rho*U*L/nu

#Skin Friction Coefficient
Cf = (2.0*math.log10(Re) - 6.5)**(-2.3)
tal_w = 0.5*rho*U**2*Cf
U_t = math.sqrt(tal_w/rho)

#Boundary Layer Height
yp = y_plus*mu/(rho*U_t)
yh = 2*yp

# Number of Layers
if y_plus > 30:
    N = 10
else:
    N = 25


#Sigma99 calculation
if Re < 5e5:
    sigma99 = 4.91*L/(Re**0.5)
else:
    sigma99 = 0.38*L/(Re**0.2)

#Growth Ratio
# Define the symbols
G = sp.Symbol('G')

# Define the equation
equation = sp.Eq(yh * (1 - G**N) / (1 - G) - sigma99, 0)

# Solve the equation for G
G_solution = sp.solve(equation, G)
# Filter out only the real solutions
real_G_solution = [sol.evalf() for sol in G_solution if sol.is_real]
G = round(real_G_solution[0], 2)