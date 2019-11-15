from __future__ import print_function
from __future__ import division

import cantera as ct
import numpy as np

print("Running Cantera Version: " + str(ct.__version__))

#Define the gas-mixutre and kinetics
#In this case, we are choosing a GRI3.0 gas
gas = ct.Solution('ALZ_F-MechC.cti')

print("Species in mech:")
for i, specie in enumerate(gas.species()):
    print(str(i) + '. ' + str(specie))

# Create a premixed mixture 
# Inlet Temperature in Kelvin and Inlet Pressure in Pascals
gas.TPX = 600.0, ct.one_atm, 'CF4:0.031, CH4:0.092, O2:0.184, N2:0.693'

# Domain width in metres
width = 0.07

# Create the flame object
flame = ct.FreeFlame(gas, width=width)

# Define tolerances for the solver
flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1)

# Define logging level
loglevel = 1

flame.solve(loglevel=loglevel, auto=True)

Su0 = flame.u[0]
print("Flame Speed is: {:.2f} cm/s".format(Su0*100))

# Import plotting modules and define plotting preference
import matplotlib.pylab as plt

plt.rcParams['axes.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['figure.figsize'] = (8,6)

plt.rcParams['figure.autolayout'] = True

plt.figure()

# Extract concentration data
X_CH4 = flame.X[13]
X_CO2 = flame.X[15]
X_H2O = flame.X[5]
X_CF4 = flame.X[59]

plt.plot(flame.grid*100, X_CH4, '-o', label=r'$CH_{4}$')
plt.plot(flame.grid*100, X_CO2, '-s', label=r'$CO_{2}$')
plt.plot(flame.grid*100, X_H2O, '-<', label=r'$H_{2}O$')
plt.plot(flame.grid*100, X_CF4, '->', label=r'$CF_{4}$')

plt.legend(loc=2)
plt.xlabel('Distance (cm)')
plt.ylabel('MoleFractions');
plt.show()
