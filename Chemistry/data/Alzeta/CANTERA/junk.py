from __future__ import print_function
from __future__ import division

import cantera as ct
import numpy as np
import os

def getCellVols(grid):
    Npts = len(grid)
    vols = []
    for i in range(Npts):
        if i==0:
            zL = grid[0]
        else:
            zL = 0.5*(grid[i-1] + grid[i])

        if i==Npts-1:
            zR = grid[i]
        else:
            zR = 0.5*(grid[i] + grid[i+1])

        vols.append(zR - zL)

    return vols


print("Running Cantera Version: " + str(ct.__version__))

#Define the gas-mixutre and kinetics
#In this case, we are choosing a GRI3.0 gas
gas = ct.Solution('ALZ_F-MechC.cti')

# print("Species in mech:")
# for i, specie in enumerate(gas.species()):
#     print(str(i) + '. ' + str(specie))

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

# vols = getCellVols(flame.grid)
# frop = flame.forward_rates_of_progress
# rrop = flame.reverse_rates_of_progress

# Npts = len(flame.grid)
# Ncells = len(vols)
# tot = 0.0
# for i in range(Ncells):
#     tot += (frop[100,i] - rrop[100,i]) * vols[i]

# print('tot ' + str(tot))

print("T " + str(gas.T))

element = "F"
diagram = ct.ReactionPathDiagram(gas, element)

diagram.title = 'Reaction path diagram following {0}'.format(element)
diagram.label_threshold = 0.01
dot_file = 'rxnpath.dot'
img_file = 'rxnpath.png'
img_path = os.path.join(os.getcwd(), img_file)

print('details: ' + str(diagram.show_details))
diagram.show_details = True
print('details: ' + str(diagram.show_details))

diagram.write_dot(dot_file)
print(diagram.get_data())
print("Wrote graphviz input file to '{0}'.".format(os.path.join(os.getcwd(), dot_file)))
os.system('dot {0} -Tpng -o{1} -Gdpi=70'.format(dot_file, img_file))
print("Wrote graphviz output file to '{0}'.".format(img_path))

# Su0 = flame.u[0]
# print("Flame Speed is: {:.2f} cm/s".format(Su0*100))

# # Import plotting modules and define plotting preference
# import matplotlib.pylab as plt

# plt.rcParams['axes.labelsize'] = 14
# plt.rcParams['xtick.labelsize'] = 12
# plt.rcParams['ytick.labelsize'] = 12
# plt.rcParams['legend.fontsize'] = 10
# plt.rcParams['figure.figsize'] = (8,6)

# plt.rcParams['figure.autolayout'] = True

# plt.figure()

# # Extract concentration data
# X_CH4 = flame.X[13]
# X_CO2 = flame.X[15]
# X_H2O = flame.X[5]
# X_CF4 = flame.X[59]

# plt.plot(flame.grid*100, X_CH4, '-o', label=r'$CH_{4}$')
# plt.plot(flame.grid*100, X_CO2, '-s', label=r'$CO_{2}$')
# plt.plot(flame.grid*100, X_H2O, '-<', label=r'$H_{2}O$')
# plt.plot(flame.grid*100, X_CF4, '->', label=r'$CF_{4}$')

# plt.legend(loc=2)
# plt.xlabel('Distance (cm)')
# plt.ylabel('MoleFractions');
# plt.show()
