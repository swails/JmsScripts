#!/usr/bin/env python 

# This program takes a file containing <frame# phi psi> and creates Ramachandran plots
# using the numpy and matplotlib modules in python 2.6
# by Asim Okur (asimokur@gmail.com) 03/10/2011



import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

kb = 1.9858775e-3 	# Boltzmann const
T = 300.0		# Simulation Temperature

name = 'dihedrals.300' # Input filename - output file name will be created based on input
phipsi = np.loadtxt('%s.dat'%name)

H, xedges, yedges = np.histogram2d(phipsi[:,3], phipsi[:,2], bins=72, range = [[-180,180],[-180,180]])
extent = (xedges[0],xedges[-1],yedges[0],yedges[-1]) #get the xmin,xmax values for 
H = -kb * T * np.log(H/H.max())    # Convert populations to free energy; highest populated bin at 0 kcal/mol

# Plot the free energy landscape:
plt.imshow(H, extent=extent, aspect='auto', origin='lower', interpolation='nearest', cmap=cm.jet) 

# Plot black contour lines
levels = np.arange(0,5,0.5) 	# Contour levels from 0 to 5 every 0.5 kcal/mol
plt.contour(H, colors = 'k', extent=extent, levels = levels, aspect='auto', origin='lower' ) 
# Fill in the surface with solid color for each contour 
plt.contourf(H, extent=extent, aspect='auto', levels = levels, origin='lower', cmap=cm.jet)

cbar = plt.colorbar(format='%4.2f')        			# Add colorbar
cbar.set_label('Free energy (kcal/mol)',fontsize='large') 	# Label for the colorbar


# Label and format the axes:
plt.ylabel('$\psi$',fontsize='large')
plt.xlabel('$\phi$',fontsize='large')
plt.xticks(np.arange(-180,180.01,60)) # from -180 to 180 every 60 degrees
plt.yticks(np.arange(-180,180.01,60)) # from -180 to 180 every 60 degrees

plt.title('Dihedral free energy profile',fontsize='x-large',fontweight='bold') # Chart title, either update manually or tie in to the input filename

plt.savefig('%s.png'%name) # save fig as png file 

