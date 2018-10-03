import numpy as np
import matplotlib.pyplot as plt
plt.close('all')

# From section 3.8.3 of wind energy explained
# Prandlt tip loss calc
B = 3 # number of blades
R = 1 # blade length
phi = np.deg2rad(10) # relative wind angle
r = np.linspace(0,R,100)
F = 2/np.pi * np.arccos(np.exp(-((B/2)*(1-(r/R)))/((r/R)*np.sin(phi))))

plt.figure(num='Tip loss for phi = %2.1f deg and %d blades' % (np.rad2deg(phi), B))
plt.plot(r,F)
plt.xlabel('Non-Dimensional Blade Radius (r/R)')
plt.ylabel('Tip Loss Factor')
