import numpy as np
import csv

class NREL_5MW(object):
    '''NREL 5MW data from Table 3.1 of "Definition of a 5-MW Reference Wind Turbine for Offshore System Development"
    '''
    def __init__(self):
        self.headers = np.array(['node (-)','r_nodes (meters)','aero_twist (degrees)','dr_nodes (meters)','chord (meters)', 'airfoil_name (-)'])
        self.node = np.linspace(1,17,17)
        self.r_nodes = np.array([2.8667, 5.6, 8.3333, 11.75, 15.85, 19.95, 24.05, 28.15, 32.25, 36.35, 40.45, 44.55, 48.65, 52.75, 56.1667, 58.9, 61.6333])
        self.aero_twist = np.array([13.308, 13.308, 13.308, 13.308, 11.48, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.37, 0.106])
        self.aero_twist = np.array([2.7333, 2.7333, 2.7333, 4.1, 4.1, 4.1, 4.1, 4.1, 4.1, 4.1, 4.1, 4.1, 4.1, 4.1, 2.7333, 2.7333, 2.7333])
        self.chord = np.array([3.542, 3.854, 4.167, 4.557, 4.652, 4.458, 4.249, 4.007, 3.748, 3.502, 3.256, 3.01, 2.764, 2.518, 2.313, 2.086, 1.419])
        self.aero_twist = np.array(['Cylinder1.dat', 'Cylinder1.dat', 'Cylinder2.dat', 'DU40_A17.dat', 'DU35_A17.dat', 'DU35_A17.dat', 'DU30_A17.dat', 'DU25_A17.dat', 'DU25_A17.dat', 'DU21_A17.dat', 'DU21_A17.dat', 'NACA64_A17.dat', 'NACA64_A17.dat', 'NACA64_A17.dat', 'NACA64_A17.dat', 'NACA64_A17.dat', 'NACA64_A17.dat'],dtype='S')
        self.file_name = 'nrel_5mw_aero_properties.h5'
        self.out_file = 'nrel_5mw_aero_properties.txt'
