author: michael lawson

This module provides the functionality to:
*Read hydrodynamic coefficients from AQWA ".lis" files
*Plot hydrodynamic coefficients
*Write hydrodynamic data to .mat files for use in WEC-SIm
*Write hydrodynamic data to HDF5 file format for use in other applications

Descrioption of files in this folder:
*aqwaio.py: python module
*aqwaio-example.py: example of how to use the awqaio module
*aqwa-example-data.lis: example data of an aqwa run. no hydrodynamic interactions for a three body flapping device that consists of a frame and two flaps
*compareBEM.m: a matlab file used during developement - normal users can ignore this file
aqwa-manually-generated.mad: a matlab file used during developement - normal users can ignore this file

This module is currently under active development, bug identification and fixed from the community are weclcome! Please post bugs and feature requests on GitHub.

Notes:
*Only works for AQWA simulations without body to body interactions turned on

