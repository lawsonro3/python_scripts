# -*- coding: utf-8 -*-
"""
Created on Mon Nov 10 11:54:47 2014

@author: mlawson
"""
import nemohio as nio

sim = nio.Nemoh(directory='/Users/mlawson/Applications/nemoh/matlabRoutines/nemohIOTest',name='nemohIOTest')


sim.mesh.meshExec = '/Users/mlawson/bin/nemohMesh'
#sim.mesh.gdfFile = '/Users/mlawson/Applications/nemoh/matlabRoutines/nemohIOTest/geometry/rbWAMIT_shifted.gdf'
#sim.mesh.stlFile = '/Users/mlawson/Applications/nemoh/matlabRoutines/nemohIOTest/geometry/box.stl'
sim.mesh.vtpFile = '/Users/mlawson/Applications/nemoh/matlabRoutines/nemohIOTest/geometry/box.vtp'
sim.mesh.cg = [0.0, 0.0, 0.0]
sim.mesh.targetNumPanels = 100
sim.mesh.numBodies = 1

sim.writeId()

#sim.mesh.readSTL()
sim.mesh.readVTP()
#sim.mesh.readGDF()
sim.mesh.writeNemohMesh()
sim.mesh.writeMeshCal()
sim.mesh.writeVtp()
sim.mesh.run()

sim.sim.writeInput()