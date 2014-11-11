"""
Created on Tue Nov  4 15:28:52 2014

@author: mlawson

This module reads data from Nemoh simulations and convert from GDF to Nemoh mesh format
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import vtk
from vtk.util.numpy_support import vtk_to_numpy

def mkVtkIdList(it):
    vil = vtk.vtkIdList()
    for i in it:
        vil.InsertNextId(int(i))
    return vil
    
class Nemoh(object):
    def __init__(self,directory,name):        
        self.baseDir = directory
        self.name = name
        self.dir = directory + os.path.sep + name
        del directory, name
        
        self.files = {}
        self.files['ID.dat'] = self.baseDir + os.path.sep + 'ID.dat'

        self.sim  =  NemohSimulation(self)   
        self.mesh =  NemohMesh(self)
        
    def writeId(self):        
        with open(self.files['ID.dat'],'w') as fid:
            fid.write(str(len(self.name)))
            fid.write('\n')
            fid.write(self.name)  

#        self.resultsFiles['RadiationCoefficients'] = self.dirs['results'] + 'RadiationCoefficieaddedMass'
#        self.resultsFiles['IRF'] = self.dirs['results'] + 'IRF.tec' 
#        self.resultsFiles['index'] = self.dirs['results'] + 'index.dat'   
#        self.resultsFiles['Forces'] = self.dirs['results'] + 'Forces.dat'   
#        self.resultsFiles['FKForceTec'] = self.dirs['results'] + 'FKForce.tec' 
#        self.resultsFiles['Fe'] = self.dirs['results'] + 'Fe.dat' 
#        self.resultsFiles['ExcitationForce'] = self.dirs['results'] + 'ExcitationForce.tec' 
#        self.resultsFiles['Diffraction'] = self.dirs['results'] + 'DiffractionForce.tec' 
#        self.resultsFiles['CM'] = self.dirs['results'] + 'CM.dat' 
#        self.resultsFiles['CA'] = self.dirs['results'] + 'CA.dat'         
    def clean(self):
        pass

class NemohMesh(object):
    def __init__(self,sim):
        self.dir = sim.dir + os.path.sep + 'mesh'
        self.baseDir = sim.dir + os.path.sep + '..'
        if os.path.exists(self.dir) is False:
            os.mkdir(self.dir)
        
        self.files = {}
        self.files['NemohMesh.tec'] = self.dir + os.path.sep + 'Mesh.tec'
        self.files['L12.dat'] = self.dir + os.path.sep + 'L12.dat'
        self.files['L10.dat'] = self.dir + os.path.sep + 'L10.dat'
        self.files['Kochin.dat'] = self.dir + os.path.sep + 'Kochin.dat'
        self.files['KH.dat'] = self.dir + os.path.sep + 'KH.dat'
        self.files['Integration.dat'] = self.dir + os.path.sep + 'Integration.dat'
        self.files['Inertia_hull.dat'] = self.dir + os.path.sep + 'Inertia_hull.dat'
        self.files['Hydrostatics.dat'] = self.dir + os.path.sep + 'Hydrostatics.dat'
        self.files['CG_hull.dat'] = self.dir + os.path.sep + 'CG_hull.dat'
        self.files['Freesurface.dat'] = self.dir + os.path.sep + 'Freesurface.dat'
        self.files['Description_Wetted.tec'] = self.dir + os.path.sep + 'Description_Wetted.tec'
        self.files['Description_Full.tec'] = self.dir + os.path.sep + 'Description_Full.tec'
        self.files['ID.dat'] = sim.dir + os.path.sep + 'ID.dat'
        self.files['Mesh.cal'] = sim.baseDir + os.path.sep + 'Mesh.cal'
        self.files['Mesh'] = None
        
        self.faces = []
        self.cords = []
        
        self._gdfFile = None
        self._stlFile = None
        self._vtpFile = None
        self._cg = None
        self._targedNumPanels = None
        self._numBodies = None
        self._prog = None
        self._meshExec = None

    @property
    def meshExec(self):
        return self._meshExec
    @meshExec.setter
    def meshExec(self,fName):
        temp1, temp2 = os.path.split(fName)
        self._meshExec = fName

    @property
    def prog(self):
        return self._prog
    @prog.setter
    def prog(self,fName):
        self._prog = fName
        if not os.path.exists(fName):
            raise Exception('The nemoh program does not exist on your operating system search path')
    
    @property
    def gdfFile(self):
        return self._gdfFile
    @gdfFile.setter
    def gdfFile(self,fName):
        self._gdfFile = fName
        self.files['GDF'] = self._gdfFile
        temp, self.meshName = os.path.split(self._gdfFile[:-3])
        self.files['vtk'] = self.dir + os.path.sep + self.meshName + 'vtp'
        self.files['Mesh'] = self.dir + os.path.sep + self.meshName + 'nemohMesh'
        self.files['MeshLog'] = self.dir + os.path.sep + self.meshName + 'nemohMesh.log'
        temp, self.files['Mesh-noPath'] = os.path.split(self.files['Mesh'])
        
    @property
    def stlFile(self):
        return self._stlFile
    @stlFile.setter
    def stlFile(self,fName):
        self._stlFile = fName
        self.files['stl'] = self._stlFile
        temp, self.meshName = os.path.split(self._stlFile[:-3])
        self.files['vtk'] = self.dir + os.path.sep + self.meshName + 'vtp'
        self.files['Mesh'] = self.dir + os.path.sep + self.meshName + 'nemohMesh'
        self.files['MeshLog'] = self.dir + os.path.sep + self.meshName + 'nemohMesh.log'
        temp, self.files['Mesh-noPath'] = os.path.split(self.files['Mesh'])
        
    @property
    def vtpFile(self):
        return self._vtpFile
    @vtpFile.setter
    def vtpFile(self,fName):
        self._vtpFile = fName
        self.files['vtp'] = self._vtpFile
        temp, self.meshName = os.path.split(self._vtpFile[:-3])
        self.files['vtk'] = self.dir + os.path.sep + self.meshName + 'vtp'
        self.files['Mesh'] = self.dir + os.path.sep + self.meshName + 'nemohMesh'
        self.files['MeshLog'] = self.dir + os.path.sep + self.meshName + 'nemohMesh.log'
        temp, self.files['Mesh-noPath'] = os.path.split(self.files['Mesh'])
        
    @property
    def cg(self):
        return self._cg        
    @cg.setter
    def cg(self,cg):
        if np.shape(cg) == (3,) is False:
            raise Exception('cg input error')
        self._cg = cg
        
    @property 
    def targedNumPanels(self):
        return self
    @targedNumPanels.setter
    def targedNumPanels(self,numPanels):
        self._targedNumPanels = numPanels
        if type(self._targedNumPanels) is not int:
            raise Exception('targedNumPanels must be an integer')
            
    @property 
    def numBodies(self):
        return self
    @numBodies.setter
    def numBodies(self,numBodies):
        self._numBodies = numBodies
        if type(self._numBodies) is not int:
            raise Exception('numBodies must be an integer')    
        if self._numBodies != 1:
            raise Exception('numBodies must be 1')  
    
    def readGDF(self):
        # Read GDF file
        with open(self.files['GDF'],'r') as fid:
            lines = fid.readlines()
            
        self.gdfLines = lines
        self.uLen = int(lines[1].split()[0])  
        self.gravity = float(lines[1].split()[1])
        self.isx = float(lines[2].split()[0])
        self.isy = float(lines[2].split()[1])
        self.numFaces = int(lines[3])
        self.numCords = self.numFaces * 4
        self.cords = np.array([temp.split() for temp in lines[4:]]).astype(np.float)

        self.cordsString = [str(temp).replace(",",'').replace('\r','') for temp in lines[4:]] # Output string for Nemoh mesh fil
        for panelNum,i in enumerate(np.arange(4,4+self.numCords,4)):
            self.faces.append(np.array([i-4,i-3,i-2,i-1]))
             
    def readSTL(self):
        reader = vtk.vtkSTLReader()
        reader.SetFileName(self.stlFile)
        reader.Update()
        self.numFaces = int(reader.GetOutput().GetNumberOfCells())
        self.numCords = self.numFaces * 3
        for i in range(self.numFaces):
            n = i*3
            self.faces.append(np.array([n,n+1,n+2,n+2]))
            self.cords.append(np.array(vtk_to_numpy(reader.GetOutput().GetCell(i).GetPoints().GetData())))
        self.cords = np.array(self.cords).reshape([self.numFaces*3,3])
        
    def readVTP(self):
        '''
        Currently assumes that all elements have 4 faces
        '''
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(self.vtpFile)
        reader.Update()
        self.numFaces = int(reader.GetOutput().GetNumberOfCells())
        self.numCords = int(reader.GetOutput().GetNumberOfPoints())
        for i in range(self.numFaces):
            n = i*4
            self.cords.append(np.array(vtk_to_numpy(reader.GetOutput().GetCell(i).GetPoints().GetData())))
            self.faces.append(np.array([n,n+1,n+2,n+3]))
        self.cords = np.array(self.cords).reshape([self.numCords,3])
        
    def writeNemohMesh(self):   
        with open(self.files['Mesh'],'w') as fid:
            fid.write(str(self.numCords))
            fid.write('\n')
            fid.write(str(self.numFaces))
            fid.write('\n')
            for i in range(np.shape(self.cords)[0]):
                fid.writelines(str(self.cords[i]).replace('[','').replace(']',''))
                fid.write('\n')
            
            for i in range(self.numFaces):
                fid.write(str(self.faces[i]).replace('[','').replace(']',''))
                fid.write('\n')
                
    def writeVtp(self):
        '''
        Function to write VTK file for visiuilization in Paravies
        '''
        mesh    = vtk.vtkPolyData()
        points  = vtk.vtkPoints()
        polys   = vtk.vtkCellArray()
        
        for i in range(self.numCords):
            points.InsertPoint(i, self.cords[i])
        for i in range(self.numFaces):
            polys.InsertNextCell( mkVtkIdList(self.faces[i]) )
            
        mesh.SetPoints(points)
        mesh.SetPolys(polys)
        
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(self.files['vtk'])
        writer.SetInputData(mesh)
        writer.SetDataModeToAscii()
        writer.Write()
        
    def writeMeshCal(self):
        with open(self.files['Mesh.cal'],'w') as fid:
            fid.write(self.files['Mesh-noPath'])
            fid.write('\n')
            fid.write('0')
            fid.write('\n')
            fid.write('0. 0.') # not sure what this is
            fid.write('\n')
            fid.write(str(self._cg).replace('[','').replace(']','').replace(',',''))
            fid.write('\n')
            fid.write(str(self.targetNumPanels))
            fid.write('\n')
            fid.write('2')
            fid.write('\n')
            fid.write('0')
            fid.write('\n')
            fid.write('1')
            
    def run(self):
        os.chdir(self.baseDir)
        if os.sys.platform == 'darwin':
            os.system(self.meshExec + '>' + self.files['MeshLog'])
        else:
            raise Exception('Only osx is supported at this time (linux and windows will prob work too with small change to this python module)')
        
class NemohSimulation(object):
    def __init__(self,sim):
        self.baseDir = sim.dir + os.path.sep + '..'
        
        self.files = {}
        self.files['input.txt'] = sim.dir + os.path.sep + 'input.txt'
#        if os.path.exists(self.dir) is False:
#            os.mkdir(self.dir)
#        
#        self.files = {}
#        self.files['ID'] = self.dir + 'ID.dat'
#        self.files['Nemoh'] = self.dir + 'Nemoh.cal'
#        self.files['Normalvelocities'] =  self.dir + 'Normalvelocities.dat'
#        self.files['input'] = self.dir + 'input.txt'
    def writeInput(self):
        with open(self.files['input.txt'],'w') as fid:
            fid.write('\n0')
        
class NemohResults(object):
    def __init__(self,outputDir,plotData=False):
        self.w = []
        self.addedMass = {}
        self.radiationDamping = {}
        self.outputDir = outputDir
        
        self.readCMCA()
        
        if plotData is True:
            self.plotAddedMassAndDamping()
        
    def readCMCA(self):
        '''
        Function to read CM.dat created by Nemoh
        '''        
        # load  data files
        with open(self.outputDir + '/results/CM.dat') as fid :
            linesCM = fid.readlines()        
        with open(self.outputDir + '/results/CA.dat') as fid :
            linesCA = fid.readlines()

        # Read the number of frequencies
        self.numFreqs = int(linesCM[0].split()[-1])

        # Read the Frequencies, the added mass matrix, and the radiation damping matrix at each frequency
        for i in xrange(self.numFreqs):
            self.w.append(float(linesCM[1+i*7].replace('\n','')))
            self.addedMass[self.w[i]] = [temp.replace('\n','') for temp in linesCM[2+i*7:8+i*7]]
            self.addedMass[self.w[i]] = np.array([np.fromstring(temp,sep=' ') for temp in self.addedMass[self.w[i]]])        
            self.radiationDamping[self.w[i]] = [temp.replace('\n','') for temp in linesCA[2+i*7:8+i*7]]
            self.radiationDamping[self.w[i]] = np.array([np.fromstring(temp,sep=' ') for temp in self.radiationDamping[self.w[i]]])
            
        self.addedMassDiag = np.array([np.diag(temp) for temp in self.addedMass.values()])
        self.radiationDampingDiag = np.array([np.diag(temp) for temp in self.radiationDamping.values()])
            
    def plotAddedMassAndDamping(self):
        plt.figure('AddedMass')
        plt.figure('RadiationDamping')        
        
        for i in xrange(6):
            plt.figure('AddedMass')
            plt.plot(self.w,self.addedMassDiag[:,i],label='Component (' + str(i) + ', ' + str(i) + ')')
            plt.title('Diagional Compinent of Added Mass Matrix')
            plt.xlabel('Wave Frequency (rad)')
            plt.ylabel('Added Mass (kg)')
            plt.legend()
            
            plt.figure('RadiationDamping') 
            plt.plot(self.w,self.radiationDampingDiag[:,i],label='Component (' + str(i) + ', ' + str(i) + ')')
            plt.title('Diagional Compinent of Radiation Damping Matrix')
            plt.xlabel('Wave Frequency (rad)')
            plt.ylabel('Radiation Damping')
            plt.legend()
            
        plt.show()

#output = NemohOutput(outputDir='/Users/mlawson/Applications/nemoh/matlabRoutines/simulationTest',plotData=True)


        
