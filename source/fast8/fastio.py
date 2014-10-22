"""
Functions for running FAST from Python and for postprocessing FAST data

Created on Tue Sep 16 16:46:12 2014

Author: mlawson
"""
import os
import getpass
from time import gmtime, strftime
import re
import matplotlib.pyplot as plt
import numpy as np

plt.close('all')
   
class FastSimulation(object):
    
    def __init__(self,fastInputFile='',fastExeFile='',fastInputData='',adInputData='',edInputData=''):
        self.fastInputs = fastInputData
        self.adInputs   = adInputData
        self.edInputs   = edInputData
        self.fastInputFile = fastInputFile
        self.exe = fastExeFile
        self._aerodynBladeNodeData = None
        
#        @property
#            def aerodynBladeNodeData(self):
#                with open(self.fastInputFile[0:-3]+'out') as fid :
#                    lines = fid.readlines()
#                for 
#            return self._aerodynBladeNodeData   
        
        # Error checking
        if not os.path.isfile(self.fastInputFile):       
            raise IOError('The fast input file does not exist')        
        if not os.path.isfile(self.exe):        
            raise IOError('The executable file does not exist')
        
        # Overwrite variables in .fst file with user inputs     
        for key,value in self.fastInputs.iteritems():
            setFastInputs(self.fastInputFile,key,value)
        
        # Get AD and ED filenames from fast input file
        try:
            self.adInputFile = self.fastInputs['AeroFile'][1:-1]
            self.edInputFile = self.fastInputs['EDFile'][1:-1]
        except IOError:
            print 'The Aerodyn and ElastoDyn file names must be defined in the fastInputData variable'

        # Error checking        
        if not os.path.isfile(self.adInputFile):       
            raise IOError('The AD input file does not exist')        
        if not os.path.isfile(self.edInputFile):        
            raise IOError('The ED input does not exist')

        # Overwrite variables in AD and ED file with user inputs       
        for key,value in self.adInputs.iteritems():
            setFastInputs(self.adInputFile,key,value)
            
        for key,value in self.edInputs.iteritems():
            setFastInputs(self.edInputFile,key,value)
            
    def run(self):
        """
        Function to run FAST
        """

        os.system(self.exe + ' ' + self.fastInputFile)
    
    def loadOutputs(self,linecolor='k',linestyle=''):
        """
        Function to load FAST .out files
        """        
        
        with open(self.fastInputFile[0:-3]+'out') as fid :
            lines = fid.readlines()
        self.output = FastOutput()
        self.output.linecolor=linecolor
        self.output.linestyle=linestyle
        names = lines[6].split()
        units = lines[7].replace("(", '').replace(")", '').split()
        data = np.array([np.fromstring(lines[i], sep=" ") for i in range(8, len(lines))])
        self.output.units = {}
        for idx,nm in enumerate(names):
            setattr(self.output,nm,data[:,idx])
            self.output.units[nm]=units[idx]
        return self.output
        
#    def plotLocalOutputs(self,varsToPlot='RotPwr'):
#        plotOutputs(self.output,varsToPlot)
        
class FastOutput(object):
    """
    FAST output data class.
    """
    def __init__(self):
        pass

class AerodynBladeNodes(object):
    """
    Aerodyn blade nodes
    """
    def __init__(self):
        pass

def readOutFile(filename):
    """
    Function to load data from a fast .out file.
    
    Inputs:
        filename: name of the the fast .out file (string) \n
        color: linecolor for plotting (string) \n
        linestyle: linetype for plottying (string)
        
    Output:
        color \n
        linestyle \n
        units \n
        data \n
    
    Adapted from lkilcher
    """
    
    with open(filename) as fid :
        lines = fid.readlines()
    out = FastOutput()
    names = lines[6].split()
    units = lines[7].replace("(", '').replace(")", '').split()
    data = np.array([np.fromstring(lines[i], sep=" ") for i in range(8, len(lines))])
    out.units = {}
    for idx,nm in enumerate(names):
        setattr(out,nm,data[:,idx])
        out.units[nm]=units[idx]
    return out

def setFastInputs(fName,varName,newValue):
    """
    This function modifies the value of FAST user inputs in FAST, AeroDyn, 
    and ElastoDyn input files.
    
    Inputs:
        fname: name of file (string) \n
        varName: Name of variable to set (string) \n
        newValue: Value of variable to be set (string) \n
    
    Outputs:
        none
    """

    with open(fName, 'r') as inFile, open('temp', 'w') as outFile:
        for lineNum, line in enumerate(inFile):
            #if varName in line:
            if re.search(r"\b" + re.escape(varName) + r"\b", line):
                outFile.write(newValue + '\t' + varName + '\t- modified by ' + getpass.getuser() + ' ' + strftime("%Y-%m-%d %H:%M:%S", gmtime()) + '\n')
                print(fName + ', line ' + str(lineNum) + '\t: '+ varName + ' = ' + newValue)
                check = {}
            else:
                outFile.write(line)
    try:
        check
    except:
        raise IOError('Error: The variable '+varName+' was not found in the file ' + fName)
    os.remove(fName)
    os.rename('temp',fName)
  
def plotOutputs(data,varNames):
    """
    Function to plot fast output data.
    
    Inputs:
        data: List of fast input files. E.g. ['run1.out','run2.out','run3.out'] \n
        varNames: Variable names to plot. E.g. ['RotPwr','RotThrust'] \n
    """

    fig, axarr = plt.subplots(np.size(varNames), sharex=True)
    for key in data:
        for i,varName in enumerate(varNames):
            axarr[i].plot(data[key].Time,getattr(data[key],varName),data[key].linestyle,color=data[key].linecolor,label=key)	
            axarr[i].set_ylabel(varName) # Note: cannot put units in label due to unicode/Python problem
    axarr[i-1].set_xlabel('Time (s)')
    axarr[0].legend()