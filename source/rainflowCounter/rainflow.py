"""
Created on Fri Apr  6 10:37:58 2012

@author: Michael Lawson

Simple rainflow counting algorythim

See "Simple rainflow counting algorithms" by Bowning and Socie, 1982
This code uses method 1 from the paper, which requires the entire signal
"""
from pylab import *
interactive(True)

class CreateData(object):
    '''
    The data object used by the rainflow counting functions defined in this program
    '''
    rawData = array([])
    cycles = array([])
    cyclesInd = list()
    peaksValueSorted = array([])
    peaksValue = list()
    peaksInd = list()    
    peaksType = list()
    
    # Load the data from a vector or a text file
    def __init__(self,dataInput=None):
        # Check the data and decide what to do
        if type(dataInput) == ndarray:
            self.rawData = dataInput    
        elif type(dataInput) == str:
            self.rawData = loadtxt(dataInput)
        else:
            raise 'Must provide a vector of data or a data file with a vector of data'
    
def Count(rf,plotData=True):
    '''
    A very memory intensive rainflow counting alrigythim. Needs to be improved
    '''    
    # Find the maximum peak and rearrange the data
    rf = FindPeaks(rf,plotData)
    maxPeakVal = max(rf.peaksValue)
    maxPeakInd = argmax(rf.peaksValue)
    rf.peaksValueSorted = append(rf.peaksValue[maxPeakInd:],rf.peaksValue[0:maxPeakInd+1])
    nPoints = size(rf.peaksValueSorted)
  
    # Perform the rainflow countt
    Eind = array([])
    for i in range(nPoints):
        if remainder(i,1000) == 0:
            print 'Processing data peak = ',i
        Eind = append(Eind,float(i))
#        print 'Eind= ',Eind
        if size(Eind) >= 3:
            X = abs(rf.peaksValueSorted[Eind[-1]]-rf.peaksValueSorted[Eind[-2]])
            Y = abs(rf.peaksValueSorted[Eind[-2]]-rf.peaksValueSorted[Eind[-3]])
            if X >= Y:
                rf.cycles = append(rf.cycles,Y)
                rf.cyclesInd.append([Eind[-2],Eind[-3]])
                Eind = append(Eind[0:-3],Eind[-1])
#                print 'X = ',X
#                print 'Y = ',Y
#                print 'rf.cyclesInd =', rf.cyclesInd
                
    if size(Eind) >= 3:
        X = abs(rf.peaksValueSorted[Eind[-1]]-rf.peaksValueSorted[Eind[-2]])
        Y = abs(rf.peaksValueSorted[Eind[-2]]-rf.peaksValueSorted[Eind[-3]])
        if X >= Y:
            rf.cycles = append(rf.cycles,Y)
            rf.cyclesInd.append([Eind[-2],Eind[-3]])
            Eind = append(Eind[0:-3],Eind[-1])
#            print 'X = ',X
#            print 'Y = ',Y
#            print 'rf.cyclesInd =', rf.cyclesInd
     
    return rf

def Hist(rf,bins):
    figure('Rainflow histogram')
    a = hist(rf.cycles,bins)
        
def FindPeaks(rf,plotData=False):
    '''
    Brute force method of finding peaks and valleys in data
    
    Input: The rainflow data object
    Output: Fills in the 
    '''
    for i in range(len(rf.rawData)):
        if i == 0 and rf.rawData[i] != rf.rawData[i+1]:
            rf.peaksValue.append(rf.rawData[i])
            rf.peaksInd.append(i)
            rf.peaksType.append('start')
        elif i == len(rf.rawData)-1 and rf.rawData[i] != rf.rawData[i-1]:
            rf.peaksValue.append(rf.rawData[i])
            rf.peaksInd.append(i)
            rf.peaksType.append('end')
        elif rf.rawData[i] > rf.rawData[i-1] and rf.rawData[i] > rf.rawData[i+1]:
            rf.peaksValue.append(rf.rawData[i])
            rf.peaksInd.append(i)
            rf.peaksType.append('peak')
        elif rf.rawData[i] < rf.rawData[i-1] and rf.rawData[i] < rf.rawData[i+1]:
            rf.peaksValue.append(rf.rawData[i])
            rf.peaksInd.append(i)
            rf.peaksType.append('valley')         
    rf.peaksValue = array(rf.peaksValue)
    rf.peaksInd = array(rf.peaksInd)
        
    if plotData is True:
        figure('Data peaks')
        plot(rf.rawData,'k')
        plot(rf.peaksInd,rf.peaksValue,'xr')
        grid(True)

    return rf