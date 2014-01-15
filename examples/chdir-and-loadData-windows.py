import os
import numpy, scipy, matplotlib

os.chdir('F:/Dropbox/currentWork/pythonExamples')
a = numpy.loadtxt('testData/formattedData.dat')
b = a[:,0]
print type(a)
print a
print b
print b.transpose()
