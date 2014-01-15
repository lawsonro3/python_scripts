"""
Created on Thu Mar 29 12:01:28 2012

@author: mlawson

This script need significant work - i'm not even sure if its right
"""
from scipy import * 
from pylab import * 

npts = 512            #Use some power of 2
t=linspace(0,1,npts+1)     # Use 2^N + 1 
dt = (t[-1]-t[0])/(len(t)-1)    # Maximum frequency is 1/2dt ?
fmax = 1/(2*dt) 
f1 = 80 
f2 = 90 
#sig = 1 + sin(2*pi*f1*t) + 1 + sin(2*pi*f2*t)  # sum of signals 
sig=(1+sin(2*pi*f1*t))*(1+sin(2*pi*f2*t))      # product of signals

figure(1)
plot(t,sig);xlabel('Time');title('Signal')

ft = fft(sig,n=npts) 
mgft=abs(ft)             #Get magnitude of fft
df = fmax/float(npts/2)
f=linspace(0,fmax,npts/2+1)
print 'fmax = ',fmax,' df = ',df,' ','\n 1st freqs = ',f[0:5]

figure(2)
plot(f,mgft[0:npts/2+1]);title('Fast Fourier Transform Magnitude')
xlabel('frequency')
ylabel('fft magnitude')

figure()
plot(linspace(0,size(mgft)-1,size(mgft)),mgft)
show()