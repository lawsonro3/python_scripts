"""
Created on Wed May  2 10:58:47 2012

@author: mlawson

A simple low pass filter for 1-D data
"""
from pylab import *
from scipy.signal import kaiserord, lfilter, firwin, freqz
interactive(True)

def LowPassFilter(x,t=None,sRate=1,db=60,width=.1,cut=.1):
    if t == None:
        t = range(len(x))
        t = array(t)
        
    # The Nyquist rate of the signal.
    nyq_rate = sRate / 2.
 
    # The desired width of the transition from pass to stop,
    # relative to the Nyquist rate.  We'll design the filter
    # with a 5 Hz transition width.
    width = width/nyq_rate
    
    # The desired attenuation in the stop band, in dB.
    ripple_db = db
    
    # Compute the order and Kaiser parameter for the FIR filter.
    N, beta = kaiserord(ripple_db, width)
    
    # The cutoff frequency of the filter.
    cutoff_hz = cut
    
    # Use firwin with a Kaiser window to create a lowpass FIR filter.
    taps = firwin(N, cutoff_hz/nyq_rate, window=('kaiser', beta))
    
    # Use lfilter to filter x with the FIR filter.
    xFilt = lfilter(taps, 1.0, x)
    
    delay = 0.5 * (N-1) / sRate
    print delay
    
    #------------------------------------------------
    # Plot the FIR filter coefficients.
    #------------------------------------------------
    figure('Filter Coefficients (%d taps)' % N)
    plot(taps, 'bo-', linewidth=2)
    grid(True)
    
    #------------------------------------------------
    # Plot the magnitude response of the filter.
    #------------------------------------------------
    figure('Frequency Response')
    clf()
    w, h = freqz(taps, worN=8000)
    plot((w/pi)*nyq_rate, absolute(h), linewidth=2)
    xlabel('Frequency (Hz)')
    ylabel('Gain')
    ylim(-0.05, 1.05)
    grid(True)
    
    
    figure('Comparison of origional and filtered signal')
    
    # Plot the original signal.
    plot(t, x)
    # Plot the filtered signal, shifted to compensate for the phase delay.
    plot(t-delay, xFilt, 'r-')
    # Plot just the "good" part of the filtered signal.  The first N-1
    # samples are "corrupted" by the initial conditions.
    plot(t[N-1:]-delay, xFilt[N-1:], 'g', linewidth=4)
    xlabel('t')
    grid(True)
        
    return xFilt[N-1:], delay