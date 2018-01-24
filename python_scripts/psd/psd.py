import numpy as np
from matplotlib.mlab import detrend_linear as detrend
fft = np.fft.fft


def psd_freq(nfft, fs):
    """
    Compute the frequency spectrum for use with psd.
    *fs* is the sampling frequency.
    """
    fs = np.float64(fs)
    return np.arange(fs / nfft, fs / 2 + fs / nfft, fs / nfft)


def cohere(a, b, nfft, window='hann', debias=True, noise=(0, 0)):
    """
    Computes the magnitude-squared coherence of *a* and *b*:

                    |S_{ab}|^2
         C_{ab} = --------------
                   S_aa * S_bb

    Here S_xy, S_xx and S_yy are the cross, and auto spectral
    densities of the signal a and b.

    *debias* specify whether to debias the signal (Benignus1969).

    The *noise* keyword may be used to specify the signals'
    noise levels (std of noise in a,b). If *noise* is a two
    element tuple or list, the first and second elements specify
    the noise levels of a and b, respectively.
    default: noise=(0,0)

    """
    nens = np.fix(2. * min(len(a), len(b)) / nfft)
    if noise.__class__ not in [list, tuple, np.ndarray]:
        noise = [noise, noise]
    elif len(noise) == 1:
        noise = [noise[0], noise[0]]
    if nens < 10:
        print "Warning: Only %d ensembles in average. " % (nens)
    if nens <= 2:
        raise Exception("Coherence must be computed from a set of ensembles.")
    # fs=1 is ok because it comes out in the normalization.  (noise
    # normalization depends on this)
    out = (cpsd(a, b, nfft, 1, window=window)) ** 2 / \
        ((psd(a, nfft, 1, window=window) - noise[0] ** 2 / np.pi) * (
         psd(b, nfft, 1, window=window) - noise[1] ** 2 / np.pi))
    if debias:
        return out * (1 + 1. / nens) - 1. / nens  # This is from Benignus1969, it seems to work (make data with different nens (nfft) agree).
    return out


def cpsd(a, b, nfft, fs, window='hann'):
    """
    Compute the cross power spectral density (CPSD) of the signals *a* and *b*.

    This performs:
    fft(a)*conj(fft(b))
    Note that this is consistent with *np.correlate*'s definition of correlation.
    (The conjugate of D.B. Chelton's definition of correlation.)

    The two signals should be the same length, and should both be real.

    See also:
    psd,cohere

    The units of the spectra is the product of the units of *a* and *b*, divided by the units of fs.
    """
    if np.iscomplexobj(a) or np.iscomplexobj(b):
        raise Exception
    auto_psd = False
    if a is b:
        auto_psd = True
    max_ind = len(a)
    nfft = min([nfft, max_ind])
    repeats = np.fix(2. * max_ind / nfft)
    fs = np.float64(fs)
    if max_ind == nfft:
        repeats = 1
    if window == 'hann':
        wind = np.hanning(nfft)
    elif window is None or window == 1:
        wind = np.ones(nfft)
    fft_inds = slice(1, np.floor(nfft / 2. + 1))
    wght = 2. / (wind ** 2).sum()
    s1 = fft(detrend(a[0:nfft]) * wind)[fft_inds]
    if auto_psd:
        pwr = np.abs(s1) ** 2
    else:
        pwr = s1 * np.conj(fft(detrend(b[0:nfft]) * wind)[fft_inds])
    if repeats - 1:
        step = np.fix((max_ind - nfft) / (repeats - 1))
        for i in range(step, max_ind - nfft + 1, step):
            s1 = fft(detrend(a[i:(i + nfft)]) * wind)[fft_inds]
            if auto_psd:
                pwr += np.abs(s1) ** 2
            else:
                pwr += s1 * \
                    np.conj(fft(detrend(b[i:(i + nfft)]) * wind)[fft_inds])
    pwr *= wght / repeats / fs
    if auto_psd:  # No need to take the abs again.
        return pwr
    return np.abs(pwr)


def psd(a, nfft, fs, window='hann'):
    """
    Compute the power spectral density (PSD).

    This function computes the one-dimensional *n*-point PSD.

    The units of the spectra is the product of the units of *a* and *b*, divided by the units of fs.

    Parameters
    ----------
    a      : array_like, currently only supports vectors.
             Input array, can be complex.
    nfft   : int
             Length of the fft to use.  Should be less then len(a).
    fs     : float
             Sampling frequency, for use in normalizing the psd, and
             for computing the frequencies.

    Returns
    -------
    pow    : complex 1d-array
             The PSD, is supposed to normalize to preserve variance.

    fre    : 1d-array
             The frequencies.

    --credit: This was copied from JN's fast_psd.m routine. --

    It detrends the data, uses a hanning window, and a
    minimum of 50% overlap.  Segments are fit so that all of the
    data is used.

    See also:
    numpy.fft

    """
    return cpsd(a, a, nfft, fs, window=window)
