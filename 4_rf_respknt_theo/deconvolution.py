# -*- coding: utf-8 -*-
"""
Module for synthetic receiver function computation,
    :Copyright:
    Author: Lili Feng and Chuanming Liu
"""
import copy
import numba
import numpy as np


# @numba.jit(numba.float32[:](numba.float32[:], numba.int32, numba.float32, numba.float32))
def _phaseshift( x, nfft, DT, TSHIFT ):
    """Add a shift to the data into the freq domain, private function for _iter_deconvolve
    """
    Xf      = np.fft.fft(x)
    # phase shift in radians
    shift_i = round(TSHIFT/DT) # removed +1 from here.
    p       = np.arange(nfft)+1
    p       = 2*np.pi*shift_i/(nfft)*p
    # apply shift
    Xf      = Xf*(np.cos(p) - 1j*np.sin(p))
    # back into time
    x       = np.real( np.fft.ifft(Xf) )/np.cos(2*np.pi*shift_i/nfft)
    return x

# @numba.jit(numba.float32[:](numba.float32[:], numba.float32[:], numba.float32))
def _FreFilter(inW, FilterW, dt ):
    """Filter input array in frequency domain, private function for _iter_deconvolve
    """
    FinW    = np.fft.fft(inW)
    FinW    = FinW*FilterW*dt
    FilterdW= np.real(np.fft.ifft(FinW))
    return FilterdW

# @numba.jit(numba.float32[:](numba.float32, numba.int32, numba.float32))
def _gaussFilter( dt, nft, f0 ):
    """
    Compute a gaussian filter in the freq domain which is unit area in time domain
    private function for _iter_deconvolve
    ================================================================================
    Input:
    dt  - sampling time interval
    nft - number freq points
    f0  - width of filter

    Output:
    gauss  - Gaussian filter array (numpy)
    filter has the form: exp( - (0.5*w/f0)^2 ) the units of the filter are 1/s
    ================================================================================
    """
    df      = 1.0/(nft*dt)
    nft21   = int(0.5*nft + 1)
    # get frequencies
    f       = df*np.arange(nft21, dtype=np.float32)
    w       = 2*np.pi*f
    w       = w/f0
    kernel  = w**2
    # compute the gaussian filter
    gauss   = np.zeros(nft, dtype=np.float32)
    gauss[:nft21]   = np.exp( -0.25*kernel )/dt
    gauss[nft21:]   = np.flipud(gauss[1:nft21-1])
    return gauss

# @numba.jit(numba.float32[:](numba.float32[:], numba.float32[:], numba.float32, numba.int32, numba.int32, numba.float32, numba.float32, numba.float32))
def _iter_deconvolve(Ztr, RTtr, dt, npts, niter, tdel, f0, minderr):
    """
    Iterative deconvolution
    """
    RMS         = np.zeros(niter, dtype = np.float32)  # RMS errors
    nfft        = 2**(npts-1).bit_length() # number points in fourier transform
    P0          = np.zeros(nfft, dtype = np.float32)# predicted spikes
    # Resize and rename the numerator and denominator
    U0          = np.zeros(nfft, dtype = np.float32) #add zeros to the end
    W0          = np.zeros(nfft, dtype = np.float32)
    U0[:npts]   = RTtr  # clear UIN;
    W0[:npts]   = Ztr   # clear WIN;
    # get filter in Freq domain
    gauss       = _gaussFilter( dt, nfft, f0 )
    # filter signals
    Wf0         = np.fft.fft(W0)
    FilteredU0  = _FreFilter(U0, gauss, dt )
    FilteredW0  = _FreFilter(W0, gauss, dt )
    R           = FilteredU0 #  residual numerator
    # Get power in numerator for error scaling
    powerU      = np.sum(FilteredU0**2)
    # Loop through iterations
    it          = 0
    sumsq_i     = 1
    d_error     = 100*powerU + minderr
    maxlag      = int(0.5*nfft)
    while( abs(d_error) > minderr  and  it < niter ):
        it          = it+1 # iteration advance
        #   Ligorria and Ammon method
        RW          = np.real(np.fft.ifft(np.fft.fft(R)*np.conj(np.fft.fft(FilteredW0))))
        sumW0       = np.sum(FilteredW0**2)
        RW          = RW/sumW0
        imax        = np.argmax(abs(RW[:maxlag]))
        amp         = RW[imax]/dt; # scale the max and get correct sign
        #   compute predicted deconvolution
        P0[imax]    = P0[imax] + amp  # get spike signal - predicted RF
        P           = _FreFilter(P0, gauss*Wf0, dt*dt ) # convolve with filter
        #   compute residual with filtered numerator
        R           = FilteredU0 - P
        sumsq       = np.sum(R**2)/powerU
        RMS[it-1]   = sumsq # scaled error
        d_error     = 100*(sumsq_i - sumsq)  # change in error
        sumsq_i     = sumsq  # store rms for computing difference in next
    # Compute final receiver function
    P   = _FreFilter(P0, gauss, dt )
    # Phase shift
    P   = _phaseshift(P, nfft, dt, tdel)
    # output first nt samples
    RFI = P[:npts]
    # output the rms values
    RMS = RMS[:it]
    if it > 1:
        return RFI, (1.0-RMS[it-1])*100.0
    else:
        return RFI, (1.-d_error)*100.



def deconvolve(Ztr, Rtr, dt, npts, tdel=0., f0=2.5, niter=200, minderr=0.0001):
    """
    Compute receiver function from synthetics waveforms (respknt/raysum)
	with iterative deconvolution algorithmn
    ========================================================================================================================
    ::: input parameters :::
	dt          - sampling interval (s)
	npts        - points of waveform
    tdel        - phase delay
    f0          - Gaussian width factor
    niter       - number of maximum iteration
    minderr     - minimum misfit improvement, iteration will stop if improvement between two steps is smaller than minderr
    ::: output :::
    rfr  		- trace for radial receiver functions (same time scale with input rf)
    ========================================================================================================================
    """
    rfr, fitness= _iter_deconvolve(Ztr, Rtr, dt, npts, niter, tdel, f0, minderr)
    if fitness < 95.:
        print 'WARNING: fittness is',fitness
    return rfr
