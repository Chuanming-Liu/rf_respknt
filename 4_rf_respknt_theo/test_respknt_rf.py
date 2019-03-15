#!/usr/bin/env python
import respknt
import numpy as np
import sys
import Plot_Module_RF as plot
import matplotlib.pyplot as plt
import timeit
import deconvolution as decv


# model setting
maxlay                       = 50
vp_in                        = np.zeros(maxlay, dtype=np.float32)
vs_in                        = np.zeros(maxlay, dtype=np.float32)
rho_in                       = np.zeros(maxlay, dtype=np.float32)
thick_in                     = np.zeros(maxlay, dtype=np.float32)
# parameter set up
slowness                     = 0.06
dt                           = 0.02 # limit .lt. 0.01
tduring                      = 30.
nl                           = 2

vp_in[:nl]                   = np.array([6.000, 8.0000], dtype = np.float32)
vs_in[:nl]                   = np.array([3.4641, 4.6188], dtype = np.float32)
rho_in[:nl]                  = np.array([2.6900, 3.3300], dtype = np.float32)
thick_in[:nl]                = np.array([30.000, 0.0000], dtype = np.float32)

start 					     = timeit.default_timer()
wavez, waver                 = respknt.respknt_interface(slowness=slowness, dt_in=dt, tduring=tduring, nl=nl, rho_in=rho_in, thick_in=thick_in, vp_in=vp_in, vs_in=vs_in)
time                         = np.arange(0, tduring+dt, dt)
npts                         = len(time)
Ttz                          = wavez[:npts]
Ttr       					 = waver[:npts]
rfr                          = decv.deconvolve(Ztr=Ttz, Rtr=Ttr, dt=dt, npts=npts, f0=2.5)

# plot seismograms
stop 						 = timeit.default_timer()
print('time: ',stop-start, 's')
plmt                         = plot.plot_wavefrom_dataset()
ax1, ax2, ax3, fig           = plmt.pic_frame(title='respknt')
xmax                         = 20
plmt.plot_wavefrom(waveform=Ttz, time=time, ax=ax1, xmax=xmax, label='Z')
plmt.plot_wavefrom(waveform=Ttr, time=time, ax=ax2, xmax=xmax, label='R')
plmt.plot_wavefrom(waveform=rfr, time=time, ax=ax3, xmax=xmax, label='RF')
plot.plot_legend()
plot.pic_save(fname='syn_wave', outdir='.')
