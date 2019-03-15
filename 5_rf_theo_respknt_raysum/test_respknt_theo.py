#!/usr/bin/env python
import respknt
import numpy as np
import sys
import Plot_Module_RF as plot
import matplotlib.pyplot as plt
import timeit
import deconvolution as decv
import theo
import pdb


# model setting
maxlay                       = 50
vp_in                        = np.zeros(maxlay, dtype=np.float32)
vs_in                        = np.zeros(maxlay, dtype=np.float32)
rho_in                       = np.zeros(maxlay, dtype=np.float32)
thick_in                     = np.zeros(maxlay, dtype=np.float32)
# parameter set up
slowness                     = 0.06
dt                           = 0.1 # limit .lt. 0.01
tduring                      = 30.
nl                           = 2

vp_in[:nl]                   = np.array([6.000, 8.0000], dtype = np.float32)
vs_in[:nl]                   = np.array([3.4641, 4.6188], dtype = np.float32)
rho_in[:nl]                  = np.array([2.6900, 3.3300], dtype = np.float32)
thick_in[:nl]                = np.array([30.000, 0.0000], dtype = np.float32)

# RF from respknt
start 					     = timeit.default_timer()
wavez, waver                 = respknt.respknt_interface(slowness=slowness, dt_in=dt, tduring=tduring, nl=nl, rho_in=rho_in, thick_in=thick_in, vp_in=vp_in, vs_in=vs_in)
time                         = np.arange(0, tduring+dt, dt)
npts                         = len(time)
Ttz                          = wavez[:npts]
Ttr       					 = waver[:npts]
Rf                           = decv.deconvolve(Ztr=Ttz, Rtr=Ttr, dt=dt, npts=npts, f0=2.5)
maxData                      = 1024
rfr                          = Rf[:maxData]
stop 						 = timeit.default_timer()
print('respknt time: ',stop-start, 's')
print rfr
#        ---------------------------              #
# RF from theo
maxsublayer                  = 1000
vsin                         = np.zeros(maxsublayer, dtype=np.float64)
vpvs                         = np.zeros(maxsublayer, dtype=np.float64)
qsin                         = 600.*np.ones(maxsublayer, dtype=np.float64)
qpin                         = 1400.*np.ones(maxsublayer, dtype=np.float64)
tthick                       = np.zeros(maxsublayer, dtype=np.float64)

nl                           = 3
vp_in[:nl]                   = np.array([6.000, 8.0000, 8.0000], dtype = np.float32)
vs_in[:nl]                   = np.array([3.64, 4.57, 4.6188], dtype = np.float32)
thick_in[:nl]                = np.array([30.000, 0.0000, 0.0000], dtype = np.float32)

# input
tthick[:nl]                  = thick_in[:nl]
vsin[:nl]                    = vs_in[:nl]
vpvs[:nl]                    = vp_in[:nl]/vs_in[:nl]
# incident angle
din                          = 180. * np.arcsin(vsin[nl-1]*vpvs[nl-1]*slowness) / np.pi
gaus                         = 2.5
minAmp                       = 0.005
t0                           = 5.
rt                           = dt
nn                           = maxData
#  maxdata = 1024
start 					     = timeit.default_timer()
rx0, rxhv, rzhv              = theo.theo(nl, vsin, tthick, vpvs, qpin, qsin, rt, din, gaus, minAmp, t0, nn)
stop 					     = timeit.default_timer()
print('theo time: ',stop-start, 's')
print rx0


# plot seismograms

plmt                         = plot.plot_wavefrom_dataset()
ax1, ax2, ax3, fig           = plmt.pic_frame(title='respknt')
xmax                         = 20
plmt.plot_wavefrom(waveform=rfr, time=np.arange(len(rfr))*dt, ax=ax1, xmax=xmax, label='respknt')
plmt.plot_wavefrom(waveform=rx0, time=np.arange(len(rx0))*dt, ax=ax2, xmax=xmax, label='theo')
# plmt.plot_wavefrom(waveform=rfr, time=time, ax=ax3, xmax=xmax, label='RF')
plot.plot_legend()
plot.pic_save(fname='syn_wave', outdir='.')
