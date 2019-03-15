# -*- coding: utf-8 -*-
"""
2019-03-14
Verify theo and respknt result
"""
import sys
path='/Users/chuanmingliu/Lance/Project/Polarization/Synthetic/Syn_investigation/common_source/Kernel_respknt_src'
sys.path.insert(0,path)
import numpy as np
import os
import Plot_Module_RF as plot
import vmodel
import Respknt_kenerl

# load model
m0                            = vmodel.model1d()
thickness                     = m0.model_1layer_30km()
# load modules
rfk                           = Respknt_kenerl.rf_kernel(m0)
rfk.init_default_layer()
# set up
t                             = 20
dt                            = 0.1
t0                            = 0.
rft                           = rfk.theo_rf(t=t, dt=dt, t0=t0)
rfr                           = rfk.respknt_rf(t=30, dt=dt, t0=t0)
rfs                           = rfk.raysum_rf(t=30, dt=dt, t0=t0)

# plot
plmt                         = plot.plot_wavefrom_dataset()
ax1, ax2, ax3, fig           = plmt.pic_frame(title='respknt')
xmax                         = 20
plmt.plot_wavefrom(waveform=rfr, time=np.arange(len(rfr))*dt, ax=ax1, xmax=xmax, label='respknt')
plmt.plot_wavefrom(waveform=rft, time=np.arange(len(rft))*dt, ax=ax2, xmax=xmax, label='theo')
plmt.plot_wavefrom(waveform=rfs, time=np.arange(len(rfs))*dt, ax=ax3, xmax=xmax, label='raysum')
plot.pic_save(fname='syn_wave', outdir='.')
