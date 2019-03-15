import theo
import numpy as np
import pdb
# set up
maxlay                       = 50
vp_in                        = np.zeros(maxlay, dtype=np.float32)
vs_in                        = np.zeros(maxlay, dtype=np.float32)
thick_in                     = np.zeros(maxlay, dtype=np.float32)
nl                           = 2
vp_in[:nl]                   = np.array([6.000, 8.0000], dtype = np.float32)
vs_in[:nl]                   = np.array([3.4641, 4.6188], dtype = np.float32)
thick_in[:nl]                = np.array([30.000, 0.0000], dtype = np.float32)
# RF from theo
maxsublayer                  = 1000
hin                          = np.zeros(maxsublayer, dtype=np.float64)
vsin                         = np.zeros(maxsublayer, dtype=np.float64)
vpvs                         = np.zeros(maxsublayer, dtype=np.float64)
qsin                         = 600.*np.ones(maxsublayer, dtype=np.float64)
qpin                         = 1400.*np.ones(maxsublayer, dtype=np.float64)
# input
hin[:nl]                     = thick_in[:nl]
vsin[:nl]                    = vs_in[:nl]
vpvs[:nl]                    = vp_in[:nl]/vs_in[:nl]
# incident angle
dt                           = 0.1 # limit .lt. 0.01
fs                           = 1./dt
slowness                     = 0.06
din                          = 180. * np.arcsin(vsin[nl-1]*vpvs[nl-1]*slowness) / np.pi
din                          = 28.
gaus                         = 2.5
minAmp                       = 0.05
t0                           = 0.
rt                           = 0.1
nn                           = 15./rt
rx0, rxhv, rzhv              = theo.theo(nl, vsin, hin, vpvs, qpin, qsin, rt, din, 2.5, 0.005, 0, nn)
print rx0
