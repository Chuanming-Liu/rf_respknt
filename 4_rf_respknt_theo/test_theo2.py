import theo
import numpy as np
import pdb
# set up


maxsublayer           = 1000
vsin                  = np.zeros(maxsublayer, dtype=np.float64)
vpvs                  = np.zeros(maxsublayer, dtype=np.float64)
qsin                  = 600.*np.ones(maxsublayer, dtype=np.float64)
qpin                  = 1400.*np.ones(maxsublayer, dtype=np.float64)
tthick                = np.zeros(maxsublayer, dtype=np.float64)

nl                    = 2
vsin[:nl]             = np.array([3.4641, 4.6188], dtype = np.float32)
vpvs[:nl]             = np.array([1.71, 1.71], dtype = np.float32)
tthick[:nl] 	      = np.array([30.000, 0.0000], dtype = np.float32)
rt                    = 0.1
din                   = 28.
nn                    = 15./rt
rx,rxhv,rzhv          = theo.theo(nl, vsin, tthick, vpvs, qpin, qsin, rt, din, 2.5, 0.005, 0, nn)
print rx
