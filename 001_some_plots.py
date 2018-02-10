import metaclass as mtc
import numpy as np

fname = 'twiss.out'
ob = mtc.twiss(fname)

import matplotlib.pyplot as pl

pl.close('all')
sp1 = pl.subplot(3,1,1)
pl.plot(ob.S, 0.5*np.float_(np.sign(ob.K0L)), '.g')
pl.plot(ob.S, -0.5*np.float_(np.sign(ob.K0L)), '.g')
pl.plot(ob.S, np.float_(np.sign(ob.K1L)), '.b')
sp1.set_ylim(-1.2, 1.2)
pl.subplot(3,1,2,sharex=sp1)
pl.plot(ob.S, ob.BETX, '-b')
pl.plot(ob.S, ob.BETY, '-r')
pl.subplot(3,1,3,sharex=sp1)
pl.plot(ob.S, ob.DX, '.-')

pl.show()