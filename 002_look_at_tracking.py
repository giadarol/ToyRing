import metaclass as mtc
import numpy as np

from scipy.constants import c as c_light

import matplotlib.pyplot as pl
pl.close('all')

for ii in range(10):
	fname = 'track.obs0001.p%04d'%(ii+1)
	ob = mtc.twiss(fname)



	pl.figure(1)

	pl.plot(ob.X)

	pl.figure(2)
	axl1 = pl.subplot(2,1,1)
	pl.plot(ob.T)
	axl2 = pl.subplot(2,1,2)
	pl.plot(ob.PT)

	pl.figure(3)
	pl.plot(ob.T/c_light*1e9, ob.PT)
	pl.grid('on')


pl.show()