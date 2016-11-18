import numpy as np
import sys
import importlib

mod = importlib.import_module(sys.argv[1])
en=np.linspace(2.0,8.0,110)
ph=np.linspace(0,1,120)
for ee in en:
    ff=mod.energy_spectrum(ee,ph)
    gg=mod.polarization_degree(ee,ph,0,0)*np.cos(2*mod.polarization_angle(ee,ph,0,0))
    gg=np.mean(gg*ff)
    ff=np.mean(ff)
    print('%g %g %g %g' % (ee, ff, gg, gg/ff))
