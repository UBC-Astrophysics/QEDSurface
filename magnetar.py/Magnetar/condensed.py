import numpy as np
from Magnetar.utils import atmosphere

class condensed_surface(atmosphere):
    def __init__(self,effective_temperature,mag_strength,mag_inclination,density,fixed_ions=True):
        self.effective_temperature=effective_temperature,mag_strength,mag_inclination
        self.mag_strength=mag_strength
        self.mag_inclination=mag_inclination
        self.blackbody_temperature=self.effective_temperature
        self.fixed_ions=fixed_ions
        self.dens=density
    def _bbfunk(self, ee):  # per mode
        return 208452.792 * ee**3 / np.expm1(ee / self.blackbody_temperature) / 2
    def emissivity_xo(dataarray):
        if (self.fixed_ions):
            return self.emissivity_xo_fixed_ions(dataarray[-3],dataarray[-2],dataarray[-1])
        else:
            return self.emissivity_xo_free_ions(dataarray[-3],dataarray[-2],dataarray[-1])
    def emissivity_xo_fixed_ions(thetak,phik,ene):
        return 1,1
    def emissivity_xo_free_ions(thetak,phik,ene):
        return 1,1
    def calcIQ(self, dataarray):
        ex,eo = self.emissivity_xo(dataarray)
        bbintens=self._bbfunk(dataarray[-1])
        return (eo+ex)*bbintens,(eo-ex)*bbintens
    def xintensity(self,dataarray):
        ex,eo = self.emissivity_xo(dataarray)
        bbintens=self._bbfunk(dataarray[-1])
        return ex*bbintens
    def ointensity(self,dataarray):
        ex,eo = self.emissivity_xo(dataarray)
        bbintens=self._bbfunk(dataarray[-1])
        return eo*bbintens
