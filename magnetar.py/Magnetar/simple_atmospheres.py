try:
    from numba import jit
except:
    def jit(**kwargs):
        def jitd(func):
            def helper(*args, **kwargs):
                return func(*args, **kwargs)
            return helper
        return jitd
    
from Magnetar.utils import atmosphere
import numpy as np

@jit(nopython=True,parallel=True)
def bbfunk( ee, tt):  # per mode
    return 208452.792 * ee**3 / np.expm1(ee / tt) / 2


def bb_atmo_f(temp):
    return atmosphere() * (lambda ee: bbfunk(ee,temp))


#
# blackbody atmosphere (convenience class with parallel structure to condensed_surface)
#
class bb_atmo(atmosphere):
    def __init__(self,teff,mag_strength,mag_inclination,*args,ofraction=0.5,**kwargs):
        self.teff  = teff
        self.ofraction = np.clip(ofraction,0,1)
        self.xtemp = (1-self.ofraction)**0.25
        self.otemp = self.ofraction**0.25
        rat=teff*(2.0/(self.xtemp**4+self.otemp**4))**0.25
        self.xtemp*=rat
        self.otemp*=rat
        self.mag_inclination = mag_inclination
        self.mag_strength=mag_strength
    def __str__(self):
        outstring='''#
# class bb_atmo
#
# effective_temperature %12g keV
# O fraction            %12g
# X temperature         %12g keV
# O temperature         %12g keV
# mag_strength          %12g Gauss [not used]
# mag_inclination       %12g radians [not used]
#\n''' % (self.teff, self.ofraction, self.xtemp, self.otemp, self.mag_strength, self.mag_inclination)
        return outstring+atmosphere.__str__(self)

    def xintensity(self, dataarray):
        if self.xtemp>0:
            return bbfunk(dataarray[-1], self.xtemp)
        else:
            return 0*dataarray[-1]

    def ointensity(self, dataarray):
        if self.otemp>0:
            return bbfunk(dataarray[-1], self.otemp)
        else:
            return 0*dataarray[-1]

#
# pure X blackbody atmosphere (convenience function with parallel structure to condensed_surface)
#
def bb_atmo_purex(teff,mag_strength,mag_inclination,*args,**kwargs):
    return bb_atmo(teff,mag_strength,mag_inclination,ofraction=0)

class Thompson_Kostenko_Magnetosphere(atmosphere):
    def __init__(self,mag_inclination,normalization=1,alpha=1):
        self.mag_inclination=mag_inclination
        self.normalization=normalization
        self.alpha=alpha
    def __str__(self):
        outstring='''#
# class modified_bb_atmo
#
# normalization %12g at 10 keV
# alpha         %12g
''' % (normalization,alpha)
        return outstring+atmosphere.__str__(self)
    @jit(nopython=True,parallel=True)
    def _ointensity(dataarray,mag_inclination,normalization,alpha):
        ee=dataarray[-1]/10.0
        sigmao=np.abs(surface_temperature/ee)**freq_power
        coskb2=(np.cos(np.radians(mag_inclination)
                           ) * np.cos(np.radians(dataarray[-3])) +
                    np.sin(np.radians(mag_inclination)
                           ) * np.sin(np.radians(dataarray[-3])
                                  ) * np.cos(np.radians(dataarray[-2])))**2
        return normalization*ee**alpha*(1.0-coskb2)
    def ointensity(self, dataarray):
        return Thompson_Kostenko_Magnetosphere._ointensity(np.array(dataarray),self.mag_inclination,self.normalization,self.alpha)
    def xintensity(self, dataarray):
        return 0.0 # np.zeros((len(dataarray[-1])))
      
        

#
# modified blackbody atmosphere (convenience class with parallel structure to condensed_surface)
#
# freq_power quantifies how the opacity depends on photon energy: kappa goes as (1/energy)**freq_power
# freq_power is 2 (E<T) to 3 (E>t) for free-free opacity and zero for scattering
#
# sigma_power quantifies how the temperature depends on column density: T goes as sigma**sigma_power
# sigma_power is (alpha+1)/(4+alpha-beta) if opacity is proportional to density**alpha temperature**beta
#
#    unmagnetized free-free : alpha=1, beta=-3.5 -> 0.23529411764705882 (4/17)
#    magnetized free-free   : alpha=1, beta=-1.5 -> 0.3076923077 (4/13)         
#    electron-scattering    : alpha=0, beta=0    -> 0.25 (4/16)
#
# Based on the power-law atmospheres for neutron stars in
#
# https://ui.adsabs.harvard.edu/abs/1998MNRAS.300..599H for the magnetized case and
#
# https://ui.adsabs.harvard.edu/abs/1984ApJ...287..244H for unmagnetized case
#
#

class modified_bb_atmo(atmosphere):
    def __init__(self,teff,mag_strength,mag_inclination,*args,freq_power=2,sigma_power=4./13.,kb_suppress=True,limb_darkening=True,**kwargs):
        self.effective_temperature  = teff
        self.freq_power = freq_power
        self.sigma_power = sigma_power
        self.surface_temperature = teff
        self.mag_inclination = mag_inclination
        self.mag_strength=mag_strength
        self.ecyc = mag_strength/4.4e13*511.
        self.kb_suppress=kb_suppress
        self.limb_darkening=limb_darkening
        self.adjust_surface_temperature()
        
    def __str__(self):
        outstring='''#
# class modified_bb_atmo
#
# effective_temperature %12g keV
# surface_temperature   %12g keV [tau=1 for upgoing O-mode with E=surface_temperature, neglecting k dot b suppression]
# mag_strength          %12g Gauss
# mag_inclination       %12g radians
# cyclotron energy      %12g keV
# k dot b supression    %12s 
# limb darkening        %12s 
# freq_power            %12g [cross-section goes as 1/freq**freqpower for O-mode]
# sigma_power           %12g [temperature goes as column-density**sigma_power]
#       
''' % (self.effective_temperature, self.surface_temperature, self.mag_strength, self.mag_inclination, self.ecyc, self.kb_suppress, self.limb_darkening, self.freq_power, self.sigma_power )
        return outstring+atmosphere.__str__(self)
    
    # Auxiliary functions defined for numba
    @jit(nopython=True,parallel=True)
    def _xintensity(dataarray,surface_temperature,freq_power,sigma_power,ecyc,limb_darkening):
        ee=dataarray[-1]
        sigmax=np.abs(surface_temperature/ee)**freq_power
        sigmax*=np.where(ee<ecyc,(ee/ecyc)**2,1)
        if limb_darkening:
            tempx=surface_temperature*np.abs(np.cos(np.radians(dataarray[-3]))/sigmax)**sigma_power
        else:
            tempx=surface_temperature*np.abs(1/sigmax)**sigma_power
        return bbfunk(ee,tempx)

    @jit(nopython=True,parallel=True)
    def _ointensity(dataarray,surface_temperature,freq_power,sigma_power,ecyc,limb_darkening,kb_suppress,mag_inclination):
        ee=dataarray[-1]
        sigmao=np.abs(surface_temperature/ee)**freq_power
        if kb_suppress:
            coskb2=(np.cos(np.radians(mag_inclination)
                       ) * np.cos(np.radians(dataarray[-3])) +
                np.sin(np.radians(mag_inclination)
                       ) * np.sin(np.radians(dataarray[-3])
                                  ) * np.cos(np.radians(dataarray[-2])))**2
            sigmao*=np.where(ee<ecyc,(1-coskb2)+coskb2*(ee/ecyc)**2,1)
        if limb_darkening:
             tempo=surface_temperature*np.abs(np.cos(np.radians(dataarray[-3]))/sigmao)**sigma_power
        else:
             tempo=surface_temperature*np.abs(1/sigmao)**sigma_power
        return bbfunk(ee,tempo)


    
    def xintensity(self, dataarray):
        return modified_bb_atmo._xintensity(np.array(dataarray),self.surface_temperature,self.freq_power,self.sigma_power,self.ecyc,self.limb_darkening)
 
    def ointensity(self, dataarray):
        return modified_bb_atmo._ointensity(np.array(dataarray),self.surface_temperature,self.freq_power,self.sigma_power,
                                            self.ecyc,self.limb_darkening,self.kb_suppress,self.mag_inclination)
