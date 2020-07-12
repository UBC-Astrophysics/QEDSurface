
import numpy as np
from scipy.integrate import simps

def stefan_boltzmann_law(effective_temperature):
    return 208452.792*np.pi**5/15*effective_temperature**4

class atmosphere:
    def __init__(self):
        self.mag_inclination = 0

    def __add__(self, o):
        return add_atmo(self, o)

    def __mul__(self, o):
        if type(o) is tuple or type(o) is list:
            if len(o) > 1:
                return rescale_atmo(self, o[0], o[1])
            else:
                return rescale_atmo(self, o[0], o[0])
        return rescale_atmo(self, o, o)

    def __div__(self, o):
        if type(o) is tuple or type(o) is list:
            if len(o) > 1:
                return apply_beaming(self, o[0], o[1])
            else:
                return apply_beaming(self, o[0], o[0])
        return apply_beaming(self, o, o)

    def __floordiv__(self, o):
        return Complete_Comptonization(self, o)

    def __invert__(self):
        return swap_xo_atmo(self)

    def __pow__(self, y):
        return partial_o_mode_comptonization(self, y)

    def loaddata(self, file, modeltype=None):
        if modeltype is None:
            ext = file.rsplit('/', 1)[-1].rsplit('.', 1)[-1]
            if (ext == 'trj' or ext == 'int' or ext == 'gnu'):
                ss = atmo_lloyd()
            else:
                ss = scatmo_caiazzo()
        else:
            if (modeltype == 'Lloyd'):
                ss = atmo_lloyd()
            elif (modeltype == 'Caiazzo'):
                ss = scatmo_caiazzo()
        return (ss.loaddata(file))

    def xintensity(self, dataarray):
        return 0.5

    def ointensity(self, dataarray):
        return 0.5

    def totalintensity(self, dataarray):
        return self.xintensity(dataarray) + self.ointensity(dataarray)

    def meanxintensity(self, angkbarray):
        return self.xintensity(angkbarray)

    def meanointensity(self, angkbarray):
        return self.ointensity(angkbarray)

    def meantotalintensity(self, angkbarray):
        return self.meanxintensity(angkbarray) + self.meanointensity(angkbarray)

    def calcIQ(self, angarray):
        oo = self.ointensity(angarray)
        xx = self.xintensity(angarray)
        return oo + xx, oo - xx

    def calcmeanIQ(self, angkbarray):
        oo = self.meanointensity(angkbarray)
        xx = self.meanxintensity(angkbarray)
        return oo + xx, oo - xx

    def fluxIQ(self, energyarray):  # outgoing only
        xval = np.array([
            -1, -0.9340014304080591343323, -0.7844834736631444186224,
            -0.565235326996205006471, -0.2957581355869393914319, 0,
            0.2957581355869393914319, 0.565235326996205006471,
            0.7844834736631444186224, 0.9340014304080591343323, 1
        ])
        wval = np.array([
            0.01818181818181818181818, 0.1096122732669948644614,
            0.187169881780305204108, 0.2480481042640283140401,
            0.2868791247790080886792, 0.3002175954556906937859,
            0.286879124779008088679, 0.2480481042640283140401,
            0.1871698817803052041081, 0.109612273266994864461,
            0.018181818181818181818180
        ])
        xval = xval / 2.0 + 0.5
        wval = wval / 2
        theta = np.degrees(np.arccos(xval))
        phi = np.linspace(5, 175, 18)
        wf0val = np.kron(xval * wval, np.ones(len(phi))) * 2 * np.pi / len(phi)
        tt, ee, pp = np.meshgrid(theta, energyarray, phi)
        dd = [tt.flatten(), pp.flatten(), ee.flatten()]
        ii, qq = self.calcIQ(dd)
        # print(ii)
        # print(wf0val)
        ii = np.reshape(ii, (len(energyarray), len(xval) * len(phi)))
        # print(ii)
        qq = np.reshape(qq, (len(energyarray), len(xval) * len(phi)))
        return np.sum(ii * wf0val, axis=1), np.sum(qq * wf0val, axis=1)

    def calcvalue(self, dd, routine):
        if (routine == 'totalintensity'):
            return self.totalintensity(dd)
        elif (routine == 'xintensity'):
            return self.xintensity(dd)
        elif (routine == 'ointensity'):
            return self.ointensity(dd)
        elif (routine == 'meantotalintensity'):
            return self.meantotalintensity(dd)
        elif (routine == 'meanxintensity'):
            return self.meanxintensity(dd)
        elif (routine == 'meanointensity'):
            return self.meanointensity(dd)
        elif (routine == 'calcIQ'):
            return self.calcIQ(dd)
        elif (routine == 'calcmeanIQ'):
            return self.calcmeanIQ(dd)
        else:
            return 0
    def adjust_surface_temperature(self):
        fluxgoal=stefan_boltzmann_law(self.effective_temperature)
        # need to update surface temperature so that the outgoing flux is the effective temperature
        # do this iteratively ... converges very quickly (typically the surface temperature is 5-20% larger than Teff)
        epeak=4*self.effective_temperature
        for ii in range(25):
            self.eeloc=np.logspace(-3,1,151)*epeak
            self.iiloc,self.qqloc=self.fluxIQ(self.eeloc)
            flux=simps(self.iiloc,self.eeloc)*1.0000036200079545**4   # small correction to account for the range in the energy integral
            epeak=self.eeloc[np.argmax(self.iiloc*self.eeloc)]
            tlast=self.surface_temperature
            self.surface_temperature/=(flux/fluxgoal)**0.25
            # print(ii,self.surface_temperature,epeak)
            if np.abs((tlast-self.surface_temperature)/tlast)<1e-3:
                break
    def __str__(self):
        try:
            outstring='''#       
#   Energy[keV]          I          Q/I\n'''
            for ee,ii,rr in zip(self.eeloc,self.iiloc,self.qqloc/self.iiloc):
                outstring='%s %12g %12g %12g\n' % (outstring,ee,ii,rr)
        except AttributeError:
            outstring=''
        return outstring+'# '+str(super())+'\n'

from Magnetar.comptonized import partial_o_mode_comptonization, Complete_Comptonization
from Magnetar.atmo_lloyd import atmo_lloyd
from Magnetar.scatmo_caiazzo import scatmo_caiazzo


#
# Applies an energy dependent rescaling to an existing atmosphere
#
class rescale_atmo(atmosphere):
    def __init__(self, atmo=None, xfunk=lambda a: 1, ofunk=lambda a: 1):
        self.mag_inclination = 0
        self.loaddata(atmo, xfunk, ofunk)
    def __str__(self):
        outstring='''#
# class rescale_atmo 
#
# magnetic inclination %g
# start of rescaled atmosphere
%s
# end of rescaled atmosphere
''' % (self.mag_inclination,str(self.atmo))
        return outstring+atmosphere.__str__(self)
    def loaddata(self, atmo, xfunk, ofunk):
        self.atmo = atmo
        self.xfunk = xfunk
        self.ofunk = ofunk
        if atmo is None:
            return self
        else:
            self.mag_inclination = atmo.mag_inclination
            return self

    def xintensity(self, dataarray):
        return self.atmo.xintensity(dataarray) * self.xfunk(dataarray[-1])

    def ointensity(self, dataarray):
        return self.atmo.ointensity(dataarray) * self.ofunk(dataarray[-1])

    def meanxintensity(self, angkbarray):
        return self.atmo.meanxintensity(dataarray) * self.xfunk(dataarray[-1])

    def meanointensity(self, angkbarray):
        return self.atmo.meanointensity(dataarray) * self.ofunk(dataarray[-1])


#
# Applies a beaming function to an existing atmosphere
#
#
# If you use the xintensity or ointensity the angle is with respect to the vertical.
#
#
# If you use the meanxintensity or meanointensity the angle is with respect to the magnetic field.
#
#
class apply_beaming(atmosphere):
    def __init__(self, atmo=None, xfunk=lambda a: 1, ofunk=lambda a: 1):
        self.mag_inclination = 0
        self.loaddata(atmo, xfunk, ofunk)

    def loaddata(self, atmo, xfunk, ofunk):
        self.atmo = atmo
        self.xfunk = xfunk
        self.ofunk = ofunk
        if atmo is None:
            return self
        else:
            self.mag_inclination = atmo.mag_inclination
            return self

    def __str__(self):
        outstring='''#
# class apply_beaming
#
# magnetic inclination %g
# start of beamed atmosphere
%s
# end of beamed atmosphere
''' % (self.mag_inclination,str(self.atmo))
        return outstring+atmosphere.__str__(self)
    def xintensity(self, dataarray):
        return self.atmo.xintensity(dataarray) * self.xfunk(dataarray[-3])

    def ointensity(self, dataarray):
        return self.atmo.ointensity(dataarray) * self.ofunk(dataarray[-3])

    def meanxintensity(self, angkbarray):
        return self.atmo.meanxintensity(dataarray) * self.xfunk(dataarray[-2])

    def meanointensity(self, angkbarray):
        return self.atmo.meanointensity(dataarray) * self.ofunk(dataarray[-2])


class swap_xo_atmo(atmosphere):
    def __init__(self, atmo):
        self.atmo = atmo
        self.mag_inclination = atmo.mag_inclination

    def __str__(self):
        outstring='''#
# class swap_xo_atmo
#
# magnetic inclination %g
# start of swapped atmosphere
%s
# end of swapped atmosphere
''' % (self.mag_inclination,str(self.atmo))
        return outstring+atmosphere.__str__(self)
    def xintensity(self, dataarray):
        return self.atmo.ointensity(dataarray)

    def ointensity(self, dataarray):
        return self.atmo.xintensity(dataarray)

    def meanxintensity(self, angkbarray):
        return self.atmo.meanointensity(angkbarray)

    def meanointensity(self, angkbarray):
        return self.atmo.meanxintensity(angkbarray)


class add_atmo(atmosphere):
    def __init__(self, atmo1, atmo2):
        self.mag_inclination = atmo1.mag_inclination
        self.loaddata(atmo1, atmo2)

    def __str__(self):
        outstring='''#
# class add_atmo
#
# magnetic inclination %g
# start of first added atmosphere
%s
# end of first added atmosphere
# start of second added atmosphere
%s
# end of second added atmosphere
''' % (self.mag_inclination,str(self.atmo1),str(self.atmo2))
        return outstring+atmosphere.__str__(self)
    def loaddata(self, atmo1, atmo2):
        self.atmo1 = atmo1
        self.atmo2 = atmo2
        return self

    def xintensity(self, dataarray):
        return self.atmo1.xintensity(dataarray) + self.atmo2.xintensity(
            dataarray)

    def ointensity(self, dataarray):
        return self.atmo1.ointensity(dataarray) + self.atmo2.ointensity(
            dataarray)

    def meanxintensity(self, angkbarray):
        return self.atmo1.meanxintensity(angkbarray) + self.atmo2.xintensity(
            angkbarray)

    def meanointensity(self, angkbarray):
        return self.atmo1.meanointensity(angkbarray) + self.atmo2.ointensity(
            angkbarray)

