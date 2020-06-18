import numpy as np
from scipy import ndimage


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


# Comptonize the O-mode to Te in keV
def Complete_Comptonization(atmo, Te):
    mu = 7
    Teloc = np.abs(Te)
    ee = np.logspace(-2, 2, 201)
    ii, qq = atmo.fluxIQ(ee)
    totscat = simps(0.5 * (ii + qq) / ee, ee)
    # print(totscat)
    ee *= Teloc
    for i in range(5):
        nphot = 208452.792 * ee**2 / np.expm1(ee / Teloc + mu) / 2 * np.pi
        if (Te < 0):
            nphot *= 2.0 / 3.0
        mu = mu - np.log(totscat / simps(nphot, ee))
        # print(mu)
    if (Te > 0):
        catmo = be_atmo(Teloc, Teloc, mu, mu)
    else:
        catmo = compton_bb_atmo(Teloc, Teloc, mu, mu)
    catmo.mag_inclination = atmo.mag_inclination
    return atmo * (lambda e: 1, lambda e: 0) + catmo * (lambda e: 0,
                                                        lambda e: 1)


class atmo_lloyd(atmosphere):
    '''
  totalintensity: 
     calculate the total intensity for the atmo model, 
     dataarray contains the coordinates, direction and energy in the form
    [[zenith_ang1, zenith_ang2, zenith_ang3],
     [azimuth1,    azimuth2,    azimuth3],
     [energy1,     energy2,     energy3]]

  meantotalintensity:  
    calculate the mean total intensity for the atmo model around the field, 
    angkbarray contains the coordinates, direction and energy in the form
    [[field_ang1, field_ang2, field_ang3],
     [energy1,     energy2,     energy3]]

    '''

    def __init__(self):
        self.t = []
        self.x = []
        self.o = []
        self.mag_inclination = 0

    def loaddata(self, file):
        header = file.rsplit('.', 1)[0]
        enerkeV = np.loadtxt(header + '.gnu', usecols=(3), unpack=True)
        self.t, self.x, self.o = np.loadtxt(
            header + '.int', usecols=(2, 3, 4), unpack=True)
        theta, phi = np.loadtxt(header + '.trj', usecols=(4, 6), unpack=True)
        self.xxarray = []
        self.yyarray = []
        self.shape = []
        for dd in (theta, phi, enerkeV):
            xx = np.unique(dd)
            self.xxarray.append(xx)
            self.yyarray.append(np.arange(len(xx)))
            self.shape.append(len(xx))
        self.t = np.reshape(self.t, self.shape)
        self.x = np.reshape(self.x, self.shape)
        self.o = np.reshape(self.o, self.shape)
        return self

    def _calcindex(self, dataarray):
        res = []
        dataarray[1] = np.remainder(dataarray[1], 360)
        dataarray[1] = np.where(dataarray[1] < 180, dataarray[1],
                                360.0 - dataarray[1])
        for xx, yy, d in zip(self.xxarray, self.yyarray, dataarray):
            res.append(np.interp(d, xx, yy))
        return res

    def _intensity(self, dataarray, datacube):
        if (len(datacube) == 0):
            return 1
        else:
            res = ndimage.map_coordinates(
                datacube, self._calcindex(dataarray), order=3, mode='nearest')
            return np.where(res > 0, res, 0)

    def _meanintensity(self, angkbarray, datacube):
        if (len(datacube) == 0):
            return 1
        else:
            angkbarray = np.array(angkbarray)
            thetaminus = self.mag_inclination - angkbarray[0]
            thetaplus = self.mag_inclination + angkbarray[0]
            dataarrayminus = [
                np.abs(thetaminus),
                np.where(thetaminus > 0, 0, 180), angkbarray[1]
            ]
            dataarrayplus = [
                np.abs(thetaplus),
                np.where(thetaplus > 0, 0, 180), angkbarray[1]
            ]
            resplus = ndimage.map_coordinates(
                datacube,
                self._calcindex(dataarrayplus),
                order=3,
                mode='nearest')
            resplus = np.where(resplus > 0, resplus, 0)
            resminus = ndimage.map_coordinates(
                datacube,
                self._calcindex(dataarrayminus),
                order=3,
                mode='nearest')
            resminus = np.where(resminus > 0, resminus, 0)
            return 0.5 * (np.where(np.abs(thetaminus) < 90, resminus, 0) +
                          np.where(np.abs(thetaplus) < 90, resplus, 0))

    def totalintensity(self, dataarray):
        return self._intensity(dataarray, self.t)

    def xintensity(self, dataarray):
        return self._intensity(dataarray, self.x)

    def ointensity(self, dataarray):
        return self._intensity(dataarray, self.o)

    def meantotalintensity(self, angkbarray):
        return self._meanintensity(angkbarray, self.t)

    def meanxintensity(self, angkbarray):
        return self._meanintensity(angkbarray, self.x)

    def meanointensity(self, angkbarray):
        return self._meanintensity(angkbarray, self.o)


def bb_atmo_f(temp):
    return atmosphere() * (lambda ee: 208452.792 * ee**3 / np.expm1(ee / temp))


#
# blackbody atmosphere
#
class bb_atmo(atmosphere):
    def __init__(self):
        self.atmo = []
        self.xtemp = 1.0
        self.otemp = 1.0
        self.mag_inclination = 0

    def loaddata(self, xtemp, otemp):
        self.xtemp = xtemp
        self.otemp = otemp
        return self

    def _bbfunk(self, ee, tt):  # per mode
        return 208452.792 * ee**3 / np.expm1(ee / tt) / 2

    def xintensity(self, dataarray):
        return self._bbfunk(dataarray[-1], self.xtemp)

    def ointensity(self, dataarray):
        return self._bbfunk(dataarray[-1], self.otemp)


#
# Fully Comptonized blackbody atmosphere (photons with fixed number in equilibrium at TX, TO)
#
class be_atmo(atmosphere):  # without beaming
    def __init__(self, xtemp, otemp, xmu, omu):
        self.atmo = []
        self.xtemp = 1.0
        self.otemp = 1.0
        self.xmu = 1.0
        self.omu = 1.0
        self.mag_inclination = 0
        self.loaddata(xtemp, otemp, xmu, omu)

    def loaddata(self, xtemp, otemp, xmu, omu):
        self.xtemp = xtemp
        self.otemp = otemp
        self.xmu = xmu
        self.omu = omu
        return self

    def _befunk(self, ee, tt, mu):  # per mode
        return 208452.792 * ee**3 / np.expm1(ee / tt + mu) / 2

    def xintensity(self, dataarray):
        return self._befunk(dataarray[-1], self.xtemp, self.xmu)

    def ointensity(self, dataarray):
        return self._befunk(dataarray[-1], self.otemp, self.omu)


#
# Fully Comptonized blackbody atmosphere (photons with fixed number in equilibrium at TX, TO)
#
class compton_bb_atmo(be_atmo):  # with costheta beaming
    def xintensity(self, dataarray):
        return self._befunk(dataarray[-1], self.xtemp, self.xmu) * np.cos(
            np.radians(dataarray[0]))

    def ointensity(self, dataarray):
        return self._befunk(dataarray[-1], self.otemp, self.omu) * np.cos(
            np.radians(dataarray[0]))


#
# Applies an energy dependent rescaling to an existing atmosphere
#
class rescale_atmo(atmosphere):
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


from scipy.integrate import simps


class partial_o_mode_comptonization(atmosphere):
    def __init__(self, atmo, y):
        self.mag_inclination = atmo.mag_inclination
        self.y = y
        self.atmo = atmo
        self.ee = np.logspace(-2, 2, 401)
        ii, qq = atmo.fluxIQ(self.ee)
        oof = (ii + qq) / 2 / self.ee**3
        self.plarge = -1.5 - (2.25 + 4 / np.abs(y))**0.5
        self.psmall = -1.5 + (2.25 + 4 / np.abs(y))**0.5
        kernelval = np.where(self.ee < 1, self.ee**self.psmall, self.ee
                             **self.plarge)
        kernelval = kernelval / np.sum(kernelval)
        self.ooout = np.convolve(
            oof, kernelval, mode='same') * self.ee**3 / np.pi
        if y < 0:
            iif, qqf = self.fluxIQ(self.ee)
            self.ooout *= simps((ii + qq) / self.ee, self.ee) / simps(
                (iif + qqf) / self.ee, self.ee)

    def xintensity(self, dataarray):
        return self.atmo.xintensity(dataarray)

    def meanxintensity(self, angkbarray):
        return self.atmo.meanxintensity(angkbarray)

    def meanointensity(self, angkbarray):
        if self.y < 0:
            return np.interp(angkbarray[-1], self.ee, self.ooout) * np.sin(
                np.radians(angkbarray[-2]))**2
        else:
            return np.interp(angkbarray[-1], self.ee, self.ooout)

    def ointensity(self, dataarray):
        if self.y < 0:
            cosangkb = np.cos(np.radians(self.mag_inclination)) * np.cos(
                np.radians(dataarray[0])
            ) + np.sin(np.radians(self.mag_inclination)) * np.sin(
                np.radians(dataarray[0])) * np.cos(np.radians(dataarray[1]))
            return np.interp(dataarray[-1], self.ee, self.ooout) * (
                1 - cosangkb**2)
        else:
            return np.interp(dataarray[-1], self.ee, self.ooout)


class partial_res_comptonization(atmosphere):
    def __init__(self, atmo, y, kTe):
        self.mag_inclination = atmo.mag_inclination
        self.y = y
        self.atmo = atmo
        self.ee = np.logspace(-2, 2, 401)
        self.ii, self.qq = atmo.fluxIQ(self.ee)
        qfrac = self.qq / self.ii
        iif = self.ii / self.ee**3
        self.plarge = -1.5 - (2.25 + 4 / np.abs(y))**0.5
        self.psmall = -1.5 + (2.25 + 4 / np.abs(y))**0.5
        kernelval = np.where(self.ee < 1, self.ee**self.psmall, self.ee
                             **self.plarge)
        kernelval = kernelval / np.sum(kernelval)
        self.iiout = np.convolve(
            iif, kernelval, mode='same') * self.ee**3 / np.pi
        self.nlndeltaenout = np.convolve(
            iif, kernelval * np.abs(np.log(self.ee)),
            mode='same') * self.ee**3 / np.pi / self.iiout
        self.meanscat = np.abs(
            self.nlndeltaenout / np.log(1 + 4 * (kTe / 511.0)))
        self.frac = np.exp(-self.meanscat)
        self.qqout = self.iiout * np.exp(-0.5 * self.meanscat) * np.interp(
            self.ee * np.exp(-self.nlndeltaenout), self.ee, qfrac)
        self.ooout = 0.5 * (self.iiout + self.qqout)
        self.xxout = 0.5 * (self.iiout - self.qqout)

    def calcIQ(self, dataarray):
        xxloc = self.xintensity(dataarray)
        ooloc = self.ointensity(dataarray)
        return ooloc + xxloc, ooloc - xxloc

    def xintensity(self, dataarray):
        if self.y < 0:
            frac = np.interp(dataarray[-1], self.ee, self.frac)
            resx = np.interp(dataarray[-1], self.ee, self.xxout)
            rat = np.interp(dataarray[-1], self.ee,
                            np.pi * self.xxout / (0.5 * (self.ii - self.qq)))
            return rat * self.atmo.xintensity(dataarray) * frac + resx * (
                1 - frac)
        else:
            return np.interp(dataarray[-1], self.ee, self.xxout)

    def ointensity(self, dataarray):
        if self.y < 0:
            cosangkb = np.cos(np.radians(self.mag_inclination)) * np.cos(
                np.radians(dataarray[0])
            ) + np.sin(np.radians(self.mag_inclination)) * np.sin(
                np.radians(dataarray[0])) * np.cos(np.radians(dataarray[1]))
            frac = np.interp(dataarray[-1], self.ee, self.frac)
            reso = np.interp(dataarray[-1], self.ee, self.ooout)
            rat = np.interp(dataarray[-1], self.ee,
                            np.pi * self.ooout / (0.5 * (self.ii + self.qq)))
            return rat * self.atmo.ointensity(dataarray) * frac + reso * (
                1 - frac) * cosangkb**2
        else:
            return np.interp(dataarray[-1], self.ee, self.ooout)


class partial_twisted_comptonization(atmosphere):
    def __init__(self, atmo, y, kTe, ehuge, phuge):
        self.mag_inclination = atmo.mag_inclination
        self.y = y
        self.atmo = atmo
        self.ee = np.logspace(-3, 3, 601)
        self.ii, self.qq = atmo.fluxIQ(self.ee)
        qfrac = self.qq / self.ii
        iif = self.ii / self.ee**3
        self.phuge = phuge
        self.ehuge = ehuge
        self.plarge = -1.5 - (2.25 + 4 / np.abs(y))**0.5
        self.psmall = -1.5 + (2.25 + 4 / np.abs(y))**0.5
        kernelval = np.where(
            self.ee < 1, self.ee**self.psmall,
            np.where(self.ee < ehuge, self.ee**self.plarge,
                     (ehuge)**self.plarge * (self.ee / ehuge)**self.phuge))
        self.kernelval = kernelval / np.sum(kernelval)
        self.iiout = np.convolve(
            iif, self.kernelval, mode='same') * self.ee**3 / np.pi
        self.nlndeltaenout = np.convolve(
            iif, self.kernelval * np.abs(np.log(self.ee)),
            mode='same') * self.ee**3 / np.pi / self.iiout
        transfunk = 0.5 + 0.5 * np.tanh((self.nlndeltaenout - np.log(ehuge)))
        self.meanscat = transfunk * 2.0 + (1 - transfunk) * np.abs(
            self.nlndeltaenout / np.log(1 + 4 * (kTe / 511.0)))
        self.frac = np.exp(-self.meanscat)
        self.qqout = self.iiout * np.exp(-0.5 * self.meanscat) * np.interp(
            self.ee * np.exp(-self.nlndeltaenout), self.ee, qfrac)
        self.ooout = 0.5 * (self.iiout + self.qqout)
        self.xxout = 0.5 * (self.iiout - self.qqout)

    def calcIQ(self, dataarray):
        xxloc = self.xintensity(dataarray)
        ooloc = self.ointensity(dataarray)
        return ooloc + xxloc, ooloc - xxloc

    def xintensity(self, dataarray):
        if self.y < 0:
            frac = np.interp(dataarray[-1], self.ee, self.frac)
            resx = np.interp(dataarray[-1], self.ee, self.xxout)
            rat = np.interp(dataarray[-1], self.ee,
                            np.pi * self.xxout / (0.5 * (self.ii - self.qq)))
            return rat * self.atmo.xintensity(dataarray) * frac + resx * (
                1 - frac)
        else:
            return np.interp(dataarray[-1], self.ee, self.xxout)

    def ointensity(self, dataarray):
        if self.y < 0:
            cosangkb = np.cos(np.radians(self.mag_inclination)) * np.cos(
                np.radians(dataarray[0])
            ) + np.sin(np.radians(self.mag_inclination)) * np.sin(
                np.radians(dataarray[0])) * np.cos(np.radians(dataarray[1]))
            frac = np.interp(dataarray[-1], self.ee, self.frac)
            reso = np.interp(dataarray[-1], self.ee, self.ooout)
            rat = np.interp(dataarray[-1], self.ee,
                            np.pi * self.ooout / (0.5 * (self.ii + self.qq)))
            return rat * self.atmo.ointensity(dataarray) * frac + reso * (
                1 - frac) * cosangkb**2
        else:
            return np.interp(dataarray[-1], self.ee, self.ooout)


class scatmo_caiazzo(atmosphere):
    '''
 
  calcIQ:  
    calculate the mean total intensity and Q for the scatmo model around the field, 
    angkbarray contains the coordinates, direction and energy in the form
    [[field_ang1, field_ang2, field_ang3],
     [energy1,     energy2,     energy3]]

    '''

    def __init__(self):
        self.Eval = [0.1, 100]
        self.i0 = [1, 1]
        self.i1 = [1, 1]
        self.i2 = [1, 1]
        self.q0 = [1, 1]
        self.q1 = [1, 1]
        self.q2 = [1, 1]
        self.mag_inclination = 0

    @staticmethod
    def f0(x):
        return 1 / 4. * 15.**(1. / 2) * (-x**2 + 1.)

    @staticmethod
    def f1(x):
        return 1. / 2 * x * 6.**(1. / 2)

    @staticmethod
    def f2(x):
        return 5. / 4 * 3.**(1. / 2) * (x**2 - 1. / 5)

    def loaddata(self, file):
        self.Eval, self.i0, self.i1, self.i2, self.q0, self.q1, self.q2 = np.loadtxt(
            file, unpack=True, usecols=range(7))

    def savedata(self, file):
        np.savetxt(file,
                   np.transpose([
                       self.Eval, self.i0, self.i1, self.i2, self.q0, self.q1,
                       self.q2
                   ]))

    #
    # Use an atmosphere model to create a scatmosphere model by calculating the various expansion coefficients
    #   https://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss%E2%80%93Lobatto_rules
    #
    def createscatmo(self, atmo, Eval, xval=None, wval=None):
        self.Eval = np.unique(Eval)
        if xval is None:
            # Use Gauss-Lobatto Quadrature with 11 points in angle
            self._xval = np.array([
                -1, -0.9340014304080591343323, -0.7844834736631444186224,
                -0.565235326996205006471, -0.2957581355869393914319, 0,
                0.2957581355869393914319, 0.565235326996205006471,
                0.7844834736631444186224, 0.9340014304080591343323, 1
            ])
            self._wval = np.array([
                0.01818181818181818181818, 0.1096122732669948644614,
                0.187169881780305204108, 0.2480481042640283140401,
                0.2868791247790080886792, 0.3002175954556906937859,
                0.286879124779008088679, 0.2480481042640283140401,
                0.1871698817803052041081, 0.109612273266994864461,
                0.018181818181818181818180
            ])
        else:
            self._xval = np.array(xval)
            if wval is None:
                self._wval = np.full(len(xval), 1.0 / len(xval))
            else:
                self._wval = np.array(wval)
        self._angkb = np.degrees(np.arccos(self._xval))
        wf0val = self.f0(self._xval) * self._wval
        wf1val = self.f1(self._xval) * self._wval
        wf2val = self.f2(self._xval) * self._wval
        tt, ee = np.meshgrid(self._angkb, self.Eval)
        dd = [tt.flatten(), ee.flatten()]
        ii, qq = atmo.calcmeanIQ(dd)
        ii = np.reshape(ii, (len(self.Eval), len(self._xval)))
        qq = np.reshape(qq, (len(self.Eval), len(self._xval)))
        self.i0 = np.sum(ii * wf0val, axis=1)
        self.i1 = np.sum(ii * wf1val, axis=1)
        self.i2 = np.sum(ii * wf2val, axis=1)
        self.q0 = np.sum(qq * wf0val, axis=1)
        self.q1 = np.sum(qq * wf1val, axis=1)
        self.q2 = np.sum(qq * wf2val, axis=1)

    def calcmeanIQ(self, angkbarray):
        ee = angkbarray[1]
        x = np.cos(np.radians(angkbarray[0]))
        f0val = self.f0(x)
        f1val = self.f1(x)
        f2val = self.f2(x)
        return np.interp(ee, self.Eval, self.i0) * f0val + np.interp(
            ee, self.Eval, self.i1) * f1val + np.interp(
                ee, self.Eval, self.i2) * f2val, np.interp(
                    ee, self.Eval, self.q0) * f0val + np.interp(
                        ee, self.Eval, self.q1) * f1val + np.interp(
                            ee, self.Eval, self.q2) * f2val

    def meantotalintensity(self, angkbarray):
        ii, qq = self.calcmeanIQ(angkbarray)
        return ii

    def meanxintensity(self, angkbarray):
        ii, qq = self.calcmeanIQ(angkbarray)
        return 0.5 * (ii - qq)

    def meanointensity(self, angkbarray):
        ii, qq = self.calcmeanIQ(angkbarray)
        return 0.5 * (ii + qq)

    def _calc_angkbarray(self, dataarray):
        angkb = np.degrees(
            np.arccos(
                np.cos(np.radians(self.mag_inclination)
                       ) * np.cos(np.radians(dataarray[0])) +
                np.sin(np.radians(self.mag_inclination)
                       ) * np.sin(np.radians(dataarray[0])
                                  ) * np.cos(np.radians(dataarray[1]))))
        angkbarray = dataarray[1:]  # remove the first angle
        angkbarray[0] = angkb  # replace the second angle with the field angle
        return angkbarray

    def calcIQ(self, dataarray):
        return self.calcmeanIQ(self._calc_angkbarray(dataarray))

    def totalintensity(self, dataarray):
        return self.meantotalintensity(self._calc_angkbarray(dataarray))

    def xintensity(self, dataarray):
        return self.meanxintensity(self._calc_angkbarray(dataarray))

    def ointensity(self, dataarray):
        return self.meanointensity(self._calc_angkbarray(dataarray))


class surface_model(atmosphere):
    def __init__(self):
        self.patches = []
        self.mcolat = []

    def loaddata(self, files, modeltype=None):
        if type(files) is str:
            files = [
                files,
            ]
        for ff in files:
            # print(ff)
            ang = float((ff.rsplit('/', 1)[-1]).rsplit('_', 2)[0])
            self.add_patch(atmosphere().loaddata(ff, modeltype), ang)
        return self.sort_patches()

    def add_patch(self, atmo, ang):
        self.mcolat.append(ang)
        hld = np.cos(np.radians(ang))**2
        hld = 4.0 * hld / (3.0 * hld + 1.0)
        atmo.mag_inclination = np.degrees(np.arccos(hld**0.5))
        self.patches.append(atmo)
        return self

    def sort_patches(self):
        if len(self.mcolat) > 1:
            dum, ii = np.unique(self.mcolat, return_index=True)
            ns = []
            mm = []
            for i in ii:
                ns.append(self.patches[i])
                mm.append(self.mcolat[i])
            self.patches = ns
            self.mcolat = mm
        return self

    def load_lloyd_data(self, files):
        return self.loaddata(files, 'Lloyd')

    def load_caiazzo_data(self, files):
        return self.loaddata(files, 'Caiazzo')

    '''
    calculate the total intensity for the surface model, dataarray contains the coordinates, direction and energy in the form
    [[
    1,     mcolat2,     mcolat3],
     [zenith_ang1, zenith_ang2, zenith_ang3],
     [azimuth1,    azimuth2,    azimuth3],
     [energy1,     energy2,     energy3]]
     
     where field_mu is the cosine of angle between field at that position and the normal,
           zenith_ang is the angle between the line of sight and the vertical at emission in degrees,
           azimuth1 is the angle between the plane containing k+normal and k+B in degrees
           energy1 is the energy in keV
           
     dataarray can contain as many columns as you want
   meantotalintensity:  
    calculate the mean total intensity for the atmo model around the field, 
    angkbarray contains the coordinates, direction and energy in the form
    [[mcolat1,     mcolat2,     mcolat3],
     [field_ang1, field_ang2, field_ang3],
     [energy1,     energy2,     energy3]]

           

    '''

    def _dointerpolate(self, dataarray, res):
        res = np.array(res)
        # print(res)

        ii = np.interp(dataarray[0], self.mcolat, np.arange(len(self.mcolat)))
        iif, iid = np.modf(ii)
        iid = iid.astype(int)
        iid = np.clip(iid, 0, len(self.mcolat) - 2)
        iif = np.where(ii <= 0, 0, np.where(ii >= len(self.mcolat) - 1, 1,
                                            iif))
        resout = 0 * iif
        cnt = np.arange(len(resout))
        resout[cnt] = res[iid, cnt] * (1 - iif) + res[iid + 1, cnt] * iif
        return resout
        '''
        print(np.shape(self.mcolat))
        print(np.shape(res))
        print(np.shape(res[:,10]))
        print(dataarray[0][10])
        resout=dataarray[0]*0
        for i in range(len(resout)):
            resout[i]=np.interp(dataarray[0][i],self.mcolat,res[:,i])
        return resout
        '''

    def _interpolate_single(self, dataarray, routine):
        if (len(self.patches) == 0):
            return 1
        elif len(self.patches) == 1:
            return self.patches[0].calcvalue(dataarray[1:], routine)
        else:
            res = []
            dd = dataarray[1:]
            for i, mu in enumerate(self.mcolat):
                res.append(self.patches[i].calcvalue(dd, routine))
            return self._dointerpolate(dataarray, res)

    def _interpolate_double(self, dataarray, routine):
        if (len(self.patches) == 0):
            return 1, 1
        elif len(self.patches) == 1:
            return self.patches[0].calcvalue(dataarray[1:], routine)
        else:
            resi = []
            resq = []
            dd = dataarray[1:]
            for i, mu in enumerate(self.mcolat):
                ii, qq = self.patches[i].calcvalue(dd, routine)
                resi.append(ii)
                resq.append(qq)
            return self._dointerpolate(dataarray, resi), self._dointerpolate(
                dataarray, resq)

    def totalintensity(self, dataarray):
        return self._interpolate_single(dataarray, 'totalintensity')

    def xintensity(self, dataarray):
        return self._interpolate_single(dataarray, 'xintensity')

    def ointensity(self, dataarray):
        return self._interpolate_single(dataarray, 'ointensity')

    def calcIQ(self, dataarray):
        return self._interpolate_double(dataarray, 'calcIQ')

    def meantotalintensity(self, angkbarray):
        return self._interpolate_single(angkbarray, 'meantotalintensity')

    def meanxintensity(self, dataarray):
        return self._interpolate_single(angkbarray, 'meanxintensity')

    def meanointensity(self, dataarray):
        return self._interpolate_single(angkbarray, 'meanointensity')

    def calcmeanIQ(self, angkbarray):
        return self._interpolate_double(angkbarray, 'calcmeanIQ')


class pfield:
    def __init__(self):
        self.data = []
        self.imean = 0
        self.qmean = 0
        self.mu = 0
        self.nu = 0
        self.mass = 0
        self.radius = 0
        self.theta = 0
        self.angkb = []
        self.ndotb = []
        self.pt = 1
        self.filename = ""
        self.ebins=[]
        self.iint=[]
        self.qint=[]        
        self.oneSided=True

    def loaddata(self, pfield_file):
        with open(pfield_file, 'r') as f:
            line = f.readline()
            a = line.split()
            try:
                self.mu = float(a[2])
                self.nu = float(a[3])
                self.mass = float(a[4])
                self.radius = float(a[5])
                self.theta = float(a[6])
            except:
                pass
        self.data = np.genfromtxt(
            pfield_file, unpack=True, names=True, skip_header=25)
        self.data = self.data[~np.isnan(self.data['s1'])]
        self.imean = np.mean(self.data['X'] + self.data['O'])
        self.qmean = np.mean(self.data['Q'] *
                             (self.data['X'] - self.data['O']))
        hld = np.cos(np.radians(self.data['mag_colat']))**2
        self.ndotb = (4 * hld / (3 * hld + 1))**0.5
        self.angkb = np.degrees(
            np.arccos((self.ndotb * np.cos(np.radians(self.data['theta'])) +
                       (1 - self.ndotb**2)**0.5 * np.sin(
                           np.radians(self.data['theta'])) * np.cos(
                               np.radians(self.data['phi'])))))
        self.filename=pfield_file
        return self

    def calcIQ(self, energies, atmo_map):
        ones = np.ones(len(energies))
        dataarray = np.array([
            np.kron(ones, self.data['theta']),
            np.kron(ones, self.data['phi']),
            np.kron(energies, np.ones(len(self.data['phi'])))
        ])
        ii, qq = atmo_map.calcIQ(dataarray)
        ii = np.reshape(ii, (len(energies), len(self.data['phi'])))
        qq = np.reshape(qq, (len(energies), len(self.data['phi'])))
        dotprod = self.data['s1'] * np.cos(
            2 * self.data['beta']) - self.data['s2'] * np.sin(
                2 * self.data['beta'])
        return np.mean(ii, axis=1), np.mean(-qq * dotprod, axis=1)

    def calcmeanIQ(self, energies, atmo_map):
        ones = np.ones(len(energies))
        dataarray = np.array([
            np.kron(ones, self.data['mag_colat']),
            np.kron(ones, self.angkb),
            np.kron(energies, np.ones(len(self.angkb)))
        ])
        ii, qq = atmo_map.calcmeanIQ(dataarray)
        ii = np.reshape(ii, (len(energies), len(self.angkb)))
        qq = np.reshape(qq, (len(energies), len(self.angkb)))
        dotprod = self.data['s1'] * np.cos(
            2 * self.data['beta']) - self.data['s2'] * np.sin(
                2 * self.data['beta'])
        return np.mean(ii, axis=1), -np.mean(qq * dotprod, axis=1)

    def recalculate(self, energy, atmo_map, gtt=1):
        dataarray = np.array([
            self.data['mag_colat'], self.data['theta'], self.data['phi'],
            self.pt * np.full(len(self.data['mag_colat']), energy / gtt)
        ])
        ii, qq = atmo_map.calcIQ(dataarray)
        self.data['X'] = 0.5 * (ii - qq)
        self.data['O'] = 0.5 * (ii + qq)
        self.imean = np.mean(ii) * gtt**3
        dotprod = self.data['s1'] * np.cos(
            2 * self.data['beta']) - self.data['s2'] * np.sin(
                2 * self.data['beta'])
        self.qmean = np.mean(-dotprod * qq) * gtt**3
        #self.qmean=np.mean(-self.data['Q']*qq)*gtt**3
        return self.imean, self.qmean

    def recalculatephiaverage(self, energy, atmo_map, gtt=1):
        dataarray = [
            self.data['mag_colat'], self.angkb,
            self.data['mag_colat'] * 0 + energy / gtt
        ]
        ii, qq = atmo_map.calcmeanIQ(dataarray)
        self.data['X'] = (ii - qq) * 0.5
        self.data['O'] = (ii + qq) * 0.5
        self.imean = np.mean(ii) * gtt**3
        dotprod = self.data['s1'] * np.cos(
            2 * self.data['beta']) - self.data['s2'] * np.sin(
                2 * self.data['beta'])
        self.qmean = np.mean(-dotprod * qq) * gtt**3
        return self.imean, self.qmean

    def calcvalues(self,surfacemodel,ebins=None,gtt=1):
        if ebins is None:
            ebins=np.logspace(-0.5,1.5,100)
        self.ebins=ebins
        self.iint=[]
        self.qint=[]
        for en in ebins:
            self.recalculate(en,surfacemodel,gtt=gtt)
            self.iint.append(self.imean)
            self.qint.append(self.qmean)
        self.iint=np.array(self.iint)
        self.qint=np.array(self.qint)
    def __str__(self):
        outstring='# Phi[rad]             Energy[keV]    I             Q/I\n'
        for ee,ii,rr in zip(self.ebins,self.iint,self.qint/self.iint):
            outstring='%s%12g %12g %12g %12g\n' % (outstring,self.theta*np.pi/180.0,ee,ii,rr)
        return outstring
        
    def plot(self,
             maxrings=15,
             size=0.08,
             plotpoints=False,
             drawedge=True,
             drawaxes=False,
             contours=[15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180],
             lev=100,
             datamap=None,
             cmap="viridis",
             ellipsecolor=[0.9, 0.0, 0.0],
             contourcolor='k',
             oneSided=True):
        import matplotlib.pyplot as plt
        from matplotlib.patches import Ellipse

        def _doellipse(bb, be, s1a, s2a, s3a):
            e = Ellipse(
                xy=[bb * np.cos(be), bb * np.sin(be)],
                width=size,
                height=size * s3a / (s1a * s1a + s2a * s2a + s3a * s3a)**0.5,
                angle=(np.arctan2(s2a, s1a) * 0.5 + be) * 180 / 3.1415)
            ax.add_artist(e)
            e.set_clip_box(ax.bbox)
            e.set_alpha(0.6)
            e.set_facecolor(ellipsecolor)
            e.set_edgecolor(ellipsecolor)
            e.set_linewidth(0.5)
            return e

        fig = plt.figure(0, figsize=(8, 8))
        ax = fig.add_subplot(111, aspect='equal')
        b, beta, s1, s2, s3, mag_colat = self.data['b'], self.data[
            'beta'], self.data['s1'], self.data['s2'], self.data[
                's3'], self.data['mag_colat']
        bmax = max(b)
        x = b * np.cos(beta) / bmax
        y = b * np.sin(beta) / bmax
        if contours is not None:
            ax.tricontour(
                np.concatenate((x, x)),
                np.concatenate((y, -y)),
                np.concatenate((mag_colat, mag_colat)),
                levels=contours,
                linewidths=0.5,
                colors=contourcolor)
        if datamap is None:
            if drawedge:
                e = Ellipse(xy=[0, 0], width=2, height=2)
                ax.add_artist(e)
                e.set_clip_box(ax.bbox)
                e.set_alpha(1.0)
                e.set_facecolor([1.0, 1.0, 1.0])
                e.set_edgecolor([0.0, 0.0, 0.0])
                e.set_linewidth(1.0)
        else:
            if oneSided:
                cntr2 = ax.tricontourf(
                    np.concatenate((x, x)),
                    np.concatenate((y, -y)),
                    np.concatenate((datamap, datamap)),
                    levels=lev,
                    cmap=cmap)
            else:
                cntr2 = ax.tricontourf(x, y, datamap, levels=lev, cmap=cmap)
            fig.colorbar(cntr2, ax=ax)

        bmax = max(b)
        nb = (int)(len(b) / 2)**0.5
        if nb > maxrings:
            skip = (int)(nb / maxrings)
            nbf = nb / skip
            nstart = np.arange(0, nbf + 1)
            istart = nstart * nstart * skip * skip * 2
            iend = (nstart * skip + 1) * (nstart * skip + 1) * 2
            ival = []
            for ii, jj in zip(istart, iend):
                ival = np.concatenate((ival, np.arange(ii + skip / 2, jj,
                                                       skip)))
            ival = ival[ival < len(b)].astype(int)
            b, beta, s1, s2, s3 = b[ival], beta[ival], s1[ival], s2[ival], s3[
                ival]

        if plotpoints:
            ax.plot(b / bmax * np.cos(beta), b / bmax * np.sin(beta), 'b.')
            if oneSided:
                ax.plot(b / bmax * np.cos(beta), -b / bmax * np.sin(beta),
                        'b.')
        for bb, be, s1a, s2a, s3a in zip(b / bmax, beta, s1, s2, np.abs(s3)):
            _doellipse(bb, be, s1a, s2a, s3a)
            if oneSided:
                _doellipse(bb, -be, s1a, -s2a, s3a)

        if drawaxes == False:
            plt.axis('off')

        return fig

class pfield_array:
    def __init__(self):
        self.mufil,self.pfi = [],[]
    def loaddata(self, files, modeltype=None):
        self.mufil,self.pfi = [],[]
        if type(files) is str:
            files = [
                files,
            ]
        for i in files:
            pfii= pfield()
            pfii.loaddata(i)
            self.pfi.append(pfii)
            self.mufil.append(pfii.theta)
        ind = np.argsort(self.mufil)
        self.mufil = [self.mufil[i] for i in ind]
        self.pfi = [self.pfi[i] for i in ind]
        return self
    def calcvalues(self,surfacemodel,ebins=None,gtt=1):
        for pf in self.pfi:
            pf.calcvalues(surfacemodel,ebins,gtt=1)
    def __str__(self):
        outstring=''
        for pf in self.pfi:
            outstring=outstring+str(pf)
        return outstring
    def recalculate(self,energy,surfacemodel,gtt=1):
        for pf in self.pfi:
            pf.recalculate(energy,surfacemodel,gtt=gtt)
    def loadfiles(self,dirpath):
        import os
        import fnmatch
        from tkinter import Tcl
        import glob
        import matplotlib.pylab as plt

        # search pfield files in input directory
        filelist = []
        for path,dirs,files in os.walk(dirpath):
            for file in files:
                if fnmatch.fnmatch(file,'pfield_angle*'):
                    filename = os.path.join(path,file)
                    filelist.append(filename)

        # sort pfield files in accending order according angle in file name
        filelist = Tcl().call('lsort', '-dict', filelist)
        #define energy array and asociated NH absorption
        self.e_pfield   = np.logspace(-0.5,2.0,30)
        eabs,sabs=np.loadtxt('2019-07-17/tbabs.dat',unpack=True)
        sabs=sabs
        ii=np.argsort(eabs)
        eabs=eabs[ii]
        sabs=sabs[ii]
        # absorption cross-section per hydrogen atom in units of 1e-24 cm^2 for our energy bins
        ssabs=np.interp(self.e_pfield,eabs,sabs)/(self.e_pfield)**3*1e-24
        # the best-fit hydrogen column is 0.52e22 /cm2
        totabs=np.exp(-0.6e22*ssabs)
        # load surface emission 2d-array (energy and theta)
        self.th_pfield  = []
        self.i_pfield2d = []
        self.iabs_pfield2d = []
        self.q_pfield2d = []
        fig, ax1 = plt.subplots()
        fig, ax2 = plt.subplots()
        for files in filelist:
            #load pfield files
            pfieldaux = pfield()
            pfieldaux.loaddata(files)
            self.th_pfield.append(pfieldaux.theta)
            #load atmophereric surface model
            surfaceaux = surface_model()
            surfaceaux.loaddata(glob.glob('../QEDSurface/doubleBB_h/*.int'))
            aa=surfaceaux.mcolat
            #add the angles
            surfaceaux.mcolat=aa+[180-i for i in aa[::-1]]
            #add the patches
            surfaceaux.patches=surfaceaux.patches+surfaceaux.patches[::-1]
            #recalculate surface emission model
            ivec_abs, ivec, qvec = [], [], []
            for ee in self.e_pfield:
                pfieldaux.recalculate(ee,surfaceaux,gtt=(1-2*2.0/10.0)**0.5)
                ivec.append(pfieldaux.imean)
                qvec.append(pfieldaux.qmean)
            #ivec = np.array(ivec)
            #qvec = np.array(qvec)
            self.i_pfield2d.append(ivec)
            self.iabs_pfield2d.append(totabs*ivec/1.7e4)
            self.q_pfield2d.append(qvec)
            #plt.loglog(self.e_pfield,ivec/self.e_pfield**2)
            if(pfieldaux.theta<=180.0):
                ax1.plot(self.e_pfield,np.array(totabs*ivec/1.7e4))
                ax2.plot(self.e_pfield,np.array(qvec)/np.array(ivec))
        #plt.ylim(1,1e4)
        ax1.set_yscale("log")
        plt.show()
        self.i_pfield2d = np.array(self.i_pfield2d)
        self.iabs_pfield2d = np.array(self.iabs_pfield2d)
        self.q_pfield2d = np.array(self.q_pfield2d)
        return self.e_pfield

    def loadfullcompton(self,dirpath):
        import os
        import fnmatch
        from tkinter import Tcl
        import glob
        import matplotlib.pylab as plt

        # search pfield files in input directory
        filelist = []
        for path,dirs,files in os.walk(dirpath):
            for file in files:
                if fnmatch.fnmatch(file,'pfield_angle*'):
                    filename = os.path.join(path,file)
                    filelist.append(filename)

        # sort pfield files in accending order according angle in file name
        filelist = Tcl().call('lsort', '-dict', filelist)
        #define energy array and asociated NH absorption
        self.e_pfield   = np.logspace(-0.5,2.0,30)
        eabs,sabs=np.loadtxt('2019-07-17/tbabs.dat',unpack=True)
        sabs=sabs
        ii=np.argsort(eabs)
        eabs=eabs[ii]
        sabs=sabs[ii]
        # absorption cross-section per hydrogen atom in units of 1e-24 cm^2 for our energy bins
        ssabs=np.interp(self.e_pfield,eabs,sabs)/(self.e_pfield)**3*1e-24
        # the best-fit hydrogen column is 0.52e22 /cm2
        totabs=np.exp(-0.6e22*ssabs)
        # load surface emission 2d-array (energy and theta)
        self.th_pfield  = []
        self.i_pfield2d = []
        self.iabs_pfield2d = []
        self.q_pfield2d = []
        fig, ax1 = plt.subplots()
        fig, ax2 = plt.subplots()
        for files in filelist:
            #load pfield files
            pfieldaux = pfield()
            pfieldaux.loaddata(files)
            self.th_pfield.append(pfieldaux.theta)
            #load atmophereric surface model
            surfaceaux = surface_model()
            surfaceaux.loaddata(glob.glob('../QEDSurface/doubleBB_h/*.int'))
            aa=surfaceaux.mcolat
            #add the angles
            surfaceaux.mcolat=aa+[180-i for i in aa[::-1]]
            #add the patches
            surfaceaux.patches=surfaceaux.patches+surfaceaux.patches[::-1]
            #repeat
            allsurface150=surface_model().loaddata(glob.glob('../QEDSurface/doubleBB_h/*.int'))
            aa=allsurface150.mcolat
            # add the angles
            allsurface150.mcolat=aa+[180-i for i in aa[::-1]]
            # add the patches
            allsurface150.patches=allsurface150.patches+allsurface150.patches[::-1]
            #apply complete comptonization on every slab
            ii150,qq150 = [],[]
            for ii, pp in enumerate(surfaceaux.patches):
                allsurface150.patches[ii] = Complete_Comptonization(pp,2.1)
                ii150int,qq150int=allsurface150.patches[ii].fluxIQ(self.e_pfield)
                ii150.append(ii150int)
                qq150.append(qq150int)

            #recalculate surface emission model
            ivec_abs, ivec, qvec = [], [], []
            for ee in self.e_pfield:
                pfieldaux.recalculate(ee,allsurface150,gtt=(1-2*2.0/10.0)**0.5)
                ivec.append(pfieldaux.imean)
                qvec.append(pfieldaux.qmean)
            #ivec = np.array(ivec)
            #qvec = np.array(qvec)
            self.i_pfield2d.append(ivec)
            self.iabs_pfield2d.append(totabs*ivec/1.7e4)
            self.q_pfield2d.append(qvec)
            #plt.loglog(self.e_pfield,ivec/self.e_pfield**2)
            if(pfieldaux.theta<=180.0):
                ax1.plot(self.e_pfield,np.array(totabs*ivec/1.7e4))
                ax2.plot(self.e_pfield,np.array(qvec)/np.array(ivec))
        #plt.ylim(1,1e4)
        ax1.set_yscale("log")
        plt.show()
        self.i_pfield2d = np.array(self.i_pfield2d)
        self.iabs_pfield2d = np.array(self.iabs_pfield2d)
        self.q_pfield2d = np.array(self.q_pfield2d)
        return self.e_pfield
    
    def loadRCS(self,dirpath):
        import numpy as np
        import matplotlib.pyplot as plt
        import glob


        # Load a surface model
        ee=np.logspace(-2,2,401)
        allsurface=surface_model().loaddata(glob.glob('../QEDSurface/B14.11T6.48/*.int'))

        # Calculate the spectra
        iilist,qqlist, mcolat = [],[],[]
        for ii, pp in enumerate(allsurface.patches):
            iiint,qqint = pp.fluxIQ(ee)
            iilist.append(iiint)
            qqlist.append(qqint)

        eo = ee/1.3
        eabs,sabs=np.loadtxt('2019-07-17/tbabs.dat',unpack=True)
        sabs=sabs
        ii=np.argsort(eabs)
        eabs=eabs[ii]
        sabs=sabs[ii]
        # absorption cross-section per hydrogen atom in units of 1e-24 cm^2 for our energy bins
        ssabs=np.interp(eo,eabs,sabs)/(eo)**3*1e-24
        # the best-fit hydrogen column is 0.52e22 /cm2
        totabs=np.exp(-0.3e22*ssabs)

        #from Magnetar import partial_res_comptonization
        new_patch150=partial_res_comptonization(allsurface.patches[3],-0.151,150)
        #comp_patch3=allsurface.patches[3]**0.151
        ii150,qq150=new_patch150.fluxIQ(ee)
        new_patch150u=partial_res_comptonization(allsurface.patches[3],0.151,150)
        #comp_patch3=allsurface.patches[3]**0.151
        ii150u,qq150u=new_patch150u.fluxIQ(ee)
        factor=4.5e-6 # 0.346 keV

        plt.semilogx(eo,qq150u/ii150u,label='H  Compt 150 u')
        plt.show()
        plt.loglog(eo,factor*totabs*(ii150u*ee),'k-.',label='H atmo + RCS')
        plt.legend()
        plt.xlabel("Energy [keV]")
        plt.ylabel(r'$E f_E$ [keV$^{-2}$ cm$^{-2}$ s$^{-1}$ keV$^{-1}$')
        #plt.loglog(ee,nout*ee**3)
        plt.xlim(0.5,1e2/1.3)
        plt.ylim(1e-4,1e-1)
        plt.show()
        
        self.e_pfield = ee
        self.th_pfield = np.linspace(0,360,4)
        self.i_pfield2d = np.concatenate([[ii150u]] * 4, axis=0)
        self.iabs_pfield2d = np.concatenate([[factor*totabs*(ii150u*ee)+1e-20]]* 4, axis=0)
        self.q_pfield2d = np.concatenate([[qq150u]] * 4, axis=0)

    def loadRCS2(self,filepath):
        import os
        import fnmatch
        from tkinter import Tcl

        # search pfield files in input directory
        filelist = []
        for path,dirs,files in os.walk(dirpath):
            for file in files:
                if fnmatch.fnmatch(file,'pfield_angle*'):
                    filename = os.path.join(path,file)
                    filelist.append(filename)

        # sort pfield files in accending order according angle in file name
        filelist = Tcl().call('lsort', '-dict', filelist)

        self.th_pfield  = []
        self.i_pfield2d = []
        self.iabs_pfield2d = []
        self.q_pfield2d = []
        fig, ax1 = plt.subplots()
        fig, ax2 = plt.subplots()
        for files in filelist:
            #load pfield files
            pfieldaux = pfield()
            pfieldaux.loaddata(files)
            self.th_pfield.append(pfieldaux.theta)
            #load atmophereric surface model
            surfaceaux = surface_model()
            surfaceaux.loaddata(glob.glob('../QEDSurface/B14.11T6.48/*.int'))
            aa=surfaceaux.mcolat
            #add the angles
            surfaceaux.mcolat=aa+[180-i for i in aa[::-1]]
            #add the patches
            surfaceaux.patches=surfaceaux.patches+surfaceaux.patches[::-1]
            #recalculate surface emission model
            ivec_abs, ivec, qvec = [], [], []
            for ee in self.e_pfield:
                pfieldaux.recalculate(ee,surfaceaux,gtt=(1-2*2.0/10.0)**0.5)
                ivec.append(pfieldaux.imean)
                qvec.append(pfieldaux.qmean)
            #ivec = np.array(ivec)
            #qvec = np.array(qvec)
      

    def lightcurve(self, inclination, field_angle, energy_lc):
        from scipy.interpolate import RectBivariateSpline
        import matplotlib.pylab as plt

        self.energy_lc = energy_lc
        self.qoveri =  self.q_pfield2d/self.i_pfield2d

        i_interp_spline = RectBivariateSpline(self.th_pfield, self.e_pfield, self.iabs_pfield2d)
        log_i_interp_spline = RectBivariateSpline(self.th_pfield, self.e_pfield, np.log10(self.iabs_pfield2d))
        q_interp_spline = RectBivariateSpline(self.th_pfield, self.e_pfield, self.q_pfield2d)
        qoveri_spline   = RectBivariateSpline(self.th_pfield, self.e_pfield, self.qoveri)

        self.inclination = np.radians(inclination)
        self.field_angle = np.radians(field_angle)
        self.phase_array = np.radians(np.linspace(0.36, 360.05730682, 10))

        self.theta_array = np.arccos(
            np.cos(self.field_angle) * np.cos(self.inclination) +
            np.sin(self.field_angle) * np.sin(self.inclination
                                              ) * np.cos(self.phase_array))

        self.eta_array = np.arcsin(np.sin(self.phase_array)*np.sin(self.field_angle)/
                                   np.sin(self.theta_array))

        outputfile = open('ixpesw-ixpeobssim-572fa0092131/ixpeobssim/config/ascii/mag_4U_2hotspots.dat', 'w')
        outputfile.write("\t%.4f\t%.4f\n" % (inclination, field_angle))
        #outputfile.write("%.4f   %.4f\n" % (inclination, field_angle))

        #fig, ax1 = plt.subplots()
        #fig, ax2 = plt.subplots()
        i_array = []
        q_array = []
        i_array2 = []
        qoveri_array =[]
        pa_array = []
        for i in range(len(self.phase_array)-1):
            i_aux   = 10.0**log_i_interp_spline(np.degrees(self.theta_array[i]),self.energy_lc)
            q_aux   =           q_interp_spline(np.degrees(self.theta_array[i]),self.energy_lc)
            qoveri_aux  =         qoveri_spline(np.degrees(self.theta_array[i]),self.energy_lc)
            pa_aux = 0.5*np.arctan2(-np.sin(2.0*self.eta_array[i])*qoveri_aux[0],np.cos(2.0*self.eta_array[i])*qoveri_aux[0])
            pa_aux = np.where(pa_aux<0,np.degrees(np.pi+pa_aux),np.degrees(pa_aux))
            i_array.append(i_aux[0])
            q_array.append(q_aux[0])
            qoveri_array.append(qoveri_aux[0])
            pa_array.append(pa_aux)
            #ax1.plot(self.energy_lc,i_aux[0])
            #ax2.plot(self.energy_lc,qoveri_aux[0])
            #save model to file
            outputfile.write("\t%d\t%.8f\t%.8f\n" % (i+1,self.phase_array[i], self.phase_array[i+1]))
            #outputfile.write("%d %.8f %.8f\n" % (i+1,self.phase_array[i], self.phase_array[i+1]))

            for j in range(len(self.energy_lc)-1):
                e0 = np.log10(self.energy_lc[j])
                e1 = np.log10(self.energy_lc[j+1])
                i0 = i_aux[0][j]/(self.energy_lc[j])**2.0
                pd = np.abs(qoveri_aux[0][j])
                pd = np.where(pd>1.0, 0.99, pd)
                pa = pa_aux[j]
                #print("\t%.6f\t%.6f\t%.7f\t%.6f\t%.4f\n" % (e0,e1,i0,pd,pa))
                outputfile.write("\t%.6f\t%.6f\t%.7e\t%.6f\t%.4f\n" % (e0,e1,i0,pd,pa))

        outputfile.close()
        #ax1.set_yscale("log")
        self.i_array = np.array(i_array)
        self.q_array = np.array(q_array)
        self.qoveri_array = np.array(qoveri_array)
        self.pa_array = np.array(pa_array)

        return 10

def oldpsicalc(mass, radius, b):
    from scipy.integrate import quad
    MoverR = mass / radius
    a2 = (1 - 2 * MoverR) * MoverR * MoverR
    x = b / radius * (1 - 2 * MoverR)**0.5
    # print(x)
    x2 = x * x
    f = (lambda u: 1 / (a2 - (1 - 2 * u) * u * u * x2)**0.5)
    return x * quad(f, 0, MoverR)[0]


# for a given impact parameter b over what angle psi does the light travel [radians]
def psicalc(mass, radius, b):
    from scipy.integrate import quad
    MoverR = mass / radius
    x = b / mass
    # print(x)
    x2 = x * x
    f = (lambda u: (1 - (1 - 2 * u) * u * u * x2)**-0.5)
    return x * quad(f, 0, MoverR)[0]


# Eq 1 of https://arxiv.org/pdf/astro-ph/0201117.pdf
def psiapprox(mass, radius, b):
    MoverR = mass / radius
    rinf = radius / (1 - 2 * MoverR)**0.5
    sinalpha = b / rinf
    cosalpha = (1 - sinalpha**2)**0.5
    cospsi = 1 - (1 - cosalpha) / (1 - 2 * MoverR)
    return np.arccos(cospsi)


# how much is a photon delayed from the surface for a given impact parameter b compared to zero impact parameter (sub-observer point) [same units as mass]
def timecalc(mass, radius, b):
    from scipy.integrate import quad
    MoverR = mass / radius
    x = b / mass
    # print(x)
    x2 = x * x
    f = (
        lambda u: ((1 - (1 - 2 * u) * u * u * x2)**-0.5 - 1) / (1 - 2 * u) / u / u
    )
    return mass * quad(f, 0, MoverR)[0]


# Approximate results inspired by https://arxiv.org/pdf/astro-ph/0201117.pdf
# good to two percent for R>5M (much better if b<0.9 rinf)
def timeapprox(mass, radius, b):
    uu = 1 - np.cos(psiapprox(mass, radius, b))
    return radius * radius * uu / (radius - mass * uu / 3.0)


# calculate both together
def calc_psiandtime(mass, radius, b):
    res1 = psiapprox(mass, radius, b)
    uu = 1 - np.cos(res1)
    return np.degrees(res1), radius * radius * uu / (radius - mass * uu / 3.0)


class rotating_pfield(pfield):
    def setrotationvector(self, rate, inclination, field_angle):
        field_angle = np.abs(field_angle)
        inclination = np.abs(inclination)
        assert(self.theta<=np.abs(field_angle+inclination) and self.theta>=np.abs(field_angle-inclination)),\
            'Magnetic field direction must lie between |inclination=-field_angle| and inclination+field_angle'

        self.rate = rate
        self.inclination = np.radians(inclination)
        self.field_angle = np.radians(field_angle)
        self.phase = 0
        self.phase_delay = []
        self.psi = []
        self.time = []
        self.energy = []
        self.spin_alt = []
        self.spin_colat = []
        self.alpha = []
        self.zeta = []
        self.spin_azim = []
        self.pt = 0
        self.pr = []
        self.pth = []
        self.pph = []
        self.ut = []
        self.ur = []
        self.uth = []
        self.uph = []
        self.delay_mag = []
        self.theta_delay = []
        self.eta_delay = []

        #=== SPIN FRAME ===#
        # calculate the value of psi and light-travel time  at each point
        self.psi, self.time = calc_psiandtime(self.mass, self.radius,
                                              self.data['b'])
        # calculate the phase using spherical trig at the time that light leaves the sub-observer point
        self.phase = np.arccos(
            (np.cos(np.radians(self.theta)) -
             np.cos(self.field_angle) * np.cos(self.inclination)) /
            (np.sin(self.field_angle) * np.sin(self.inclination)))
        # calculate the colatitude of each surface elemente relative to the spin axis
        self.eta = np.arccos(
            (np.cos(self.field_angle) -
             np.cos(self.inclination) * np.cos(np.radians(self.theta))) /
            (np.sin(self.inclination) * np.sin(np.radians(self.theta))))
        self.spin_colat = np.arccos(
            np.cos(self.inclination) * np.cos(np.radians(self.psi)) +
            np.sin(self.inclination) * np.sin(np.radians(self.psi)) * np.
            cos(np.pi - self.data['beta'] - self.eta))
        # zeta angle between n-b plane and n-los plane
        beta = np.where(self.data['beta'] > np.pi, self.data['beta'] - np.pi,
                        self.data['beta'])
        sin_zeta = np.sin(np.radians(self.theta)) * np.sin(
            np.pi - beta) / np.sin(np.radians(self.data['mag_colat']))
        zeta = np.where(sin_zeta > 1.0, np.pi / 2.0, np.arcsin(sin_zeta))
        self.zeta = np.where(
            np.cos(np.radians(self.theta)) <
            np.cos(np.radians(self.data['mag_colat'])) * np.cos(
                np.radians(self.psi)), np.pi - zeta, zeta)
        # alpha angle between n-s plane and n-los plane
        xi = np.where((self.data['beta'] < (np.pi - self.eta)),
                      np.pi - self.data['beta'] - self.eta,
                      self.eta + self.data['beta'] - np.pi)
        xi = np.where((self.data['beta'] > (2.0 * np.pi - self.eta)),
                      2.0 * np.pi + (np.pi - self.data['beta'] - self.eta), xi)
        sin_alpha = np.sin(xi) * np.sin(self.inclination) / np.sin(
            self.spin_colat)
        alpha = np.where(sin_alpha > 1.0, np.pi / 2.0, np.arcsin(sin_alpha))
        self.alpha = np.where(
            np.cos(self.inclination) <
            np.cos(self.spin_colat) * np.cos(np.radians(self.psi)),
            np.pi - alpha, alpha)
        # calculate the azimuth in the rotating frame
        #        aux_azim1 = np.pi - (
        #            self.alpha + np.radians(self.data['phi']) + self.zeta - np.pi)
        #        aux_azim2 = np.pi - (
        #            self.alpha + np.radians(self.data['phi']) - self.zeta)

        aux_azim1 = (
            self.alpha + np.radians(self.data['phi']) + self.zeta - np.pi)
        aux_azim2 = (self.alpha + np.radians(self.data['phi']) - self.zeta)
        self.spin_azim = np.where(self.data['mag_colat'] < 90.0, aux_azim1,
                                  aux_azim2)
        self.spin_alt = np.pi / 2.0 - np.radians(self.data['theta'])

        #=== DOPPLER SHIFT ===#
        # Diagonal components (+,-,-,-) of the metric tensor g_ii
        gtt = 1 - 2 * self.mass / self.radius
        grr = -1.0 / gtt
        gthth = -self.radius**2.0
        gphph = -(self.radius * np.sin(self.spin_colat))**2.0
        # calculate the four velocity u^i= (ut, ur, uth, uph).
        # self.rate  with 1/c factor to convert [hz] to [1/cm]
        c = 3.0e+10
        self.ut = 1.0 / np.sqrt(gtt - gphph * (self.rate / c)**2)
        self.ur = 0
        self.uth = 0
        self.uph = self.ut * (self.rate / c)
        # calculate the photon four momentum p_i=(pt,pr,pth,pph) in the spin frame using p_t=1
        self.pt = 1.0
        #self.pr   = not needed
        #self.pth  = not needed
        self.pph = np.sqrt((-self.pt**2 * gphph / gtt) /
                           ((np.tan(self.spin_alt))**2 *
                            (1 + 1 / (np.tan(self.spin_azim))**2) + 1 /
                            (np.tan(self.spin_azim))**2 + 1))
        # defines where pph is positive or negative on the sphere
        self.pph = np.where((self.data['beta'] < np.pi - self.eta) |
                            (self.data['beta'] > 2.0 * np.pi - self.eta),
                            -self.pph, self.pph)
        #self.pph  = np.where(self.spin_azim>np.pi/2.0,-self.pph,self.pph)
        # calculate the photon energy E = p_i * u^i at the star surface
        self.energy = self.pt * self.ut + self.pph * self.uph

        #=== TIME DELAY ===#
        # because of the time delay  and the star rotation, for every b>0 there is a slight change in the direction of the
        # mag axis,  which should be recomputed from phase_delay. One option is to assign self.phase to b=0, then for b>=0
        self.phase_delay = self.phase - self.rate * self.time / c
        # compute the theta_delay angle between the magnetic axis and LOS
        cos_theta_delay = np.cos(self.field_angle) * np.cos(
            self.inclination) + np.sin(self.field_angle) * np.sin(
                self.inclination) * np.cos(self.phase_delay)
        self.theta_delay = np.arccos(cos_theta_delay)
        # same as the lines before but for the eta angle
        cos_eta_delay = (
            np.cos(self.field_angle) -
            np.cos(self.inclination) * np.cos(self.theta_delay)) / (
                np.sin(self.inclination) * np.sin(self.theta_delay))
        self.eta_delay = np.arccos(cos_eta_delay)
        # cartesians components of the spin axis
        sx = np.sin(self.inclination) * np.cos(self.eta)
        sy = np.sin(self.inclination) * np.sin(self.eta)
        sz = np.cos(self.inclination)
        # the components of the spin axis and delayed magnetic axis satisfy the system
        #     (1) bx*sx + by*sy + bz*sz = cos(field_angle)
        #     (2) by**2 + by**2 = (sin(theta_delay))**2
        #     (3) bz = cos(theta_delay)
        # this yields an equation a1*bx**2 + a2*bx + c = 0, whose coeficients are
        a1 = sx**2.0 + sy**2.0
        a2 = -2.0 * sx * (
            np.cos(self.field_angle) - np.cos(self.theta_delay) * sz)
        a3 = (np.cos(self.field_angle) - np.cos(self.theta_delay) * sz
              )**2.0 - (np.sin(self.theta_delay) * sy)**2.0
        # calculate the solution of the polynomial equation and the components (bx,by,bz) of the delayed magnetic axis
        delta1 = a2**2.0 - 4.0 * a1 * a3
        self.bx = np.where(delta1 < 0, -a2 / (2.0 * a1),
                           (-a2 + np.sqrt(delta1)) / (2.0 * a1))
        delta2 = (np.sin(self.theta_delay))**2.0 - self.bx**2.0
        self.by = np.where(delta2 < 0, 0, -np.sqrt(delta2))
        self.bz = np.cos(self.theta_delay)
        # cartesians components of the surface element
        nx = np.sin(np.radians(self.psi)) * np.cos(np.pi - self.data['beta'])
        ny = np.sin(np.radians(self.psi)) * np.sin(np.pi - self.data['beta'])
        nz = np.cos(np.radians(self.psi))
        # calculate the magnetic colatitude using cos(mag_col) = b*n
        self.mag_delay = np.arccos(self.bx * nx + self.by * ny + self.bz * nz)
##############################################################################################################################
        #calculate  zeta delay. I need theta delay and mag_delay
        # zeta angle between n-b plane and n-los plane
#        beta = np.where(self.data['beta'] > np.pi, self.data['beta'] - np.pi,
#                        self.data['beta'])
#        beta = np.where(self.data['beta'] > np.pi, self.data['beta'],
#                        self.data['beta'])
#        sin_zeta_delay = np.sin(self.theta_delay) * np.sin(
#            np.pi - beta) / np.sin(self.mag_delay)
#        zeta = np.where(sin_zeta_delay > 1.0, np.pi / 2.0, np.arcsin(sin_zeta_delay))
#        self.zeta_delay = np.where(np.cos(self.theta_delay)<np.cos(self.mag_delay) * np.cos(np.radians(self.psi)),
#                                   np.pi - zeta, zeta)
        cos_zeta_delay  = (np.cos(self.theta_delay)-np.cos(self.mag_delay)*np.cos(np.radians(self.psi)))/(np.sin(self.mag_delay)*np.sin(np.radians(self.psi)))
        self.zeta_delay = np.arccos(cos_zeta_delay)

##############################################################################################################################
        #=== LIGHT ABERRATION ===#
        zen = np.radians(self.data['theta'])
        # photon azimuth angle in the spin frame. Notice that self.spin_azim is between [0,pi].
        # The azimuth angle is redefined below between [-pi,pi], consistent with the
        # range of the output of the np.arctan2 function used later.
        azim = np.where((self.data['beta'] < np.pi - self.eta) |
                        (self.data['beta'] > 2.0 * np.pi - self.eta),
                        self.spin_azim, -self.spin_azim)
        self.spin_azim2 = azim
        # calculate v/c of the surface element. A correction due to gravitational redshit migth be needed
        self.beta_spin = self.rate * self.radius * np.sin(self.spin_colat) / c
        # aberration in the spin frame
        self.spin_zen_aberr = np.arccos(
            np.cos(zen) * np.sqrt(1 - self.beta_spin**2.0) /
            (1 - np.sin(zen) * np.sin(azim) * self.beta_spin))
        # arctan2 output between [-pi,pi]. the abs function narrow it to [0,pi]
        #        self.spin_azim_aberr = np.abs(np.arctan2(
        #            np.sin(zen) * np.sin(azim) - self.beta_spin,
        #            np.sin(zen) * np.cos(azim) * np.sqrt(1 - self.beta_spin**2.0)))

        self.spin_azim_aberr = (np.arctan2(
            np.sin(zen) * np.sin(azim) - self.beta_spin,
            np.sin(zen) * np.cos(azim) * np.sqrt(1 - self.beta_spin**2.0)))

        # calculate aberration back to the magnetic frame for different regions over the sphere
        mag_azim_aberr1 =  self.alpha + self.zeta + self.spin_azim_aberr
        mag_azim_aberr2 = -self.alpha + self.zeta + self.spin_azim_aberr
        mag_azim_aberr3 = -self.alpha - self.zeta + self.spin_azim_aberr
        mag_azim_aberr4 =  self.alpha - self.zeta + self.spin_azim_aberr
        mag_azim_aberr5 = -self.alpha + self.zeta + self.spin_azim_aberr + 2.0 * np.pi
        mag_azim_aberr6 = -self.alpha - self.zeta + self.spin_azim_aberr
        mag_azim_aberr7 = -self.alpha + self.zeta + self.spin_azim_aberr
        mag_azim_aberr8 = -self.alpha - self.zeta + self.spin_azim_aberr + 2.0 * np.pi

        #regions over the sphere
        region1 = (self.data['beta'] > np.pi - self.eta) & (self.data['beta'] < 2.0 * np.pi - self.eta)
        region2 =  self.data['beta'] < np.pi
        region3 =  self.data['beta'] < np.pi - self.eta
        region4 =  self.spin_azim_aberr>0
        region5 =  np.radians(self.psi)>np.pi/2.0

        mag_azim_aberr = np.where(region1,
                                  np.where(region2,mag_azim_aberr1,mag_azim_aberr4),
                                  np.where(region3,
                                           np.where(region4,
                                                    mag_azim_aberr2,
                                                    np.where(region5,mag_azim_aberr7,mag_azim_aberr5)),
                                           np.where(region4,
                                                    mag_azim_aberr3,
                                                    np.where(region5,mag_azim_aberr8,mag_azim_aberr6))))

        mag_azim_aberr = np.where(np.abs(mag_azim_aberr) > np.pi,
                                  np.where(self.data['beta'] < np.pi, mag_azim_aberr - 2.0 * np.pi, 2.0 * np.pi + mag_azim_aberr),
                                  mag_azim_aberr)

        self.mag_azim_aberr = mag_azim_aberr
        return self

#       condition4 = (self.data['beta'] < np.pi - self.eta) | (self.data['beta'] > 2.0 * np.pi - self.eta)
#       condition5 =  self.data['beta'] < np.pi - self.eta
#       condition6 =  self.data['beta'] > 2.0 * np.pi - self.eta
#        mag_azim_aberr = np.where(condition1,
#                                  np.where(condition2,mag_azim_aberr1,mag_azim_aberr4),
#                                  np.where(condition3,mag_azim_aberr2,mag_azim_aberr3))

#        mag_azim_aberr = np.where(condition4,
#                                  np.where(self.spin_azim_aberr < 0,
#                                           np.where(self.data['beta'] < np.pi,mag_azim_aberr5,mag_azim_aberr6),
#                                           mag_azim_aberr),
#                                  mag_azim_aberr)
#        mag_azim_aberr = np.where(condition5,
#                                  np.where(self.spin_azim_aberr < 0,
#                                           np.where(np.radians(self.psi) > np.pi / 2.0,mag_azim_aberr7,mag_azim_aberr),
#                                           mag_azim_aberr),
#                                  mag_azim_aberr)
#        mag_azim_aberr = np.where(condition6,
#                                  np.where(self.spin_azim_aberr < 0,
#                                           np.where(np.radians(self.psi) > np.pi / 2.0,mag_azim_aberr8,mag_azim_aberr), 
#                                           mag_azim_aberr),
#                                  mag_azim_aberr)
