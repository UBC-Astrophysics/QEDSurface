import numpy as np

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

from scipy.ndimage import map_coordinates

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
            res = map_coordinates(
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
            resplus = map_coordinates(
                datacube,
                self._calcindex(dataarrayplus),
                order=3,
                mode='nearest')
            resplus = np.where(resplus > 0, resplus, 0)
            resminus = map_coordinates(
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

