import numpy as np
from Magnetar.utils import atmosphere


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
        self.file=file
        self.Eval, self.i0, self.i1, self.i2, self.q0, self.q1, self.q2 = np.loadtxt(
            file, unpack=True, usecols=range(7))
    def __str__(self):
        outstring='''#
# class scatmo_caiazzo 
#
# file                 %s
# magnetic inclination %g
''' % (self.file,self.mag_inclination)
        return outstring+atmosphere.__str__(self)

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

