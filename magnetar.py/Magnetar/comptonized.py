import numpy as np
from Magnetar.utils import *

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
