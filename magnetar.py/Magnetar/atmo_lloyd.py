#from numba import jit
import numpy as np
from .utils import atmosphere
from scipy.ndimage import map_coordinates
try:
    from .fast_interp import interp3d, interp2d
except:
    pass

def rreplace(s, old, new, occurrence):
    li = s.rsplit(old, occurrence)
    return new.join(li)

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

    def __init__(self,file=None):
        self.t = []
        self.x = []
        self.o = []
        self.file = ""
        self.mag_inclination = 0
        if file is not None:
            self.loaddata(file)
    def loaddata(self, file):
        # header = file.rsplit('.', 1)[0]
        enerkeV = np.loadtxt(rreplace(file,'.int','.gnu',1), usecols=(3), unpack=True)
        self.t, self.x, self.o = np.loadtxt(file, usecols=(2, 3, 4), unpack=True)
        theta, phi = np.loadtxt(rreplace(file,'.int','.trj',1), usecols=(4, 6), unpack=True)
        self.xxarray = []
        self.yyarray = []
        self.shape = []
        self.file=file
        for dd in (theta, phi, enerkeV):
            xx = np.unique(dd)
            self.xxarray.append(xx)
            self.yyarray.append(np.arange(len(xx)))
            self.shape.append(len(xx))
        self.t = np.reshape(self.t, self.shape)
        self.x = np.reshape(self.x, self.shape)
        self.o = np.reshape(self.o, self.shape)

        try:
            if self.shape[1]==1:
                tempshape=[self.shape[0],self.shape[2]]
                tfunk= interp2d([0,0],tempshape,[1,1],np.reshape(self.t,tempshape))
                xfunk= interp2d([0,0],tempshape,[1,1],np.reshape(self.x,tempshape))
                ofunk= interp2d([0,0],tempshape,[1,1],np.reshape(self.o,tempshape))
                self.tfunk= (lambda a,b,c : tfunk(a,c))
                self.xfunk= (lambda a,b,c : xfunk(a,c))
                self.ofunk= (lambda a,b,c : ofunk(a,c))
            else:
                self.tfunk= interp3d([0,0,0],self.shape,[1,1,1],self.t)
                self.xfunk= interp3d([0,0,0],self.shape,[1,1,1],self.x)
                self.ofunk= interp3d([0,0,0],self.shape,[1,1,1],self.o)
        except:
            self.tfunk = None
        return self
    def __str__(self):
        outstring='''#
# class atmo_lloyd 
#
# file                 %s
# magnetic inclination %g
''' % (self.file,self.mag_inclination)
        return outstring+atmosphere.__str__(self)
    # @jit(nopython=True,parallel=True)
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
            if self.tfunk is None:
                res = map_coordinates(
                    datacube, self._calcindex(dataarray), order=3, mode='nearest')
                return np.where(res > 0, res, 0)
            else:
                ci=self._calcindex(dataarray)
                #print(ci[2])
                if datacube is self.o:
                    res=self.ofunk(ci[0],ci[1],ci[2])
                elif datacube is self.x:
                    res=self.xfunk(ci[0],ci[1],ci[2])
                else:
                    res=self.tfunk(ci[0],ci[1],ci[2])
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
            resplus = _intensity(dataarrayplus,datacube)
            resminus = _intensity(dataarrayminus,datacube)
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

