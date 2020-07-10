import numpy as np
from Magnetar.utils import atmosphere

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

def dipole_model(cval,tpole,bpole,*args,**kwargs):
    return cval(tpole,bpole,0.0,*args,**kwargs)    