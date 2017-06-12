import numpy as np
import sys

from scipy.interpolate import interp1d

ang, lpol = np.loadtxt(sys.argv[1],unpack=True,delimiter=',')

enindex=(np.floor(ang/100)).astype(int)
angtrue=ang-enindex*100

lpolsmooth=[]
aout=np.linspace(0,90,91)
for ii in range(6):
    x=angtrue[enindex==ii]
    x=np.insert(x,0,0)
    x=np.append(x,90)
    y=lpol[enindex==ii]
    y=np.insert(y,0,0)
    y=np.append(y,y[-1])
    f=interp1d(x,y)
    lpolsmooth.append(f(aout))
    
np.savetxt(sys.argv[2],np.transpose(np.vstack((aout,lpolsmooth))),fmt='%10.6f',header='  Angle      Pi(1keV)  Pi(2keV)   Pi(5keV)  Pi(10keV)  Pi(15keV)  Pi(30keV)')

