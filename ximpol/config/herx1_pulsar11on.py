#!/usr/bin/env python
#
# Copyright (C) 2016, the ximpol team.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU GengReral Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

"""References:
Slowikowska A. et al., MNRAS 397, 103-23

Valid range (MJD)       : 52944--52975
Epoch, t0 (MJD)         : 52960.000000296
nu0 (Hz)                : 29.8003951530036
nudot(10^-10 Hz s^-1)   : -3.73414
nudddot (10^-20 Hz s^-2): 1.18
"""

import numpy
import scipy.signal
import scipy.interpolate
from scipy.optimize import curve_fit
import os

try: 
    from ximpol import XIMPOL_CONFIG
    from ximpol.core.rand import xUnivariateGenerator
    from ximpol.core.spline import xInterpolatedUnivariateSpline
    from ximpol.srcmodel.roi import xPeriodicPointSource, xEphemeris, xROIModel
    from ximpol.srcmodel.spectrum import power_law
    ximpol_loaded=True
except ImportError:
    ximpol_loaded=False

def _full_path(file_name):
    """Convenience function to retrieve the relevant files.
    """
    if ximpol_loaded:
        return os.path.join(XIMPOL_CONFIG, 'ascii', file_name)
    else:
        return os.path.join('..', 'ascii', file_name)

# Grab all the relevant files.

# radius angle energy phot_flux (x-o)/(x+o)
# rad, inclination, energy, phot_flux, ratio = numpy.loadtxt(_full_path("11qedon.txt"),usecols=range(5),
#unpack=True)

dtype=[('incl',float),('ener',float),('flux',float),('ratio',float)]
arr=numpy.loadtxt(_full_path("11qedon.txt"),usecols=(1,2,3,4),dtype=dtype)
arr=numpy.sort(arr,order=('ener','incl'))
inclination=arr['incl']
energy=arr['ener']
phot_flux=arr['flux']
ratio=arr['ratio']
inclination=numpy.radians(inclination)
incllist=numpy.unique(inclination)
enerlist=numpy.unique(energy)
energy_spectrum_inclination=numpy.vectorize(scipy.interpolate.RectBivariateSpline(enerlist,incllist,
                                                                  (phot_flux*energy).reshape((len(enerlist),len(incllist))),kx=3,ky=3))
ratio_inclination=numpy.vectorize(scipy.interpolate.RectBivariateSpline(enerlist,incllist,
                                                                  (ratio).reshape((len(enerlist),len(incllist))),kx=3,ky=3))
# energy_spectrum_inclination=numpy.vectorize(scipy.interpolate.interp2d(energy,inclination,phot_flux*energy,kind='linear'))
# ratio_inclination=numpy.vectorize(scipy.interpolate.interp2d(energy,inclination,ratio,kind='linear'))

# geometry of dipole
alpha=numpy.radians(85)
beta=numpy.radians(5)
# Mind that you have to wrap this into a function to be used.
# This could be done in a more general way with some sort of library function.

#
def inclination(t):
    phi_phase=numpy.radians(t*360)
    return numpy.arccos(numpy.abs(numpy.cos(alpha)*numpy.cos(alpha-beta)+
                        numpy.sin(alpha)*numpy.sin(alpha-beta)*numpy.cos(phi_phase)))


def polarization_degree(E, t, ra, dec):
    return numpy.abs(ratio_inclination(E,inclination(t)))


# And, again, this needs to be wrapped into a function.
def polarization_angle(E, t, ra, dec):
    phi_phase=numpy.radians(t*360)
    rat=ratio_inclination(E,inclination(t))
    ang=numpy.arctan(numpy.sin(alpha)*numpy.sin(phi_phase)/
                        (numpy.sin(alpha+beta)*numpy.cos(alpha)-
                         numpy.cos(alpha+beta)*numpy.sin(alpha)*
                         numpy.cos(phi_phase)))
    return numpy.where(rat>0,ang+1.5707963268,ang)

def energy_spectrum(E,t):
    return energy_spectrum_inclination(E,inclination(t))

if ximpol_loaded:
    
    # Build the PL normalization as a function of the phase.
    fmt = dict(xname='Pulsar phase', yname='Flux',
               yunits='keV cm$^{-2}$ s$^{-1}$')

    _phi=numpy.linspace(0,1,23)
    enbin=numpy.linspace(2,10,100)
    esum=_phi*0
    for ii,pp in enumerate(_phi):
        esum[ii]=numpy.trapz(energy_spectrum(enbin,pp),x=enbin)

    pl_normalization_spline = xInterpolatedUnivariateSpline(_phi, esum, k=3, **fmt)

    # Build the polarization angle as a function of the phase.
    _pol_angle = polarization_angle(0,_phi,0,0)
    fmt = dict(xname='Pulsar phase', yname='Polarization angle [rad]')
    pol_angle_spline = xInterpolatedUnivariateSpline(_phi, _pol_angle, k=1, **fmt)


    # Build the polarization degree as a function of the phase.
    _pol_degree = polarization_degree(0,_phi,0,0)
    fmt = dict(xname='Pulsar phase', yname='Polarization degree')
    pol_degree_spline = xInterpolatedUnivariateSpline(_phi, _pol_degree, k=1, **fmt)


    ROI_MODEL = xROIModel(254.457625,  35.3423888889)
    herx1_ephemeris = xEphemeris(0., 0.8064516129,  0, 0)
    herx1_pulsar = xPeriodicPointSource('Her X-1', ROI_MODEL.ra, ROI_MODEL.dec,
                                        energy_spectrum, polarization_degree,
                                        polarization_angle, herx1_ephemeris)
    ROI_MODEL.add_source(herx1_pulsar)


if __name__ == '__main__':
    print(ROI_MODEL)
    from ximpol.utils.matplotlib_ import pyplot as plt
    plt.figure()
    pl_index_spline.plot(show=False)
    plt.figure()
    pl_normalization_spline.plot(show=False)
    plt.figure()
    pol_angle_spline.plot(show=False)
    plt.figure()
    pol_degree_spline.plot(show=False)
    plt.show()
