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

from ximpol import XIMPOL_CONFIG
from ximpol.core.rand import xUnivariateGenerator
from ximpol.core.spline import xInterpolatedUnivariateSpline
from ximpol.srcmodel.roi import xPeriodicPointSource, xEphemeris, xROIModel
from ximpol.srcmodel.spectrum import power_law


def _full_path(file_name):
    """Convenience function to retrieve the relevant files.
    """
    return os.path.join(XIMPOL_CONFIG, 'ascii', file_name)


# geometry of dipole
alpha=numpy.radians(30)
beta=numpy.radians(15)
t0=0.1
# Mind that you have to wrap this into a function to be used.
# This could be done in a more general way with some sort of library function.
def polarization_degree(E, t, ra, dec):
    return 1.0+0*t


#
def inclination(t):
    phi_phase=numpy.radians((t-t0)*360)
    return numpy.arccos(numpy.cos(alpha)*numpy.cos(alpha-beta)+
                        numpy.sin(alpha)*numpy.sin(alpha-beta)*numpy.cos(phi_phase))


pdegree_ang=numpy.radians(numpy.linspace(0,90,19))
pdegree_data1=2*numpy.array([0,0.032765,0.0490032,0.065200,0.070770,0.086773,0.1025,0.116778,0.130228,0.142785,0.154202,0.164246,0.172673,0.180627,0.187338,0.192776,0.196445,0.198123,0.198696])

pdegree_data2=numpy.array([0,0.14017,0.485033,0.787875,0.862695,0.919405,0.934439,0.95239,0.965822,0.969983,0.974817,0.979925,0.980247,0.984679,0.9794,0.988191,0.984482,0.987218,0.985947])

pdegree_funk1=(scipy.interpolate.interp1d(pdegree_ang,pdegree_data1,kind='cubic'))
pdegree_funk2=(scipy.interpolate.interp1d(pdegree_ang,pdegree_data2,kind='cubic'))


def polarization_degree_noqed(E, t, ra, dec):
    efactor=numpy.tanh(2*(E-4))
    efactor=(efactor+1)/2
    return pdegree_funk1(inclination(t))*(1-efactor)+pdegree_funk2(inclination(t))*(efactor)

# And, again, this needs to be wrapped into a function.
def polarization_angle(E, t, ra, dec):
    phi_phase=numpy.radians((t-t0)*360)
    return numpy.arctan(numpy.sin(alpha)*numpy.sin(phi_phase)/
                        (numpy.sin(alpha+beta)*numpy.cos(alpha)-
                         numpy.cos(alpha+beta)*numpy.sin(alpha)*
                         numpy.cos(phi_phase)))+1.5707963268




# Build the actual energy spectrum.
# energy_spectrum = power_law(pl_normalization_spline, pl_index_spline)
fdata=[]
for fname in ("model02.dat","model24.dat","model46.dat","model68.dat","model80.dat"):
    en, f = numpy.loadtxt(_full_path(fname),usecols=(0,3),unpack=True,skiprows=3)
    fdata.append(f)

phdata=numpy.array([-0.1,0.1,0.3,0.5,0.7,0.9,1.1])
fdata.insert(0,fdata[-1])
fdata.append(fdata[1])
energy_spectrum=numpy.vectorize(scipy.interpolate.interp2d(en,phdata,fdata,kind='quintic'))

# Build the flux as a function of the phase.
fmt = dict(xname='Pulsar phase', yname='Flux ( 2 - 10 keV )',
           yunits='keV cm$^{-2}$ s$^{-1}$')


enbin=numpy.linspace(2,10,100)
esum=phdata*0
for ii,pp in enumerate(phdata):
    esum[ii]=numpy.trapz(energy_spectrum(enbin,pp),x=enbin)
pl_normalization_spline = xInterpolatedUnivariateSpline(phdata, esum, k=3, **fmt)


fmt = dict(xname='Pulsar phase', yname='Polarization angle [rad]')
_phi=numpy.linspace(0,1,23)
_pol_angle = polarization_angle(0,_phi,0,0)
pol_angle_spline = xInterpolatedUnivariateSpline(_phi, _pol_angle, k=1, **fmt)


_pol_degree = polarization_degree(0,_phi,0,0)
fmt = dict(xname='Pulsar phase', yname='Polarization degree')
pol_degree_spline = xInterpolatedUnivariateSpline(_phi, _pol_degree, k=1, **fmt)

_pol_degree_noqed = polarization_degree_noqed(0,_phi,0,0)
pol_degree_spline_noqed = xInterpolatedUnivariateSpline(_phi,
                                                        _pol_degree_noqed,
                                                        k=1, **fmt)

ROI_MODEL = xROIModel(26.59342,  61.75078)
fouru_ephemeris = xEphemeris(0., 0.115088121,  -2.64E-14, 0)
fouru_pulsar = xPeriodicPointSource('4U 0142+61', ROI_MODEL.ra, ROI_MODEL.dec,
                                   energy_spectrum, polarization_degree_noqed,
                                   polarization_angle, fouru_ephemeris)
ROI_MODEL.add_source(fouru_pulsar)


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
