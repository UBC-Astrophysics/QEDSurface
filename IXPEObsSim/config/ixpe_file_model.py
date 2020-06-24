#!/usr/bin/env python
#
# Copyright (C) 2016--2019, the ixpeobssim team.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
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

from __future__ import print_function, division

import numpy
import os

### BEGIN KEY DEFINITIONS

## Source position and spin rate
source_name = '4U 0142+61'
ra, dec = 26.5933625, 61.7508861111
nu0 = 0.11509211554
nudot = -2.6783806e-14
nuddot = 0.

## File with the model
# should be in the config/ascii directory
# columns should be named in the first row
# EnergykeV Phirad EnergykeV QI I
# the rows and columns can be in any order
filename="Double_Blackbody.txt"

# are the intensities or fluxes in the file in energy units? (IXPEObsSim wants photon units)
# this is done after other renormalizations
intensity_energy_units=True

## Geometry of dipole
# alpha is the angle between the line of sight and rotation axis
alpha=numpy.radians(45)
# beta is the angle between the magnetic axis and line of sight
beta=numpy.radians(85)

## Renormalize the phase-averaged flux
# different renormalizations for the phase averaged flux
normflux=1e-10 # total flux from 2-8 keV in (counts or erg)/cm2/s
#
#normflux=enerlist**(-2)*1e-2 # normalize by an array
#
#normflux=(lambda x: 1e-3*x**-2) # normalize by a function of energy
#
#normflux='renorm.txt' # normalize using data in the file
#  the first row should name the columns including EnergykeV and I

## Value of NH in per cm2 (either give a value or a value and a file from the config/ascii directory
# ROIModel can also include the absorption
# NH=1e24 
# NH='1e22;tbabs.dat'
# the first row should name the columns including Energy and sigma
#     sigma is the cross section times (E/keV)^2 / (1e-24 cm^2)
#     ssabs=numpy.interp(enerlist,abarr['Energy'],abarr['sigma']/(enerlist)**3*1e-24


# final band renorm (you can renormalize the phase-average flux in the band 2-8 keV after the absorption
# finalnorm=1e-2 # total flux from 2-8 keV in (counts or erg)/cm2/s after absorption

### END KEY DEFINITIONS



from scipy.interpolate import RectBivariateSpline
try:
    from ixpeobssim import IXPEOBSSIM_CONFIG_ASCII
    from ixpeobssim.config import file_path_to_model_name
    from ixpeobssim.core.spline import xInterpolatedBivariateSpline
    from ixpeobssim.srcmodel.roi import xPeriodicPointSource, xEphemeris, xROIModel
    from ixpeobssim.utils.matplotlib_ import plt, setup_gca, last_line_color
    ixpe_loaded=True
except ImportError:
    ixpe_loaded=False

def _full_path(file_name):
    """Convenience function to retrieve the relevant files.
    """
    if ixpe_loaded:
        return os.path.join(IXPEOBSSIM_CONFIG_ASCII, file_name)
    else:
        return os.path.join('.', 'ascii', file_name)

if ixpe_loaded:
    __model__ = file_path_to_model_name(__file__)
else:
    __model__ = 'ixpe_file_model'


#dtype=[('incl',float),('ener',float),('flux',float),('ratio',float)]

arr=numpy.genfromtxt(_full_path(filename),names=True)
# sort in the order that we need
arr=numpy.sort(arr,order=('EnergykeV','Phirad'))

inclination=arr['Phirad']
energy=arr['EnergykeV']
ratio=arr['QI']
incllist=numpy.unique(inclination)
enerlist=numpy.unique(energy)
fluxmod=arr['I'].reshape((len(enerlist),len(incllist)))

#
# if the intensity is an energy flux, then convert to make a photon flux
#

if intensity_energy_units:
    fluxmod=numpy.transpose(numpy.transpose(fluxmod)/(enerlist*1.60217662e-9)) # 1.60217662e-9 erg = 1 keV

# build the spectrum and polarization as a function of energy and inclination 
energy_spectrum_inclination=numpy.vectorize(RectBivariateSpline(enerlist,incllist,
                                                                  fluxmod,kx=3,ky=3))
ratio_inclination=numpy.vectorize(RectBivariateSpline(enerlist,incllist,
                                                                  (ratio).reshape((len(enerlist),len(incllist))),kx=3,ky=3))

# energy_spectrum_inclination=numpy.vectorize(scipy.interpolate.interp2d(energy,inclination,flux*energy,kind='linear'))
# ratio_inclination=numpy.vectorize(scipy.interpolate.interp2d(energy,inclination,ratio,kind='linear'))

# inclination as a function of phase
def inclination(t):
    phi_phase=numpy.radians(t*360)
    if (incllist[-1] > numpy.pi/2):
        return numpy.arccos(numpy.cos(alpha)*numpy.cos(alpha-beta)+
                                      numpy.sin(alpha)*numpy.sin(alpha-beta)*numpy.cos(phi_phase))
    else:
        return numpy.arccos(numpy.abs(numpy.cos(alpha)*numpy.cos(alpha-beta)+
                                      numpy.sin(alpha)*numpy.sin(alpha-beta)*numpy.cos(phi_phase)))

# polarization degree as a function of energy and phase
def pol_deg(E, t, ra=None, dec=None):
    return numpy.abs(ratio_inclination(E,inclination(t)))


# polarization angle as a function of energy and phase
def pol_ang(E, t, ra=None, dec=None):
    phi_phase=numpy.radians(t*360)
    rat=ratio_inclination(E,inclination(t))
    ang=numpy.arctan(numpy.sin(alpha)*numpy.sin(phi_phase)/
                        (numpy.sin(alpha+beta)*numpy.cos(alpha)-
                         numpy.cos(alpha+beta)*numpy.sin(alpha)*
                         numpy.cos(phi_phase)))
    return numpy.where(rat>0,ang+1.5707963268,ang)

# energy spectrum as a function of energy and phase
def spec(E,t):
    return energy_spectrum_inclination(E,inclination(t))

# set up phase bins
phase=numpy.linspace(0,1,101)
#    print(phase,inclination(phase))
#    ee,tt=numpy.meshgrid(enerlist,phase)
tt,ee=numpy.meshgrid(phase,enerlist)
flux=spec(ee,tt)
meanflux=numpy.mean(flux,axis=-1)

# check if we have to renormalize
try:
    normflux
except NameError:
    normflux = None

# perform the renormalization
if normflux is not None:
    if type(normflux) is float:
        from scipy.integrate import simps
        # renormalize the total flux over the band 2-8 keV
        ok=(enerlist>2) & (enerlist<8)
        flux=flux/simps(meanflux[ok],enerlist[ok])*normflux
    elif type(normflux) is str:
        narr=numpy.genfromtxt(_full_path(normflux),names=True)
        narr=numpy.sort(narr,order=('EnergykeV'))
        flux=numpy.transpose((numpy.interp(enerlist,narr['EnergykeV'],narr['I'])/meanflux)*numpy.transpose(flux))
    elif callable(normflux):
        # assume norm flux is a function of energy 
        flux=numpy.transpose((normflux(enerlist)/meanflux)*numpy.transpose(flux))
    else:
        # assume norm flux is an array of fluxes at the same energies as enerlist
        flux=numpy.transpose((normflux/meanflux)*numpy.transpose(flux))


# check if there is an NH value
try:
    NH
except NameError:
    NH = None

if NH is not None:
    if type(NH) is float:
        abfilename='tbabs.dat'
    elif type(NH) is str:
        aa=NH.split(';')
        abfilename=aa[1]
        NH=float(aa[0])
    else: 
        raise ValueError('The value of NH should be a float or a string with value;filename.')
    abarr=numpy.genfromtxt(_full_path(abfilename),names=True)
    abarr=numpy.sort(abarr,order=('Energy'))
    ssabs=numpy.interp(enerlist,abarr['Energy'],abarr['sigma'])/(enerlist)**3*1e-24
    flux=numpy.transpose(numpy.exp(-NH*ssabs)*numpy.transpose(flux))

# check if we have to do a final normalization
try:
    finalnorm
except NameError:
    finalnorm = None

if finalnorm is not None:
    if type(finalnorm) is float:
        meanflux=numpy.mean(flux,axis=-1)
        from scipy.integrate import simps
        # renormalize the total flux over the band 2-8 keV
        ok=(enerlist>2) & (enerlist<8)
        flux=flux/simps(meanflux[ok],enerlist[ok])*finalnorm
    else:
        raise ValueError('The value of finalnorm should be a float.')

#
# if the intensity is an energy flux, then renormalize to make a photon flux
# KEY LINE: intensity units
#
if intensity_energy_units:
    flux=numpy.transpose(numpy.transpose(flux)/(enerlist*1.60217662e-9)) # 1.60217662e-9 erg = 1 keV
        
if ixpe_loaded:

#    print(phase,inclination(phase))
#    ee,tt=numpy.meshgrid(enerlist,phase)

    # setup the splines which should be faster than doing the geometry everytime
    # also they are useful for the plot and include the spectral renormalization
    fmt = dict(xlabel='Energy [keV]', ylabel='Phase',
               zlabel='Flux [cm$^{-2}$ s$^{-1}$ keV$^{-1}$]')

    spec_spline = xInterpolatedBivariateSpline(enerlist, phase, flux, kx=3, ky=3, **fmt)

    fmt = dict(xlabel='Energy [keV]', ylabel='Phase', zlabel='Polarization degree')

    pol_deg_data=pol_deg(ee,tt)
    pol_deg_spline = xInterpolatedBivariateSpline(enerlist, phase, pol_deg_data,kx=3, ky=3, **fmt)

    fmt = dict(xlabel='Energy [keV]', ylabel='Phase',zlabel='Polarization angle [deg]')

    pol_ang_data=numpy.degrees(pol_ang(ee,tt))
    pol_ang_spline = xInterpolatedBivariateSpline(enerlist, phase, pol_ang_data,kx=3, ky=3, **fmt)

    # Move on to the actual ROI model.
    ROI_MODEL = xROIModel(ra, dec)
    ephem = xEphemeris(0., nu0, nudot, nuddot)
    src = xPeriodicPointSource(source_name, ra, dec, spec_spline, pol_deg_spline, pol_ang_spline, ephem)
    ROI_MODEL.add_source(src)

def display_spectrum(emin=1.1, emax=12., phase_indices=[10, 40, 60, 80]):
    """Display the energy spectrum.
    """
    # Full, 2-d energy spectrum.
    plt.figure('%s spectrum' % __model__)
    spec_spline.plot()

    # Slices of the energy spectrum at different pulse-phase values.
    plt.figure('%s spectrum phase slices' % __model__)
    for i in phase_indices:
        _phase = phase[i]
        slice_ = spec_spline.hslice(_phase, k=3)
        slice_.plot(label='Pulse phase = %.2f' % _phase)
        plt.plot(enerlist, flux[:,i], 'o', color=last_line_color())
    setup_gca(xmin=emin, xmax=emax, logx=True, logy=True, grids=True,
              ymin=1.e-5, ymax=2.e-2, legend=True)


def display_pol_deg(emin=1.1, emax=12., phase_indices=[10, 40, 60, 80],
                    energy_indices=[4, 7, 10, 13]):
    """Display the polarization degree.
    """
    # Polarization degree 2-d plot.
    plt.figure('%s polarization degree' % __model__)
    pol_deg_spline.plot()

    # Slices of the polarization degree at different pulse-phase values.
    plt.figure('%s polarization degree phase slices' % __model__)
    for i in phase_indices:
        _phase = phase[i]
        slice_ = pol_deg_spline.hslice(_phase, k=3)
        slice_.plot(label='Pulse phase = %.2f' % _phase)
        plt.plot(enerlist, pol_deg_data[:,i], 'o', color=last_line_color())
    setup_gca(xmin=emin, xmax=emax, legend=True)

    # Slices of the polarization degree at different energies.
    plt.figure('%s polarization degree energy slices' % __model__)
    for i in energy_indices:
        _energy = enerlist[i]
        slice_ = pol_deg_spline.vslice(_energy, k=3)
        slice_.plot(label='Energy = %.2f keV' % _energy)
        plt.plot(phase, pol_deg_data[i,:], 'o', color=last_line_color())
    setup_gca(xmin=0., xmax=1., legend=True)


def display_pol_ang(emin=1.1, emax=12., phase_indices=[10, 40, 60, 80],
                    energy_indices=[4, 7, 10, 13]):
    """Display the polarization angle.
    """
    # Polarization angle 2-d plot.
    plt.figure('%s polarization angle' % __model__)
    pol_ang_spline.plot()

    # Slices of the polarization angle at different pulse-phase values.
    plt.figure('%s polarization angle phase slices' % __model__)
    for i in phase_indices:
        _phase = phase[i]
        slice_ = pol_ang_spline.hslice(_phase, k=3)
        slice_.plot(label='Pulse phase = %.2f' % _phase)
        plt.plot(enerlist, pol_ang_data[:,i], 'o', color=last_line_color())
    setup_gca(xmin=emin, xmax=emax, legend=True)

    # Slices of the polarization angle at different energies.
    plt.figure('%s polarization angle energy slices' % __model__)
    for i in energy_indices:
        _energy = enerlist[i]
        slice_ = pol_ang_spline.vslice(_energy, k=3)
        slice_.plot(label='Energy = %.2f keV' % _energy)
        plt.plot(phase, pol_ang_data[i,:], 'o', color=last_line_color())
    setup_gca(xmin=0., xmax=1., legend=True)


def display(emin=1.1, emax=12., phase_indices=[10, 40, 60, 80],
            energy_indices=[4, 7, 10, 13]):
    """Display the source model.
    """
    display_spectrum(emin, emax, phase_indices)
    display_pol_deg(emin, emax, phase_indices, energy_indices)
    display_pol_ang(emin, emax, phase_indices, energy_indices)

if __name__ == '__main__':
    from ixpeobssim.config import bootstrap_display
    bootstrap_display()
