#!/usr/bin/env python
#
# Copyright (C) 2016, the ximpol team.
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


import os
import pyregion
import numpy

from ximpol import XIMPOL_CONFIG, XIMPOL_DATA, XIMPOL_EXAMPLES
from ximpol import xpColor
from ximpol.utils.logging_ import logger
from ximpol.core.pipeline import xPipeline
from ximpol.evt.binning import xBinnedMap, xBinnedModulationCube
from ximpol.srcmodel.img import xFITSImage
from ximpol.utils.matplotlib_ import pyplot as plt

from ximpol.config.herx1_pulsar import pol_degree_spline, pol_angle_spline, pl_normalization_spline

from ximpol.utils.os_ import rm
from ximpol.utils.system_ import cmd


"""Script-wide simulation and analysis settings.
"""
model_type = '10on'

base_name = 'herx1_%s_onebin'%model_type

CFG_FILE = os.path.join(XIMPOL_CONFIG, 'herx1_pulsar.py')
OUT_FILE_PATH_BASE = os.path.join(XIMPOL_DATA, base_name)

MCUBE_FILE_PATH = '%s_mcube.fits'%OUT_FILE_PATH_BASE
ANALYSIS_FILE_PATH = '%s_analysis.txt' % OUT_FILE_PATH_BASE
EVT_FILE_PATH = '%s.fits' % OUT_FILE_PATH_BASE
SIM_DURATION = 1000000.

#2.0 - 3.0 keV, 3.0 - 4.5 keV,  4.5 keV - 6.0 keV, 6.0 - 8.0 

E_BINNING = [2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.8, 4.0, 4.5, 6.0, 8.0]
"""Main pipeline object.
"""
PIPELINE = xPipeline(clobber=False)

#Added this method in run() so that we can simulate several runs with different seeds and merge the output files. This is needed for bright sources which have large amount of events in the output file.

def run(repeat=2):
    #First simulate the events
    file_list = []
    for i in range(repeat):
        output_file_path = EVT_FILE_PATH.replace('.fits', '_%d.fits' % i)
        file_list.append(output_file_path)
        PIPELINE.xpobssim(configfile=CFG_FILE, duration=SIM_DURATION,
                          outfile=output_file_path, seed=i)
    file_list = str(file_list).strip('[]').replace('\'', '').replace(' ', '')
    if PIPELINE.clobber:
        rm(EVT_FILE_PATH)
    cmd('ftmerge %s %s' % (file_list, EVT_FILE_PATH))

    PIPELINE.xpbin(EVT_FILE_PATH, algorithm='MCUBE', ebinalg='LIST',
                       ebinning=E_BINNING)
    
def analyze():

    """Analyze the data.Testing this method, but I must be missing something, it does not work yet.
    """
    logger.info('Opening output file %s...' % ANALYSIS_FILE_PATH)
    analysis_file = open(ANALYSIS_FILE_PATH, 'w')
    _mcube = xBinnedModulationCube(MCUBE_FILE_PATH)
    _mcube.fit()
    for j, fit in enumerate(_mcube.fit_results):
    #_fit_results = _mcube.fit_results[0]
        _pol_deg = fit.polarization_degree
        _pol_deg_err = fit.polarization_degree_error
        _pol_angle = fit.phase
        _pol_angle_err = fit.phase_error
        _energy_mean = _mcube.emean[j]
        _emin = _mcube.emin[j]
        _emax = _mcube.emax[j]
        _data = (_energy_mean,_emin, _emax,_pol_deg, _pol_deg_err, _pol_angle, _pol_angle_err)
        print _data
        _fmt = ('%.4e   ' * len(_data)).strip()
        _fmt = '%s\n' % _fmt
        _line = _fmt % _data
        analysis_file.write(_line)
    analysis_file.close()



def view():
    #_energy_mean,_emin, _emax, _pol_deg, _pol_deg_err, _pol_angle, \
    #    _pol_angle_err = \
     #   #                numpy.loadtxt(ANALYSIS_FILE_PATH, unpack=True)
    
    _mcube = xBinnedModulationCube(MCUBE_FILE_PATH)
    _mcube.fit()
    _fit_results = _mcube.fit_results[0]
    plt.figure('Polarization degree')
    _mcube.plot_polarization_degree(show=False, color='blue')
    pol_degree_spline.plot(color='lightgray',label='Model %s'%model_type, show=False)
    plt.figtext(0.2, 0.85,'XIPE %s ks'%(SIM_DURATION/1000.),size=18)
    #plt.errorbar(_energy_mean, _pol_deg, yerr=_pol_deg_err, color='blue',marker='o')
    
    plt.legend()

    plt.figure('Polarization angle')
    _mcube.plot_polarization_angle(show=False, color='blue', degree=False)
    pol_angle_spline.plot(color='lightgray',label='Model %s'%model_type, show=False)
    plt.figtext(0.2, 0.85,'XIPE %s ks'%(SIM_DURATION/1000.),size=18)
    #plt.errorbar(_energy_mean,_pol_angle, yerr= _pol_angle_err,color='blue',marker='o')

    plt.legend()
    plt.figure('MDP %s'%base_name)
    mdp = _mcube.mdp99
    emean = _mcube.emean
    emin =  _mcube.emin
    emax =  _mcube.emax
    width = (emax-emin)/2.
    plt.errorbar(emean,mdp,xerr=width, label='MDP99',marker='o',linestyle='--')
    plt.figtext(0.2, 0.85,'XIPE %s ks'%(SIM_DURATION/1000.),size=18)
    plt.xlim([1,10])
    plt.ylabel('MPD 99\%')
    plt.xlabel('Energy (keV)')
    #plt.legend()
    plt.show()

if __name__ == '__main__':
    run()
    analyze()
    view()