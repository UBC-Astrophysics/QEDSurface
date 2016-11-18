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
import numpy

from ximpol import XIMPOL_CONFIG, XIMPOL_DATA, XIMPOL_DOC
from ximpol.utils.logging_ import logger
from ximpol.core.pipeline import xPipeline
from ximpol.evt.binning import xBinnedModulationCube, xEventBinningBase
from ximpol.evt.event import xEventFile
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.utils.matplotlib_ import save_current_figure
from ximpol.config.herx1_pulsar10on import pol_degree_spline, \
    pol_angle_spline, pl_normalization_spline

from ximpol.utils.os_ import rm
from ximpol.utils.system_ import cmd


"""Script-wide simulation and analysis settings.
"""
CFG_FILE_PATH = os.path.join(XIMPOL_CONFIG, 'herx1_pulsar10on.py')
OUT_FILE_PATH_BASE = os.path.join(XIMPOL_DATA, 'herx1_pulsar_twobin100ks')
EVT_FILE_PATH = '%s.fits' % OUT_FILE_PATH_BASE
ANALYSIS_FILE_PATH = '%s_analysis.txt' % OUT_FILE_PATH_BASE
SIM_DURATION = 100000.
NUM_PHASE_BINS = 2
EQP_BINNING = False
PHASE_BINNING = None
E_BINNING = [1., 10.]
OUTPUT_FOLDER = os.path.join(XIMPOL_DOC, 'figures', 'showcase')


"""Main pipeline object.
"""
PIPELINE = xPipeline(clobber=False)


def _sel_file_path(i):
    """Return the path to the i-th xpselct output file.
    """
    return '%s_phase%04d.fits' % (OUT_FILE_PATH_BASE, i)

def _mcube_file_path(i):
    """Return the path to the i-th xpbin MCUBE output file.
    """
    return '%s_phase%04d_mcube.fits' % (OUT_FILE_PATH_BASE, i)

def _pha1_file_path(i):
    """Return the path to the i-th xpbin PHA1 output file.
    """
    return '%s_phase%04d_pha1.fits' % (OUT_FILE_PATH_BASE, i)

def _phase_binning():
    """Read the input event file and create an equipopulated binning in the
    pulsar phase.
    """
    if EQP_BINNING:
        evt_file = xEventFile(EVT_FILE_PATH)
        phase = evt_file.event_data['PHASE']
        return xEventBinningBase.equipopulated_binning(NUM_PHASE_BINS, phase,
                                                       0., 1.)
    else:
        return numpy.linspace(0., 1., NUM_PHASE_BINS)


def generate(count=10):
    """Generate the events.
    """
    #First simulate the events
    PIPELINE.xpobssim(configfile=CFG_FILE_PATH, duration=SIM_DURATION,
                      outfile=EVT_FILE_PATH)
'''
    file_list = []
    for i in range(count):
        output_file_path = EVT_FILE_PATH.replace('.fits', '_%d.fits' % i)
        file_list.append(output_file_path)
        PIPELINE.xpobssim(configfile=CFG_FILE_PATH,
                          duration=SIM_DURATION/count,
                          outfile=output_file_path, seed=i)
    file_list = str(file_list).strip('[]').replace('\'', '').replace(' ', '')
    if PIPELINE.clobber:
        rm(EVT_FILE_PATH)
    cmd('ftmerge %s %s' % (file_list, EVT_FILE_PATH))
'''
#2.0 - 3.2 keV, 3.2 - 4.5 keV,  4.5 keV - 6.5 keV, 6.0 - 8.0 

E_BINNING = [2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 6.0, 8.0]

    
def prepare():
    """Prepare the event data for the actual analysis.
    """
    for i, (_min, _max) in enumerate(zip(PHASE_BINNING[:-1],
                                         PHASE_BINNING[1:])):
        PIPELINE.xpselect(EVT_FILE_PATH, phasemin=_min, phasemax=_max,
                          outfile=_sel_file_path(i))
        PIPELINE.xpbin(_sel_file_path(i), algorithm='MCUBE', ebinalg='LIST',
                       ebinning=E_BINNING, outfile=_mcube_file_path(i))
        PIPELINE.xpbin(_sel_file_path(i), algorithm='PHA1',
                       outfile=_pha1_file_path(i))

def analyze():
    """Analyze the data.
    """

    for j in range(len(E_BINNING)):
        if (j==len(E_BINNING)-1):
            logger.info('Opening output file %s...' % ANALYSIS_FILE_PATH)
            analysis_file = open(ANALYSIS_FILE_PATH, 'w')
        else:
            logger.info('Opening output file %s_%d ...' % (ANALYSIS_FILE_PATH,j))
            analysis_file = open('%s_%d' % (ANALYSIS_FILE_PATH,j), 'w')
        for i, (_min, _max) in enumerate(zip(PHASE_BINNING[:-1],
                                             PHASE_BINNING[1:])):
            _mcube = xBinnedModulationCube(_mcube_file_path(i))
            _mcube.fit()
            _fit_results = _mcube.fit_results[j]
            _emean = _mcube.emean[j]
            _emin = _mcube.emin[j]
            _emax = _mcube.emax[j]
            print("%g - %g - %g\n" % (_emin,_emean,_emax))
            _phase = 0.5*(_min + _max)
            _phase_err = 0.5*(_max - _min)
            _pol_deg = _fit_results.polarization_degree
            _pol_deg_err = _fit_results.polarization_degree_error
            _pol_angle = _fit_results.phase
            _pol_angle_err = _fit_results.phase_error
            _spec_fitter = PIPELINE.xpxspec(_pha1_file_path(i), plot=False,model='powerlaw')
            (_index, _index_err), (_norm, _norm_err) = _spec_fitter.fit_parameters()
            # The division by the phase interval is a workaround and we should
            # keep track of that in xpselect.
            _norm /= (_max - _min)
            _norm_err /= (_max - _min)
            _data = (_phase, _phase_err, _pol_deg, _pol_deg_err, _pol_angle,
                     _pol_angle_err, _index, _index_err, _norm, _norm_err)
            _fmt = ('%.4e   ' * len(_data)).strip()
            _fmt = '%s\n' % _fmt
            _line = _fmt % _data
            analysis_file.write(_line)
        analysis_file.close()

def plot(save=False):
    """Plot the stuff in the analysis file.
    """
    for j in range(len(E_BINNING)):
        if (j==len(E_BINNING)-1):
            analysis_file = ANALYSIS_FILE_PATH
            emin=E_BINNING[0]
            emax=E_BINNING[-1]
        else:
            analysis_file = '%s_%d' % (ANALYSIS_FILE_PATH,j)
            emin=E_BINNING[j]
            emax=E_BINNING[j+1]
            
        sim_label = 'XIPE %s ks' % (SIM_DURATION/1000.)
        mod_label = 'Input model'
        lc_label = 'Light curve'
        _phase, _phase_err, _pol_deg, _pol_deg_err, _pol_angle,\
            _pol_angle_err, _index, _index_err, _norm,\
            _norm_err = numpy.loadtxt(analysis_file, unpack=True)
        plt.figure('Polarization degree (%g-%g keV)' % (emin,emax))
        plt.title('Her X-1 (%g-%g keV)' % (emin,emax))
        pl_normalization_spline.plot(scale=20, show=False, color='lightgray',
                                     label=lc_label)
        plt.errorbar(_phase, _pol_deg, xerr=_phase_err, yerr=_pol_deg_err, fmt='o',
                     label=sim_label)
        pol_degree_spline.plot(show=False, label=mod_label)
        plt.axis([0., 1., 0., 1.1])
        plt.legend(bbox_to_anchor=(0.4, 0.5))
        if save:
            save_current_figure('herx1_polarization_degree_%d' % j, OUTPUT_FOLDER, False)
        plt.figure('Polarization angle (%g-%g keV)' % (emin,emax))
        plt.title('Her X-1 (%g-%g keV)' % (emin,emax))
        pl_normalization_spline.plot(scale=60, offset=0, show=False,
                                     color='lightgray', label=lc_label)
        plt.errorbar(_phase, _pol_angle, xerr=_phase_err, yerr=_pol_angle_err,
            fmt='o', label=sim_label)
        pol_angle_spline.plot(show=False, label=mod_label)
        plt.axis([0., 1., 0, 3.15])
        plt.legend(bbox_to_anchor=(0.4, 0.3))
        if save:
            save_current_figure('herx1_polarization_angle_%d' %j, OUTPUT_FOLDER, False)
        plt.figure('Flux (%g-%g keV)' % (emin,emax))
        plt.title('Her X-1 (%g-%g keV)' % (emin,emax))
        plt.errorbar(_phase,
                     0.21*_norm/(2-_index)*(10.**(2.0-_index)-2.**(2.-_index)), xerr=_phase_err,
                     yerr=0.21*_norm_err/(2-_index)*(10.**(2-_index)-2.**(2.-_index)), fmt='o',
                     label=sim_label)
        pl_normalization_spline.plot(scale=1,show=False, label=mod_label)
        plt.axis([0., 1., None, None])
        plt.legend(bbox_to_anchor=(0.75, 0.95))
        if save:
            save_current_figure('herx1_flux_%d' % j, OUTPUT_FOLDER, False)
#        plt.show()

def run(save_plots=False):
    """Run all the tasks.
    """
    if os.path.exists(ANALYSIS_FILE_PATH):
        logger.info('%s exists, delete it if you want to recreate it.' %\
                    ANALYSIS_FILE_PATH)
    else:
        generate()
        global PHASE_BINNING
        PHASE_BINNING = _phase_binning()
        prepare()
        analyze()
    plot(save_plots)


if __name__ == '__main__':
    run(save_plots=True)
