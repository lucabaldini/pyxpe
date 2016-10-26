#!/usr/bin/env python
# Copyright (C) 2007--2016 the X-ray Polarimetry Explorer (XPE) team.
#
# For the license terms see the file LICENSE, distributed along with this
# software.
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

import numpy
import matplotlib.pyplot as plt

from pyxpe.recon.injection import compute_mean_rms
from pyxpe.recon.daq import load_run_info
from pyxpe.recon.binio import xpeBinaryFileWindowed
import os, sys

# if you want a specific, single run, pass it as argument
try:
    runIn = int(sys.argv[1])
except IndexError: 
    runIn = None

# notice that event on the edge can have a smaller roi
SINGLE_PXL_WINDOW = [396, 324] 

def get_ci_run_noise_stat(run_id, verbose = False):
    """ Noisy channel analysis of charge injection runs
    """
    print "Processing run Id %d" % run_id
    run_info = load_run_info('/data1/xpe/xpedata/033_%07d' % run_id)
    thr_dac = run_info['config']['threshold_dacs'][0]
    if verbose:
        print run_info
    nMinWin    = 0
    nLargerWin = 0
    RoiDict = {}
    file_path = run_info['data_file_path']
    if not os.path.isfile(file_path):
        # No event found, upper limits not implemented yet
        return ( thr_dac, 0, 0, -1, {})
    for event in xpeBinaryFileWindowed(file_path):
        # check event has minimum window of 396
        if event.num_pixels() in SINGLE_PXL_WINDOW:
            nMinWin +=1
            # Eval ROI and check if already found, if not add it
            roi = event.roi()
            if RoiDict.has_key(roi):
                RoiDict[roi] +=1
            else: 
                RoiDict[roi] =1
        else:
            nLargerWin +=1
            if verbose:
                print "******************* ROI different than Minimum:"
                print event.roi(),  event.num_pixels()
    timeLastEvt = event.microseconds
    # dac, nGoodRoi, nBadRoi, timeLastEvt, RoiDict
    return ( thr_dac, nMinWin, nLargerWin, timeLastEvt, RoiDict)


def thr_dac2mV(dac):
    """ Convert threshold DAC to mV value
    """
    mV = 0.8*dac +2
    mV = int(round(mV))
    return mV


#
# DO THE ACTUAL WORK:
#


#rList = [585]#578, 579, 580, 581, 582, 583, 584]

# run list for XPOL_PI_NR_4
rList = [570, 569, 568, 567, 571]
label = "XPOL_PI_NR_4"

# run list for XPOL_PI_NR_3
rList = [597, 596, 595, 594, 593, 592, 591, 590] # 592, 591, 590 with 0 counts
label = "XPOL_PI_NR_3"


th   = []
rate = []
nch  = []
if runIn !=None:
   rList = [runIn] 
   label = "%d" %runIn

for run in rList:
    (thr_dac, nRoi, nRoiBad, tLastEvt, RoiDict) = get_ci_run_noise_stat(run)

    print "----------------------------------------------"
    print "Run summary for id %d" % run
    print "Threshold DAC = %d (%d mV)" % (thr_dac, thr_dac2mV(thr_dac))
    print "Number of events %d (roi != than min = %d)" % \
        (nRoi +nRoiBad, nRoiBad)
    ave_rate = (nRoi +nRoiBad)/(tLastEvt/1e6)
    print "Average Rate = %f Hz" % ave_rate
    n_roi = len(RoiDict.keys())
    print "Number of different ROI = %d " %  n_roi
    print "-----------------------------"
    print "Time of last evt %f us" % tLastEvt
    print "ROI printout: "
    n = 0
    for rr in RoiDict.keys():
        n += RoiDict[rr]
        print  rr, RoiDict[rr]
    print "Total event in Dict = %d"%n
    if nRoi>0: # upper limits not implemented yet
        th.append(thr_dac2mV(thr_dac))
        rate.append(ave_rate)
        nch.append(n_roi)
    if False and runIn ==None:
        raw_input('...inspect...')


# do some plot if it's the case
if len(th)>0:
    th   = numpy.array(th)
    rate = numpy.array(rate)
    nch = numpy.array(nch)
    print th
    plt.plot(th, rate, 'o', label='rate [Hz]')
    plt.plot(th, nch, 'o', label='# noisy ch')
    plt.grid()
    plt.xlabel('Threshold [mV]')
    plt.title('Threshold Scan %s' %label)
    plt.legend()
    plt.yscale('log')
    plt.axis([160, 360, 0.1, 1000])

    plt.show()
