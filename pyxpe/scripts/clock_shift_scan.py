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
from scipy.sparse import lil_matrix

from pyxpe.recon.injection import analyze_run

if __name__ == '__main__':
    import argparse
    formatter = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=formatter)
    parser.add_argument('folder', type=str,
                        help='the data folder (containing all the run folders)')
    parser.add_argument('-s', '--start_run_number', type=int,
                        help='initial run number')
    parser.add_argument('-o', '--output_folder', type=str,
                        help='Output folder')
    parser.add_argument('-l', '--label', type=str,
                        help='A custom label appended to the output file names')                                                
    args = parser.parse_args()
    startRunNumber = args.start_run_number
    dataFolder = args.folder
    numFrequency = 4
    numShift = 23
    machineId = 2
    shifts = numpy.zeros(shape=(numFrequency,numShift))
    displacements = numpy.zeros(shape=(numFrequency,numShift))
    figName  = "DispVsClshift"
    plt.figure(figName)
    plt.xlabel("Clock shift [ns]")
    plt.ylabel("Displacement [pixels]")
    _x0 = 0
    _y0 = 0
    for freqId in range (0, numFrequency):
        for shiftId in range (0, numShift):           
            _run = '{0:03d}_'.format(machineId)
            _run += '{0:07d}'.format(startRunNumber + shiftId 
                                     + numShift * freqId)
            _folder = dataFolder + "/" + _run
            _x0, _y0, _freq, _shift, _xbarycenter = analyze_run(_folder)
            shifts[freqId][shiftId] = _shift
            displacements[freqId][shiftId] = _xbarycenter - _x0
        plt.plot(shifts[freqId, : ], displacements[freqId, : ],
                 label="%.2f MHz"%_freq)
    plt.title("Pixel (%d , %d)" % (_x0, _y0))                
    plt.xlim([220., 850.])
    ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(handles, labels, loc =  9)
    if args.output_folder is None:
        plt.show()
    else:
        outFileNameBase = figName + '_%d_%d'%(_x0, _y0)
        if args.label is not None:
            outFileNameBase = outFileNameBase + "_" + args.label
        outFilePath = args.output_folder + '/' + outFileNameBase
        plt.savefig(outFilePath + '.png')
        plt.ylim([-0.06, 0.06])
        plt.xlim([400., 800.])
        plt.savefig(outFilePath + '_zoomed.png')
    
    
