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
import os
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
    parser.add_argument('-o', '--output_file_path', type=str,
                        help='Output file path')
    parser.add_argument('-k', dest='clobber', action='store_true',
                        help='Override existing file')
    args = parser.parse_args()
    startRunNumber = args.start_run_number
    dataFolder = args.folder
    numFrequency = 4
    numShift = 23
    machineId = 2
    if (args.clobber or not os.path.exists(args.output_file_path)):
        f = open(args.output_file_path, 'w')
    else:
        f = open(args.output_file_path, 'a')
    for freqId in range (0, numFrequency):
        for shiftId in range (0, numShift):           
            _run = '{0:03d}_'.format(machineId)
            _run += '{0:07d}'.format(startRunNumber + shiftId 
                                     + numShift * freqId)
            _folder = dataFolder + "/" + _run
            _x0, _y0, _freq, _shift, _charge, _xbarycenter =\
                                                         analyze_run(_folder)
            _displacement = _xbarycenter - _x0
            line = '{} {} {} {} {} {}\n'.format(_x0, _y0, _freq, _shift,
                    _charge, _displacement)
            f.write(line)
    f.close()
