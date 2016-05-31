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

from pyxpe.binio import xpeBinaryFileWindowed



def analyze_injection_run(file_path, num_events):
    """
    """
    assert(file_path.endswith('.mdat'))
    mean = lil_matrix((300, 352))
    rms = lil_matrix((300, 352))
    event_id = 0
    for event in xpeBinaryFileWindowed(file_path):
        print event
        for (col, row) in numpy.ndindex(event.num_columns(), event.num_rows()):
            _col = col + event.xmin
            _row = row + event.ymin
            _adc = event.adc_values[col, row]
            mean[_col, _row] += _adc
            #rms[_col, _row] += _adc**2
        event_id += 1
    mean /= float(event_id)
    print mean
    



if __name__ == '__main__':
    import argparse
    formatter = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=formatter)
    parser.add_argument('binfile', type=str,
                        help='the input binary file')
    parser.add_argument('-n', '--num_events', type=int, default=10,
                        help = 'number of events to be processed')
    args = parser.parse_args()
    analyze_injection_run(args.binfile, args.num_events)
