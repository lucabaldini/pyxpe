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


import os

from pyxpe.utils.logging_ import logger



class xpeDetectorConfiguration(dict):

    """Container encapsulating the detector configuration.
    """

    CLOCK_FREQ_DICT = {96: 1.25, 64: 2.5, 32: 5., 0: 10.}

    def __init__(self, file_path=None, **kwargs):
        """
        """
        if file_path is not None:
            self.__load_from_file(file_path) 
        for key, value in kwargs.items():
            self[key] = value

    def __load_from_file(self, file_path):
        """Load the configuration from file.
        """
        logger.info('Loading configuration from %s...' % file_path)
        key = None
        for line in open(file_path).readlines():
            line = line.rstrip()
            if line.startswith('//'):
                if key is not None:
                    if len(value) == 1:
                        value = value.pop()
                    self[key] = value
                key = line.strip('/ #').lower().replace(' ', '_')
                value = []
            else:
                value.append(int(line))

    def clock_shift_ns(self):
        """Return the clock shift in ns.
        """
        return self['clock_shift']*25

    def clock_freq_mhz(self):
        """Return the clock frequency in MHz.
        """
        return self.CLOCK_FREQ_DICT[self['clock_frequency_code']]
    


def load_run_info(folder_path):
    """Load all the relevant run information from a folder written by the
    DAQ.
    """
    folder_path = os.path.normpath(folder_path)
    folder_name =  os.path.basename(folder_path)
    logger.info('Loading run info from %s...' % folder_path)
    station_id, run_id = [int(item) for item in folder_name.split('_')]
    data_file_path = os.path.join(folder_path, '%s_data.mdat' % folder_name)
    config_file_path = os.path.join(folder_path, '%s_detector.cfg' %\
                                    folder_name)
    config = xpeDetectorConfiguration(config_file_path)
    info = {}
    info['station_id'] = station_id
    info['run_id'] = run_id
    info['data_file_path'] = data_file_path
    info['config_file_path'] = config_file_path
    info['config'] = config
    return info


if __name__ == '__main__':
    folder_path = '/data/work/xpedata/timing/033_0000432/'
    print(load_run_info(folder_path))
