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

# Main class for running a simulation job (an `experiment')
# follows roughly the concept of xpesim/TExperiment

import xperandom
from math import sqrt
from gas import gasmix
import logging
logging.basicConfig(format='%(module)s:%(levelname)s:%(message)s',\
                    level=logging.INFO)

class experiment:
    """Experiment object
    """

    def __init__(self, rndseed = None):
        """ The experiment parameter should be passed here
        """
        self.rnd = xperandom.xperandom()
        self.rnd.set_seed(rndseed) # time-based seed

    def set_gas(self, mixid, pressure):
        """ Set gas properties, mainly via mixture id
        """
        self.gas = gasmix(mixid, pressure)
        

    def set_detector(self, tickness):
        """
        """
        self.tickness = tickness


    def set_source(self, pol_angle, pol_degree, cntr_x, cntr_y):
        """
        """
        self.pol_angle  = pol_angle
        self.pol_degree = pol_degree
        self.cntr_x     = cntr_x
        self.cntr_y     = cntr_y



        

if __name__ == '__main__':
     e = experiment()
     # simple setup
     e.set_gas(12, 1.0) #Gas ID: 12 (HeDME 20-80); Pressure (atm): 1.0
     e.set_detector(10.0) # Thickness (cm): 1.0
     e.set_source(0,0,0,0)
