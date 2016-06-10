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

from pyxpe.simulation.xperandom import xperandom
from pyxpe.simulation.gas       import gasmix
from pyxpe.simulation.xpetrack  import xpetrack
from pyxpe.logging_             import logger 


class experiment:
    """Experiment object
    """

    def __init__(self, rndseed = None):
        """ The experiment parameter should be passed here
        """
        self.rnd = xperandom()
        self.rnd.set_seed(rndseed) # time-based seed

    def set_gas(self, mixid, pressure):
        """ Set gas properties, mainly via mixture id
        """
        self.gas = gasmix(mixid, pressure)
        

    def set_detector(self, tickness):
        """ For now just keep relevant parameters.
        In the (near) future implement a detector + GEM class with all 
        the right parameters
        """
        self.tickness = tickness
        self.gem_z    = 0 #let's put the GEM at Z=0


    def set_source(self, pol_angle, pol_degree, cntr_x, cntr_y):
        """ For now just the position in an arbitrary reference frame.
        In the (near) future, a source class is needed.
        """
        self.energy     = 5.9 # keV
        self.pol_angle  = pol_angle
        self.pol_degree = pol_degree
        self.cntr_x     = cntr_x
        self.cntr_y     = cntr_y

    def extract_photon(self):
        """ Extract a random photon on the surface of the detector.
        should be from Source object, now just 3 fixed numbers
        """
        energy = self.energy
        x = self.cntr_x
        y = self.cntr_y
        return (energy, x, y)

    def process_event(self, i):
        """
        """
        # get event on the surface of the detector
        energy, x, y = self.extract_photon() 
        # extract the conversion along the depth of the gas chamber
        # and repeat this step untill the conversion is in the active volume
        #
        # Lambda must be in mm, but X-sec*density=[cm2/g][g/cm3]= [1/cm]
        Lambda = 0.1*self.gas.GetPhotoelectricCrossSection(energy)*self.gas.Density
        z = -10
        while z<self.gem_z: # Z=0 is GEM position
            z = self.tickness - self.rnd.exp(1./Lambda) #rnd->Exp(1./Lambda)
        print 1./Lambda, z
        # now that the event converted, init a track object
        self.track = xpetrack(self.gas, self.rnd)
        self.track.set_photon(energy, x, y, z)

        

if __name__ == '__main__':
    e = experiment()
    # simple setup
    e.set_gas(12, 1.0) #Gas ID: 12 (HeDME 20-80); Pressure (atm): 1.0
    e.set_detector(10.0) # Thickness (cm): 1.0
    e.set_source(0,0,0,0)
    # run events
    #print 10./(e.gas.GetPhotoelectricCrossSection(e.energy)*e.gas.Density)
    for i in xrange(1):
        e.process_event(i)
        
