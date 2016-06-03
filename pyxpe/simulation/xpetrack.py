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

# This class take care of physical process of the track in the GPD
# a track is both the photoelectron and the auger electron
# this class combine the TPhoton and TTrack class of xpesim

from pyxpe.simulation.xperandom import xperandom
from pyxpe.simulation.gas       import gasmix
from pyxpe.logging_             import logger
import numpy as np

class xpetrack:
    """Experiment object
    """

    def __init__(self, gas, rnd ):
        """ init basic objects
        """
        self.gas = gas # gas mixture
        self.rnd = rnd # random number generator

    def set_polarization(self, angle, level):
        """ Set polarization angle and degree
        This is a beam property, not a single photon one!
        """
        self.pol_angle = np.deg2rad(angle) # in rad
        self.pol_level = level               # 0 to 1
    
    def set_photon(self, energy, x, y, z):
        """ Set conversion point of the xray 
        """
        self.ph_energy = energy
        self.ph_x = x # conversion point
        self.ph_y = y
        self.ph_z = z
        logger.debug("Photon with E=%.3f (x,y,z)=(%f,%f,%f)" %\
                     (self.ph_energy,self.ph_x,self.ph_y,self.ph_z ))

    def set_detector_bounds(self, xmin, xmax, ymin, ymax, zmin, zmax):
        """
        """
        pass

    def get_photoelectron_theta(self):
        """
        """
        #ELECTRON_MASS = 511. 
        beta = np.sqrt(1.0 - np.power((self.ph_energy/511. + 1),-2.))
        pass
        
    def get_photoelectron_phi(self):
        """
        """
        pass
    
    def extract_phelectron(self):
        """
        """
        # first select converting element
        self.conv_element = self.gas.GetConvertingElement(self.ph_energy,
                                                          self.rnd.random())
        logger.debug("Photon absopt in %s (k-edge %f) "%
                     (self.conv_element.ChemicalSymbol,
                      self.conv_element.kEdge))
        self.res_energy  =  self.ph_energy - self.conv_element.kEdge
        # then extract photoelectron direction

        
        #PhotoelectronTheta = Photon->GetPhotoelectronTheta();
        #PhotoelectronPhi   = Photon->GetPhotoelectronPhi();
        #CX = sin(PhotoelectronTheta)*cos(PhotoelectronPhi);
        #CY = sin(PhotoelectronTheta)*sin(PhotoelectronPhi);
        #CZ = cos(PhotoelectronTheta);

  
if __name__ == '__main__':
    g = gasmix(12, 1.0)
    r = xperandom()
    r.set_seed(666) # diabolic seed
    # test one track
    t = xpetrack(g,r)
    t.set_polarization(30., 1.)
    t.set_photon(5.9, 0, 0, 0.711544053318)
    t.extract_phelectron()

    #
    #ll = []
    #for i in xrange(1):
    #    ll.append(t.extract_phelectron())
    #    with return self.conv_element.AtomicNuber
    #print len(ll), ll.count(8), ll.count(6), ll.count(2), ll.count(1)
    #2.0 10000 6059 3927 10 4
    #3.0 10000 6174 3819 7 0
    #5.9 10000 6304 3689 5 2
    #8.0 10000 6382 3612 5 1
    #10.0 10000 6417 3579 4 0
