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
import numpy as np

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
        # one the source is defined, we can init the track object
        self.track = xpetrack(self.gas, self.rnd)
        
        # maybe we can move this in photon extraction to get
        # different polarization for different photons..
        self.track.set_polarization(self.pol_angle, self.pol_degree)

    def extract_photon(self):
        """ Extract a random photon on the surface of the detector.
        should be from Source object, now just 3 fixed numbers
        """
        energy = self.energy
        x = self.cntr_x
        y = self.cntr_y
        return (energy, x, y)

    def process_event(self, i):
        """ Do all the calculation in this function
        """
        # get event on the surface of the detector
        energy, x, y = self.extract_photon() 
        # extract the conversion along the depth of the gas chamber
        # and repeat this step untill the conversion is in the active volume
        #
        # Lambda must be in mm, but X-sec*density=[cm2/g][g/cm3]= [1/cm]
        Lambda = 0.1*self.gas.GetPhotoelectricCrossSection(energy)*\
                 self.gas.Density
        z = -10
        while z<self.gem_z: # Z=0 is GEM position
            z = self.tickness - self.rnd.exp(1./Lambda) #rnd->Exp(1./Lambda)

        # now that the event converted, eval the track
        self.track.reset()
        self.track.set_photon(x, y, z, energy)
        self.track.extract_phelectron()
        self.track.propagate_track()

    def plot_event(self):
        """ for DEBUG, plot the event
        """
        import matplotlib as mpl
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pyplot as plt
        n = len(self.track.phe_scattering_v)
        logger.info("Found %d elements in photo-electron track" % n)
        x = np.zeros(n)
        y = np.zeros(n)
        z = np.zeros(n)
        an = len(self.track.aug_scattering_v)
        logger.info("Found %d elements in auger electron track" % an)
        ax = np.zeros(an)
        ay = np.zeros(an)
        az = np.zeros(an)
        for i in xrange(n):
            x[i] = self.track.phe_scattering_v[i].x
            y[i] = self.track.phe_scattering_v[i].y
            z[i] = self.track.phe_scattering_v[i].z
        for i in xrange(an):
            ax[i] = self.track.aug_scattering_v[i].x
            ay[i] = self.track.aug_scattering_v[i].y
            az[i] = self.track.aug_scattering_v[i].z
        fig = plt.figure()
        axx = fig.add_subplot(221, projection='3d')
        axx.set_xlabel('x [mm]')
        axx.set_ylabel('y [mm]')
        axx.set_zlabel('z [mm]')
        plt.plot(x, y, z)
        plt.plot(ax, ay, az)
        ion_list = self.track.get_ion_pairs()
        logger.info("Found %d ion pairs" % len(ion_list[0])+\
                    " - %d from photoelectron" % self.track.get_ion_stats()[0])
        plt.plot(ion_list[0], ion_list[1], ion_list[2], 'ro')
        axx1= fig.add_subplot(222)
        axx1.set_xlabel('x [mm]')
        axx1.set_ylabel('x [mm]')
        plt.grid(color='gray')
        plt.plot(x, y)
        plt.plot(ax, ay)
        plt.plot(ion_list[0], ion_list[1], 'ro')
        plt.show()
        
if __name__ == '__main__':
    logger.setLevel(20) # INFO
    e = experiment(666)
    # simple setup
    e.set_gas(12, 1.0) #Gas ID: 12 (HeDME 20-80); Pressure (atm): 1.0
    e.set_detector(10.0) # Thickness (cm): 1.0
    e.set_source(0,0,0,0) # pol_angle, pol_degree, cntr_x, cntr_y
    # run events
    #print 10./(e.gas.GetPhotoelectricCrossSection(e.energy)*e.gas.Density)
    import time
    t0 = time.time()
    for i in xrange(5):
        logger.info("*********** EVENT %d" % i)
        e.process_event(i)
        e.plot_event()
    t0 = time.time()-t0
    print i, t0, (i+1)/t0
    

    
