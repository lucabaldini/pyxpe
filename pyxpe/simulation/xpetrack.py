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
from pyxpe.simulation.gas       import MIN_ANALYTIC_CROSS_SECTIONS_ENERGY
from pyxpe.logging_             import logger
import numpy as np


class xpepoint:
    """ simple class for a single point in the photoelectron track.
    It is identified by a position (xyz), an energy E, 
    a direction with the 3 cosine dir (cx, cy, cz)
    """
    def __init__(self, x,y,z, E, cx, cy, cz):
        """ just a collection of number
        """
        self.x = x
        self.y = y
        self.z = z
        self.E = E
        self.cx = cx
        self.cy = cy
        self.cz = cz

    def __str__(self):
        return "xpepoint (x,y,z,E):(%.3f,%.3f,%.3f,%.4f)" %\
            (self.x, self.y, self.z, self.E)\
            + " dir(x,y,z): (%.2f,%.2f,%.2f)" %\
            (self.cx, self.cy, self.cz)

    def __sub__(self, other):
        return xpepoint(self.x-other.x, self.y-other.y, self.z-other.z,
                        self.E-other.E,
                        self.cx-other.cx, self.cy-other.cy, self.cz-other.cz)
        

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
        self.pol_level = level             # 0 to 1
    
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

    
    def extract_phelectron(self):
        """ First step in building a track: extract photoelectron
        """
        # first select converting element
        self.conv_element = self.gas.GetConvertingElement(self.ph_energy,
                                                          self.rnd.random())
        logger.debug("Photon absopt in %s (k-edge %f) "%
                     (self.conv_element.ChemicalSymbol,
                      self.conv_element.kEdge))
        self.res_energy  =  self.ph_energy - self.conv_element.kEdge
        # then extract photoelectron direction
        #ELECTRON_MASS = 511. 
        beta = np.sqrt(1.0 - np.power((self.ph_energy/511. + 1),-2.))
        self.phe_theta = self.rnd.photoelectron_theta(beta)
        self.phe_phi   = self.rnd.photoelectron_phi(self.pol_angle,
                                                    self.pol_level)
        self.__cx = np.sin(self.phe_theta)*np.cos(self.phe_phi)
        self.__cy = np.sin(self.phe_theta)*np.sin(self.phe_phi)
        self.__cz = np.cos(self.phe_theta)

    def propagate_track(self):
        """ Propagate photoelectron track in the gas.
        This is the most time consuming part, need to find a way to
        avoid this loop....
        """
        self.phe_scattering_v = [xpepoint(self.ph_x, self.ph_y, self.ph_z,
                                          self.res_energy,
                                          self.__cx, self.__cy, self.__cz)]
        # propagate ultil x-section permits or
        # track outside the detector (TBD)
        self.__total_ion_pair = 0
        self.__ion_pair_x = np.array([])
        self.__ion_pair_y = np.array([])
        self.__ion_pair_z = np.array([])
        while self.phe_scattering_v[-1].E > MIN_ANALYTIC_CROSS_SECTIONS_ENERGY:
            self.phe_scattering_v.append(self.eval_next_point(
                self.phe_scattering_v[-1]))
        # Last electrons are created in the coordinates of the last collision.
        if(self.phe_scattering_v[-1].E <= MIN_ANALYTIC_CROSS_SECTIONS_ENERGY
           and len(self.phe_scattering_v)>0):
            nPairs =  self.get_npairs(self.phe_scattering_v[-1].E)
            self.__ion_pair_x = np.append(self.__ion_pair_x, \
                                          np.array(nPairs*[self.phe_scattering_v[-1].x]))
            self.__ion_pair_y = np.append(self.__ion_pair_y, \
                                          np.array(nPairs*[self.phe_scattering_v[-1].y]))
            self.__ion_pair_z = np.append(self.__ion_pair_z, \
                                          np.array(nPairs*[self.phe_scattering_v[-1].z]))

        
    def eval_next_point(self, phe_point):
        """\brief Evaluate coordinates of the next step in photoelectron path. 
        The formula are taken from Joy's book - 
        note that there is an error in the book: V1=AM*sin(Phi)->V1=AN*sin(Phi)
        Eval also ionization here and update global variables.
        """
        logger.debug("Propagate from %s" %phe_point)
        # get mean free path and extract a random number for path
        lambd = self.gas.GetElasticMeanFreePath(phe_point.E, "MOTT") # cm.
        path  = self.rnd.exp(lambd)
        # eval next scattering element and direction 
        # should we use energy at the end of the track?
        re0 = self.rnd.random()
        elem = self.gas.GetScatteringElement(phe_point.E, "MOTT", re0)
        ra0 = self.rnd.random()
        ra1 = self.rnd.random()
        phi = elem.GetScatteringAngle(phe_point.E, "MOTT", ra0, ra1);
        psi = self.rnd.uniform(0, 2*np.pi);
        if phe_point.cz==0:
            V1   = 0;
            V2   = np.sin(phi);
        else:
            AN   = -phe_point.cx/phe_point.cz;
            AM   = 1.0/np.sqrt(1 + AN*AN);
            V1   = AM*np.sin(phi);
            V2   = AN*AM*np.sin(phi);
        #
        V3   = np.cos(psi);
        V4   = np.sin(psi);
        CA   = phe_point.cx*np.cos(phi) + V1*V3 + phe_point.cy*V2*V4;
        CB   = phe_point.cy*np.cos(phi) + V4*(phe_point.cz*V1 -phe_point.cx*V2)
        CC   = phe_point.cz*np.cos(phi) + V2*V3 - phe_point.cy*V1*V4;
        # eval energy loss - must not exceed residual energy
        energy_loss = min(path*self.gas.GetStoppingPower(phe_point.E),
                          phe_point.E);
        # build next point
        next_point = xpepoint(phe_point.x + path*CA,
                              phe_point.y + path*CB,
                              phe_point.z + path*CC,
                              phe_point.E - energy_loss,
                              CA, CB, CC )
        logger.debug("Propagate to %s" %next_point)

        #
        # Eval ionization, stored in 'ion_pair' variables
        #
        nPairs =  self.get_npairs(energy_loss)
        self.__total_ion_pair += nPairs
        Position = self.rnd.uniform(0,1,nPairs);
        self.__ion_pair_x = np.append(self.__ion_pair_x, \
                                      phe_point.x + path*CA*Position)
        self.__ion_pair_y = np.append(self.__ion_pair_y, \
                                      phe_point.y + path*CB*Position)
        self.__ion_pair_z = np.append(self.__ion_pair_z, \
                                      phe_point.z + path*CC*Position)
         
        return next_point

    def get_ion_pairs(self):
        """ 
        """
        return (self.__ion_pair_x, self.__ion_pair_y, self.__ion_pair_z)
        

    def get_npairs(self, energy):
        """ \brief Returns the number of e-ion pairs generated between 
        two collisions.
        This version uses Compound Poisson distribution
        """
        MeanNumberSecondary =  1000.*energy/self.gas.WIonization
        NumberPrimary = self.rnd.poisson(MeanNumberSecondary/3.)
        NumberSecondaries = sum(self.rnd.poisson(3., NumberPrimary))
        return NumberSecondaries
        
            
def test_theta_phi():
    g = gasmix(12, 1.0)
    r = xperandom()
    r.set_seed(666) # diabolic seed
    t = xpetrack(g,r)
    t.set_polarization(30., 0.5)
    t.set_photon(5.9, 0, 0, 0.7) # E, z,y,z
    import matplotlib.pyplot as plt
    tl = []
    pl = []
    for i in xrange(1000):
        t.extract_phelectron()
        tl.append(np.rad2deg(t.phe_theta))
        pl.append(np.rad2deg(t.phe_phi))
    plt.figure(1)
    plt.subplot(211)
    th = plt.hist(np.array(tl), bins=100)
    plt.subplot(212)
    ph = plt.hist(np.array(pl), bins=100)
    plt.show()  


def plot_track(pts_list, ion_list = None):
    import matplotlib as mpl
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    n = len(pts_list)
    x = np.zeros(n)
    y = np.zeros(n)
    z = np.zeros(n)
    for i in xrange(n):
        x[i] = pts_list[i].x
        y[i] = pts_list[i].y
        z[i] = pts_list[i].z

        
    fig = plt.figure()
    ax = fig.add_subplot(211, projection='3d')
    ax.set_xlabel('x [cm]')
    ax.set_ylabel('y [cm]')
    ax.set_zlabel('z [cm]')
    plt.plot(x, y, z)
    if ion_list !=None:
        plt.plot(ion_list[0], ion_list[1], ion_list[2], 'ro')
    fig.add_subplot(212)
    plt.grid(color='gray')
    plt.plot(x, y)
    if ion_list !=None:
        plt.plot(ion_list[0], ion_list[1], 'ro')
    plt.show()
    
    

    
if __name__ == '__main__':
    g = gasmix(12, 1.0)
    r = xperandom()
    r.set_seed(666) # diabolic seed
    # test one track
    logger.setLevel(20) # INFO
    t = xpetrack(g,r)
    t.set_polarization(30., 0.5)
    t.set_photon(5.9, 0, 0, 0.7) # E, x,y,z

    for j in xrange(10):
        print ">>>>>>>>>>>>>>>>>", j
        t.extract_phelectron()
        t.propagate_track()
        plot_track(t.phe_scattering_v, t.get_ion_pairs())
    
    #test_theta_phi()
    
    
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
