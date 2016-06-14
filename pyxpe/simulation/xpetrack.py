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
    """Track object: 
    photoelectron (path and ionization) & Auger electron(path and ionization)

    Coordinates are in mm, angles in rad, energies in keV

    INput quantities:
    set_polarization(angle, level) : polarization information of the photon
    set_photon(x, y, z, energy)    : photon conversion point and energy

    OUTput quantities:
    phe_theta : Theta angle of emission of photoelectron
    phe_phi   : Phi angle of emission of photoelectron;
                Contains the MC Polarization information.
    [phe|aug]_scattering_v : vector of xpepoint objects with information on the
                             photoelectron|auger electron path. 
                             The first item contains the actual ele. emission
    get_ion_pairs() : return the 3 list of (x,y,z) position of ionization pairs;
    get_ion_stats() : return the number of items ions_pairs lists due to
                      photoelectron and auger electrons, 
                      to separate the 2 contributions
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
    
    def set_photon(self, x, y, z, energy):
        """ Set conversion point of the xray 
        """
        self.ph_energy = energy
        self.ph_x = x # conversion point
        self.ph_y = y
        self.ph_z = z
        logger.info("Photon with E=%.3f (x,y,z)=(%f,%f,%f)" %\
                     (self.ph_energy,self.ph_x,self.ph_y,self.ph_z ))

    def set_detector_bounds(self, xmin, xmax, ymin, ymax, zmin, zmax):
        """
        """
        pass

    def reset(self):
        """ RESET all relevant quantities
        """
        # kill photon
        self.ph_energy = 0
        self.ph_x = 0
        self.ph_y = 0
        self.ph_z = 0
        self.phe_theta = None
        self.phe_phi   = None
        self.__cx = 0
        self.__cy = 0
        self.__cz = 0
        # kill track segment
        self.phe_scattering_v = []
        self.aug_scattering_v = []
        # kill ionization
        self.__total_ion_pair  = 0
        self.__total_ion_phe   = 0
        self.__total_ion_auger = 0
        self.__ion_pair_x = np.array([])
        self.__ion_pair_y = np.array([])
        self.__ion_pair_z = np.array([])

    
    def extract_phelectron(self):
        """ First step in building a track: extract photoelectron
        """
        # first select converting element
        self.conv_element = self.gas.GetConvertingElement(self.ph_energy,
                                                          self.rnd.random())
        logger.info("Photon absopt in %s (k-edge %f) "%
                     (self.conv_element.ChemicalSymbol,
                      self.conv_element.kEdge))
        self.res_energy  =  self.ph_energy - self.conv_element.kEdge
        
        # then extract photoelectron direction #ELECTRON_MASS = 511. 
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
        # propagate phe until x-section permits or
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
            nPairs =  self.__get_npairs(self.phe_scattering_v[-1].E)
            self.__total_ion_pair += nPairs
            self.__ion_pair_x = np.append(self.__ion_pair_x, \
                                          np.array(nPairs*[self.phe_scattering_v[-1].x]))
            self.__ion_pair_y = np.append(self.__ion_pair_y, \
                                          np.array(nPairs*[self.phe_scattering_v[-1].y]))
            self.__ion_pair_z = np.append(self.__ion_pair_z, \
                                          np.array(nPairs*[self.phe_scattering_v[-1].z]))
        #
        # Now take care of Auger electron (extract and propagate)
        #
        self.__total_ion_phe   = self.__total_ion_pair # count phe ionization
        self.__total_ion_auger = 0
        # Choose between Auger and fluorescence:
        if self.rnd.random()<= self.conv_element.FluorescenceYield:
            logger.debug("NO Auger electron - No fluorescence photon simulated")
            self.aug_scattering_v = None
            return None
        # Extract Auger electron 
        logger.debug("Auger electron with energy %f" % self.conv_element.kEdge)
        auger_theta = np.arccos(self.rnd.uniform(-1,1))
        auger_phi   = self.rnd.uniform(0, 2*np.pi)
        auger_cx = np.sin(auger_theta)*np.cos(auger_phi);
        auger_cy = np.sin(auger_theta)*np.sin(auger_phi);
        auger_cz = np.cos(auger_theta);
        self.aug_scattering_v = [xpepoint(self.ph_x, self.ph_y, self.ph_z,
                                          self.conv_element.kEdge,
                                          auger_cx, auger_cy, auger_cz)]
        # propagate phe until x-section permits or
        # track outside the detector (TBD)
        while self.aug_scattering_v[-1].E > MIN_ANALYTIC_CROSS_SECTIONS_ENERGY:
            self.aug_scattering_v.append(self.eval_next_point(
                self.aug_scattering_v[-1]))
        # Last electrons are created in the coordinates of the last collision.
        if(self.aug_scattering_v[-1].E <= MIN_ANALYTIC_CROSS_SECTIONS_ENERGY):
           nPairs =  self.__get_npairs(self.aug_scattering_v[-1].E)
           self.__total_ion_pair += nPairs
           self.__ion_pair_x = np.append(self.__ion_pair_x, \
                                         np.array(nPairs*[self.aug_scattering_v[-1].x]))
           self.__ion_pair_y = np.append(self.__ion_pair_y, \
                                         np.array(nPairs*[self.aug_scattering_v[-1].y]))
           self.__ion_pair_z = np.append(self.__ion_pair_z, \
                                         np.array(nPairs*[self.aug_scattering_v[-1].z]))
        # Update number of Auger electron
        # (redundant since is len(ion_pairs)-ion_phe)
        self.__total_ion_auger = self.__total_ion_pair - self.__total_ion_phe
        return None
    
    def eval_next_point(self, phe_point):
        """\brief Evaluate coordinates of the next step in photoelectron path. 
        The formula are taken from Joy's book - 
        note that there is an error in the book: V1=AM*sin(Phi)->V1=AN*sin(Phi)
        Eval also ionization here and update global variables.
        Valid for photoelectron and auger track
        """
        logger.debug("Propagate from %s" %phe_point)
        # get mean free path and extract a random number for path
        # GetElasticMeanFreePath in cm -> mm
        lambd = 10.*self.gas.GetElasticMeanFreePath(phe_point.E, "MOTT")
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
        # GetStoppingPower in keV/cm. -> keV/mm
        energy_loss = min(0.1*path*self.gas.GetStoppingPower(phe_point.E),
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
        nPairs =  self.__get_npairs(energy_loss)
        self.__total_ion_pair += nPairs
        Position = self.rnd.uniform(0,1,nPairs);
        self.__ion_pair_x = np.append(self.__ion_pair_x, \
                                      phe_point.x + path*CA*Position)
        self.__ion_pair_y = np.append(self.__ion_pair_y, \
                                      phe_point.y + path*CB*Position)
        self.__ion_pair_z = np.append(self.__ion_pair_z, \
                                      phe_point.z + path*CC*Position)
        return next_point

    def __get_npairs(self, energy):
        """ \brief Returns the number of e-ion pairs generated between 
        two collisions.
        This version uses Compound Poisson distribution
        """
        MeanNumberSecondary =  1000.*energy/self.gas.WIonization #1000*keV/eV
        NumberPrimary = self.rnd.poisson(MeanNumberSecondary/3.)
        NumberSecondaries = sum(self.rnd.poisson(3., NumberPrimary))
        return NumberSecondaries
    
    def get_ion_pairs(self):
        """ Return the lists of ion pairs position (x,y,z).
        They contains ALL the ions, from photoelectron (at the beginning) 
        and Auger electron (the ending part).
        To divide the list in the 2 contribution, see get_ion_stats()
        """
        return (self.__ion_pair_x, self.__ion_pair_y, self.__ion_pair_z)

    def get_ion_stats(self):
        """ Return the number of photoelectron pairs in the ions_pairs list
        and the number of Auger electron pair in the same list.
        Useful to separate the 2 contributions.
        The sum must be the len of the list.
        """
        return (self.__total_ion_phe, self.__total_ion_auger)

    def get_photoelectron_phi(self):
        """ Get direct info on the photoelectron emission angle
        """
        return self.phe_phi

        
            
def test_theta_phi():
    g = gasmix(12, 1.0)
    r = xperandom()
    r.set_seed(666) # diabolic seed
    t = xpetrack(g,r)
    t.set_polarization(30., 0.5)
    t.set_photon(0, 0, 0.7, 5.9) # z,y,z, E
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
    n = []
    x = []
    y = []
    z = []
    for pl in pts_list:
        n.append(len(pl))
        x.append(np.zeros(n[-1]))
        y.append(np.zeros(n[-1]))
        z.append(np.zeros(n[-1]))
        for i in xrange(n[-1]):
            x[-1][i] = pl[i].x
            y[-1][i] = pl[i].y
            z[-1][i] = pl[i].z
            
    fig = plt.figure()
    ax = fig.add_subplot(211, projection='3d')
    ax.set_xlabel('x [mm]')
    ax.set_ylabel('y [mm]')
    ax.set_zlabel('z [mm]')
    for i in xrange(len(n)):
        plt.plot(x[i], y[i], z[i])
        
    if ion_list !=None:
        plt.plot(ion_list[0], ion_list[1], ion_list[2], 'ro')
    ax1= fig.add_subplot(212)
    ax1.set_xlabel('x [mm]')
    ax1.set_ylabel('x [mm]')
    plt.grid(color='gray')
    for i in xrange(len(n)):
        plt.plot(x[i], y[i])
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
    t.set_photon(0, 0, 0.7, 5.9) # x,y,z, E

    for j in xrange(1):
        print ">>>>>>>>>>>>>>>>>", j
        t.extract_phelectron()
        t.propagate_track()
        print 
        plot_track([t.phe_scattering_v, t.aug_scattering_v], t.get_ion_pairs())

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
