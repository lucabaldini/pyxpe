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


# Wrapper lib for random number generation
# to change engine (python, root, heprep), by changing only this file

# old choice was python pseudo-random generator
# https://docs.python.org/2.7/library/random.html
#import random

# current default is numpy
# http://docs.scipy.org/doc/numpy/reference/routines.random.html
import numpy as np 
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt


class xpeUnivariateGenerator(InterpolatedUnivariateSpline):

    """
    """

    def __init__(self, rv, pdf, w=None, bbox=[None, None], k=1):
        """Constructor.
        """
        self.__rv = rv
        self.__pdf = pdf
        InterpolatedUnivariateSpline.__init__(self, rv, pdf, w, bbox, k)
        self.ppf = self.build_ppf()

    def norm(self):
        """Return the integral over the entire spline domain.
        """
        return self.integral(self.__rv[0], self.__rv[-1])

    def build_ppf(self):
        """Create the percent point function (or inverse of the cdf).
        """
        _y = self.__rv
        _xmin = self.__rv[0]
        _x = np.array([self.integral(_xmin, _xp) for _xp in _y])
        _x, _mask = np.unique(_x, return_index=True)
        _x/= self.norm()
        _y = _y[_mask]
        return InterpolatedUnivariateSpline(_x, _y)
    
    def pdf(self, rv):
        """Return the pdf value(s) at the point(s) rv.
        """
        return self(rv)

    def rvs(self, size=1):
        """Return random variates of arbitrary size.
        """
        return self.ppf(np.random.sample(size))

    def plot(self, show=True):
        """
        """
        plt.plot(self.__rv, self.__pdf)
        if show:
            plt.show()
        


class xperandom:
    """Wrapper random class
    """

    def __init__(self):
        self.engine = np.random

    def get_engine(self):
        """ Return the engine, in order to use it directly
        """
        return self.engine
    
    def set_seed(self, seed=None):
        """ Set the seed of the pseudo-random generator
        """
        self.engine.seed(seed)

    def random(self):
        """ Main engine call
        Return the next random floating point number in the range [0.0, 1.0).
        """
        return self.engine.random()

    def get_state(self):
        """
        """
        return self.engine.get_state()

    def set_state(self, state):
        """
        """
        return self.engine.set_state(state)
    
    def exp(self, l, n = 1):
        """ e**(x/l) as in ROOT:TRandom:Exp
        docs.scipy.org/doc/numpy/reference/
        generated/numpy.random.exponential.html
        """
        return self.engine.exponential(l, n)

    def uniform(self, a, b, n=1):
        """ http://docs.scipy.org/doc/numpy/reference/
        generated/numpy.random.uniform.html
        """
        return self.engine.uniform(a, b, n)

    def poisson(self, m, n=1):
        """ http://docs.scipy.org/doc/numpy/reference/
        generated/numpy.random.poisson.html
        """
        return self.engine.poisson(m, n)

    def multigauss(self, m, s, n = 1):
        """ Multivariate Gauss (e.g. for multiple gauss with different sigma)
        http://docs.scipy.org/doc/numpy/reference/
        generated/numpy.random.multivariate_normal.html
        """
        return  self.engine.multivariate_normal(m, s, n)

    def photoelectron_theta(self, beta, nevt = 1):
        """ Get photoelectron theta
        """
        return self.__photoelectron_theta_1(beta, nevt)
    
    def photoelectron_phi(self, pol_angle, pol_level, nevt = 1):
        """ Get photoelectron phi for beam polarized at level 
        pol_level (from 0 to 1) and angle pol_angle (deg)
        """
        return self.__photoelectron_phi_1(pol_angle, pol_level, nevt)
    
    def __photoelectron_theta_1(self, beta, nevt = 1):
        """ Get photoelectron theta - spline 1d
        """
        rv = np.linspace(0, np.pi, 100)
        pdf = np.sin(rv)**3/((1-beta*np.cos(rv))**4)
        generator = xpeUnivariateGenerator(rv, pdf)
        values = np.pi - generator.rvs(nevt) # ref system with +Z up
        if nevt==1:
            return values[0]
        return values

    def __photoelectron_phi_1(self, pol_angle, pol_level, nevt = 1):
        """ Get photoelectron phi for beam polarized at level 
        pol_level (from 0 to 1) and angle pol_angle (rad) - spline 1d
        """
        rv = np.linspace(0, 2*np.pi, 100)
        pdf = (1-pol_level) + (2*pol_level)*np.cos(rv+ pol_angle)**2
        generator = xpeUnivariateGenerator(rv, pdf)
        values = generator.rvs(nevt)
        if nevt==1:
            return values[0]
        return values
    
def test_random():
    """
    """
    r = xperandom()

    print "Test direct vs wrapped call:"
    r.get_engine().seed(2)
    print("Direct engine call:", \
          r.get_engine().random(), \
          r.get_engine().random())
    r.set_seed(2)
    print("Wrapped call:", r.random(), r.random())

    print "Test seeding:"
    b = np.sqrt(1.0 - np.power(((5.9 - 0.5)/511. + 1),-2.))
    r.set_seed(666) # diabolic seed
    print r.random()
    print r.photoelectron_theta(b)
    print r.photoelectron_theta(b)
    r.set_seed(666) 
    print r.random()
    print r.photoelectron_theta(b)
    print r.photoelectron_theta(b)

    print "Test state dump:"
    currentstate = r.get_state()
    #print ("Get state", currentstate)
    print("NextRandom", r.random(), r.photoelectron_theta(b))
    r.set_state(currentstate)
    print("PrevRandom", r.random(), r.photoelectron_theta(b))

        
    print "Test exp:", r.exp(2, 10)
    print "Test uniform", r.uniform(2,5.2, 10)
    

def test_univariate(num_events=100000, num_bins=100):
    """
    """
    rv = np.linspace(0, 2*np.pi, 100)
    pdf = np.cos(rv)**2/np.pi
    generator = xpeUnivariateGenerator(rv, pdf)
    values = generator.rvs(num_events)
    h = plt.hist(values, bins=num_bins)
    bin_width = np.pi*2/num_bins
    plt.plot(rv, pdf*num_events*bin_width)
    plt.show()

def compare_with_root(N = 1000, energy = 5.9):
    r = xperandom()
    b = np.sqrt(1.0 - np.power(((energy - 0.5)/511. + 1),-2.))
    # beta da 0.02 a 0.2

    # test vs ROOT
    import ROOT
    ThetaDist = ROOT.TF1("ThetaDist","[0]*pow(sin(x), 3.0)/pow(1.0-[1]*cos(x), 4.0)", 0.0, ROOT.TMath.Pi());
    ThetaDist.SetParameter(0, 1.0);

    h1 = ROOT.TH1F("h1", "h ROOT", 100,0,ROOT.TMath.Pi())
    h2 = ROOT.TH1F("h2", "h np", 100,0,ROOT.TMath.Pi())
    h2.SetLineColor(2)
    print "Extract %d evts" %N
    from time import time
    
    # extract ROOT
    t0 = time()
    for i in range(N):
        ThetaDist.SetParameter(1, b);
        h1.Fill(ThetaDist.GetRandom())
    t0 = time()-t0
    print "root time", t0
    
    # exctract np
    t0 = time()
    for i in range(N):
        h2.Fill(r.photoelectron_theta(b))
    t0 = time()-t0
    print "numpy time", t0
    
    ccc = ROOT.TCanvas()
    h1.Draw()
    h2.Draw("sames")
    raw_input("Enter to close canvas")


if __name__ == '__main__':
    test_random()
    #test_univariate()
    
    
