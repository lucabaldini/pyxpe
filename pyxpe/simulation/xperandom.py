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

# default now is python pseudo-random generator
# https://docs.python.org/2.7/library/random.html
import random

import numpy
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
        _x = numpy.array([self.integral(_xmin, _xp) for _xp in _y])
        _x, _mask = numpy.unique(_x, return_index=True)
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
        return self.ppf(numpy.random.sample(size))

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
        self.engine = random

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
        return self.engine.getstate()

    def set_state(self, state):
        """
        """
        return self.engine.setstate(state)
    
    def exp(self, l):
        """
        """
        return self.engine.expovariate(l)



def test_random():
    """
    """
    r = xperandom()
     
    r.get_engine().seed(2)
    print("Direct engine call:", \
          r.get_engine().random(), \
          r.get_engine().random())
    
    r.set_seed(2)
    print("Wrapped call:", r.random(), r.random())
    
    currentstate = r.get_state()
    #print ("Get state", currentstate)
    print("NextRandom", r.random())
    r.set_state(currentstate)
    print("PrevRandom", r.random())

def test_univariate(num_events=100000, num_bins=100):
    """
    """
    rv = numpy.linspace(0, 2*numpy.pi, 100)
    pdf = numpy.cos(rv)**2/numpy.pi
    generator = xpeUnivariateGenerator(rv, pdf)
    values = generator.rvs(num_events)
    h = plt.hist(values, bins=num_bins)
    bin_width = numpy.pi*2/num_bins
    plt.plot(rv, pdf*num_events*bin_width)
    plt.show()

    
if __name__ == '__main__':
    test_univariate()
    
