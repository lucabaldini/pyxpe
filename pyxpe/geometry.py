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



class xpePoint2d(numpy.ndarray):

    """Small numpy ndarray subclass representing a point in two dimensions.
    """

    def __new__(cls, x, y, units='mm'):
        """Look here
        http://docs.scipy.org/doc/numpy-1.10.1/user/basics.subclassing.html
        as to way we need __new__, as opposed to __init__, here.
        """
        obj = numpy.ndarray.__new__(cls, shape=(2,), dtype='d',
                                    buffer=numpy.array([x, y]))
        obj.units = units
        return obj

    def __array_finalize__(self, obj):
        """Again, look at
        http://docs.scipy.org/doc/numpy-1.10.1/user/basics.subclassing.html
        """
        if obj is None:
            return
        self.units = getattr(obj, 'units', None)

    def x(self):
        """Return the first element of the array (i.e., the x coordinate).
        """
        return self[0]
        
    def y(self):
        """Return the second element of the array (i.e., the x coordinate).
        """
        return self[1]

    def draw(self, **kwargs):
        """
        """
        kwargs['s'] = kwargs.get('s', None) or 75
        plt.scatter(self.x(), self.y(), **kwargs)

    def __str__(self):
        """String formatting
        """
        return '(%.3f, %.3f) %s' % (self.x(), self.y(), self.units)
    


class xpeRay2d:

    """Class representing a 2-dimensional ray.
    """

    def __init__(self, p0, phi):
        """Constructor.
        """
        self.p0 = p0
        self.phi = phi

    def at(self, d):
        """
        """
        _x = self.p0.x() + d*numpy.cos(self.phi)
        _y = self.p0.y() + d*numpy.sin(self.phi)
        return xpePoint2d(_x, _y)

    def draw(self, r=1., **kwargs):
        """
        """
        x1 = self.p0.x() + r*numpy.cos(self.phi)
        x2 = self.p0.x() - r*numpy.cos(self.phi)
        y1 = self.p0.y() + r*numpy.sin(self.phi)
        y2 = self.p0.y() - r*numpy.sin(self.phi)
        plt.plot([x1, x2], [y1, y2], 'k-', **kwargs)
