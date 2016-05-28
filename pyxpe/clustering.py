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

from pyxpe.event import pXpeBinaryFileWindowed
from pyxpe.xpol import XPOL_MATRIX


class p2dPoint(numpy.ndarray):

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

    def __str__(self):
        """String formatting
        """
        return '(%.3f, %.3f) %s' % (self.x(), self.y(), self.units)
    


class pCluster:

    """Class representing a cluster of pixels, i.e., a photoelectron track.
    """

    def __init__(self, x, y, adc_values):
        """Constructor.
        """
        assert len(x) == len(y) == len(adc_values)
        self.x = x
        self.y = y
        self.adc_values = adc_values
        self.__compute_pulse_height()
        self.__compute_baricenter()
        self.__do_moments_analysis()

    def size(self):
        """Return the cluster size.
        """
        return len(self.adc_values)

    def __compute_pulse_height(self):
        """Calculate the pulse height of the cluster.
        """
        self.pulse_height = self.adc_values.sum()

    def __compute_baricenter(self):
        """Calculate the baricenter of the cluster.
        """
        _x = numpy.sum(self.x*self.adc_values)/self.pulse_height
        _y = numpy.sum(self.y*self.adc_values)/self.pulse_height
        self.baricenter = p2dPoint(_x, _y)

    def __do_moments_analysis(self):
        """Do the first-step moments analysis.
        """
        dx = (self.x - self.baricenter.x())
        dy = (self.y - self.baricenter.y())
        num = 2*numpy.sum(dx*dy*self.adc_values)
        den = numpy.sum((dy**2. - dx**2.)*self.adc_values)
        self.phi0 = -0.5*numpy.arctan(num/den)
        dxp = numpy.cos(self.phi0)*dx + numpy.sin(self.phi0)*dy
        dyp = -numpy.sin(self.phi0)*dx + numpy.cos(self.phi0)*dy
        self.mom2_long = numpy.sum(dxp**2*self.adc_values)/self.pulse_height
        self.mom2_trans = numpy.sum(dyp**2*self.adc_values)/self.pulse_height


    
def single_clustering(event, zero_suppression=10):
    """Dummy single-clustering algorithm for testing purposes.

    This takes all the pixels above the zero-suppression threshold in the window
    and returns the corresponing cluster object.
    """
    adc_values =  event.adc_counts[event.adc_counts >= zero_suppression]
    col, row = numpy.where(event.adc_counts >= zero_suppression)
    x, y = XPOL_MATRIX.pixel2world_recon(event.xmin + col, event.ymin + row)
    return pCluster(x, y, adc_values)



def test(filePath, num_events):
    """
    """
    input_file = pXpeBinaryFileWindowed(filePath)
    for i in xrange(num_events):
        event = input_file.next()
        print event
        cluster = single_clustering(event, 9)
        print cluster.size(), cluster.pulse_height, cluster.baricenter,\
            cluster.phi0, cluster.mom2_long, cluster.mom2_trans
        event.draw()


        
if __name__ == '__main__':
    import argparse
    formatter = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=formatter)
    parser.add_argument('infile', type=str,
                        help='the input binary file')
    parser.add_argument('--num_events', type=int, default=10,
                        help = 'number of events to be processed')
    args = parser.parse_args()
    test(args.infile, args.num_events)
