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

from matplotlib import collections, transforms
from pyxpe.xpol import XPOL_COLUMN_PITCH



class xpe2dPoint(numpy.ndarray):

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

    def draw(self, show=True):
        """
        """
        plt.plot(self.x(), self.y(), 'o')
        if show:
            plt.show()

    def __str__(self):
        """String formatting
        """
        return '(%.3f, %.3f) %s' % (self.x(), self.y(), self.units)
    


class xpeCluster:

    """Class representing a cluster of pixels, i.e., a photoelectron track.
    """

    def __init__(self, x, y, adc_values):
        """Constructor.
        """
        assert len(x) == len(y) == len(adc_values)
        self.x = x
        self.y = y
        self.adc_values = adc_values
        # Calculate the pulse height.
        self.pulse_height = self.adc_values.sum()
        # Calculate the baricenter position.
        _x = numpy.sum(self.x*self.adc_values)/self.pulse_height
        _y = numpy.sum(self.y*self.adc_values)/self.pulse_height
        self.baricenter = xpe2dPoint(_x, _y)
        self.__do_moments_analysis()

    def num_pixels(self):
        """Return the cluster size.
        """
        return len(self.adc_values)

    def __cmp__(self, other):
        """
        """
        return other.pulse_height - self.pulse_height

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

    def __str__(self):
        """String formatting.
        """
        return 'Cluster @ %s, %d pixels, pulse height = %d ADC counts' %\
            (self.baricenter, self.num_pixels(), self.pulse_height)


    
def single_clustering(event, zero_suppression, coordinate_system):
    """Dummy single-clustering algorithm for testing purposes.

    This takes all the pixels above the zero-suppression threshold in the window
    and returns the corresponing cluster object.

    Args
    ----
    event : pXpeEventBase instance
        The underlying event object.

    zero_suppression : float or array
        The zero suppression threshold.
    """
    x, y, adc_values = event.hit_data(zero_suppression, coordinate_system)
    return [xpeCluster(x, y, adc_values)]



def hierarchical_clustering(event, zero_suppression, coordinate_system,
                            method='single', metric='euclidean',
                            criterion='distance',
                            max_distance=1.1*XPOL_COLUMN_PITCH):
    """Lightweight wrapper over the scipy.cluster.hierarchy module.

    This is essentially calling scipy.cluster.hierarchy.linkage and
    scipy.cluster.hierarchy.fcluster, returning a list of xpeCluster objects
    sorted by pulse height.

    The default parameters in the method signature are those producing the
    behaviour you would naively expect, i.e., contiguous pixels are gruped
    together in the same cluster (the clustering is done using the
    pixel-to-pixel euclidean distance and the distance cut is placed just above
    the longest among the readout pitches in the two directions).

    Warning
    -------
    Do not mess around with the function arguments unless you know what you're
    doing---scipy is capable of producing a surprising variety of different
    behaviours, and the problem we're trying to solve here is fairly simple.

    Args
    ----
    event : pXpeEventBase instance
        The underlying event object.

    zero_suppression : float or array
        The zero suppression threshold.

    method : str (default 'single')
        The clustering method passed to scipy.cluster.hierarchy.linkage

    metric : str (default 'euclidean')
        The metric passed to scipy.cluster.hierarchy.linkage

    criterion: str (default 'distance')
        The criterion passed to scipy.cluster.hierarchy.fcluster

    max_distance : float
        The maximum distance used by scipy.cluster.hierarchy.fcluster

    Return
    ------
    a list of xpeCluster objects sorted by pulse height.
    """
    import scipy.cluster.hierarchy
    x, y, adc_values = event.hit_data(zero_suppression, coordinate_system)
    data = numpy.vstack((x, y),).transpose()
    Z = scipy.cluster.hierarchy.linkage(data, method, metric)
    cluster_ids = scipy.cluster.hierarchy.fcluster(Z, max_distance, criterion)
    cluster_list = []
    for i in xrange(1, max(cluster_ids) + 1):
        _mask = numpy.where(cluster_ids == i)
        cluster_list.append(xpeCluster(x[_mask], y[_mask], adc_values[_mask]))
    cluster_list.sort()
    return cluster_list



def test(filePath, num_events, zero_suppression=9, coordinate_system='pixy'):
    """
    """
    from pyxpe.binio import xpeBinaryFileWindowed
    input_file = xpeBinaryFileWindowed(filePath)
    for i in xrange(num_events):
        event = input_file.next()
        print event
        cluster_list = hierarchical_clustering(event, zero_suppression,
                                               coordinate_system)
        cluster = cluster_list[0]
        print cluster
        print cluster.phi0, cluster.mom2_long, cluster.mom2_trans
        event.draw(zero_suppression, show=False)
        cluster.baricenter.draw()


        
if __name__ == '__main__':
    import argparse
    formatter = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=formatter)
    parser.add_argument('infile', type=str,
                        help='the input binary file')
    parser.add_argument('--num_events', type=int, default=9,
                        help = 'number of events to be processed')
    args = parser.parse_args()
    test(args.infile, args.num_events)
