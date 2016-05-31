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

from pyxpe.xpol import XPOL_COLUMN_PITCH
from pyxpe.geometry import xpePoint2d, xpeRay2d
from pyxpe.xpol import xpeHexagonCollection, adc2colors



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
        self.baricenter = xpePoint2d(_x, _y)
        # Run a first moments analysis with all the pixels.
        if self.num_pixels() > 2:
            self.phi0, self.mom2_long,\
                self.mom2_trans = self.do_moments_analysis(self.baricenter)
            self.axis = xpeRay2d(self.baricenter, self.phi0)
            self.mom3_principal = self.mom3(self.baricenter, self.phi0)

    def num_pixels(self):
        """Return the cluster size.
        """
        return len(self.adc_values)

    def __cmp__(self, other):
        """Comparison operator (sort the clusters by pulse height).
        """
        return other.pulse_height - self.pulse_height

    def do_moments_analysis(self, pivot, weights=None):
        """Run a two-dimensional moments analysis on the cluster.
        
        Args
        ----
        pivot : xpePoint2d instance
            The pivot point for the moments analysis (by defaults is the
            baricenter of the cluster).

        weights : array
            A set of weights for the moments analysis.

        Warning
        -------
        Shall we be using atan2, here?
        """
        # Calculate the offsets with respect to the pivot.
        dx = (self.x - pivot.x())
        dy = (self.y - pivot.y())
        # Solve for the angle of the principal axis.
        num = 2*numpy.sum(dx*dy*self.adc_values)
        den = numpy.sum((dy**2. - dx**2.)*self.adc_values)
        phi = -0.5*numpy.arctan(num/den)
        # Calculate the eigenvalues of the inertia tensor.
        dxp = numpy.cos(phi)*dx + numpy.sin(phi)*dy
        dyp = -numpy.sin(phi)*dx + numpy.cos(phi)*dy
        mom2_long = numpy.sum(dxp**2*self.adc_values)/self.pulse_height
        mom2_trans = numpy.sum(dyp**2*self.adc_values)/self.pulse_height
        # Need to swap the eigrnvalues?
        if mom2_trans > mom2_long:
            mom2_long, mom2_trans = mom2_trans, mom2_long
            phi += 0.5*numpy.pi
        # Return the results of the moments analysis.
        return phi, mom2_long, mom2_trans

    def projection1d(self, pivot, phi):
        """Project the charge distribution on the ray passing by the pivot
        point at an angle phi, and return the corresponding one-dimensional
        array of coordinates.
        """
        return numpy.cos(phi)*(self.x - pivot.x()) +\
            numpy.sin(phi)*(self.y - pivot.y())

    def mom3(self, pivot, phi):
        """
        """
        xp = self.projection1d(pivot, phi)
        return numpy.sum((xp**3)*self.adc_values)/self.pulse_height

    def fit_spline(self):
        """
        """
        from scipy.interpolate import UnivariateSpline
        weights = (self.adc_values/float(self.adc_values.max()))**0.5
        s = UnivariateSpline(self.x, self.y, w=weights, s=0.5)
        _x = numpy.linspace(self.x.min(), self.x.max(), 25)
        _y = s(_x)
        plt.plot(_x, _y, '-', lw=2, color='black')

    def draw(self, coordinate_system, color_map='Reds', text=True, show=True):
        """Draw the cluster.
        """
        hit_positions = numpy.vstack((self.x, self.y),).transpose()
        colors = adc2colors(self.adc_values, 0, color_map)
        if coordinate_system == 'pixy':
            angle = numpy.pi/2.
        else:
            angle = 0
        hex_col = xpeHexagonCollection(offsets=hit_positions, rotation=angle,
                                       edgecolors='gray', facecolors=colors)
        fig = hex_col.figure
        if text:
            adc_ref = 0.5*self.adc_values.max()
            for x, y, val in zip(self.x, self.y, self.adc_values):
                if val < adc_ref:
                    col = 'black'
                else:
                    col = 'white'
                plt.text(x, y, '%s' % val, horizontalalignment='center',
                         verticalalignment='center', size=8, color=col)
        plt.xlabel('x [mm]')
        plt.ylabel('y [mm]')
        self.baricenter.draw()
        self.axis.draw()
        self.fit_spline()
        if show:
            plt.show()

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
                            max_distance=1.01*XPOL_COLUMN_PITCH):
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



def test(filePath, num_events, zero_suppression=9, coordinate_system='xpedaq'):
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
        print cluster.phi0, cluster.mom2_long, cluster.mom2_trans,\
            cluster.mom3_principal
        cluster.draw(coordinate_system)


        
if __name__ == '__main__':
    import argparse
    formatter = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=formatter)
    parser.add_argument('infile', type=str,
                        help='the input binary file')
    parser.add_argument('-n', '--num_events', type=int, default=10,
                        help = 'number of events to be processed')
    parser.add_argument('-z', '--zero-suppression', type=int, default=9,
                        help = 'zero-suppression threshold')
    parser.add_argument('-c', '--coordinate-system', type=str, default='pixy',
                        help = 'coordinate system for the clustering')
    args = parser.parse_args()
    test(args.infile, args.num_events, args.zero_suppression,
         args.coordinate_system)
