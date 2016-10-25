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


import ROOT
import numpy

from matplotlib.patches import Ellipse

from pyxpe.recon.binio import xpeBinaryFileWindowed
from pyxpe.recon.rootio import xpePixyTree, xpeReconTree
from pyxpe.recon.clustering import hierarchical_clustering
from pyxpe.recon.geometry import xpePoint2d, xpeRay2d
from pyxpe.utils.logging_ import logger



class xpeMomentsAnalysis:

    """Small container representing the output of a moments analysis.
    """

    def __init__(self, cluster=None, pivot=None, weights=1.):
        """Constructor.
        """
        self.pivot = None
        self.phi = None
        self.mom2_long = None
        self.mom2_trans = None
        if cluster is not None:
            self.run(cluster, pivot, weights)

    def successful(self):
        """Return true if the moments analysis has been run successfully on a
        cluster.
        """
        return self.phi != None

    def run(self, cluster, pivot=None, weights=1.):
        """Run the actual moments analysis on a cluster.
        
        Args
        ----
        cluster : a pyxpe.recon.clustering.xpeCluster object
            The actual cluster object the moments analysis is run onto.

        pivot : xpePoint2d instance
            The pivot point for the moments analysis (by defaults is the
            baricenter of the cluster).

        weights : array
            A set of weights for the moments analysis (default to 1---these
            are multiplied by the pixel ADC counts).
        """
        if pivot is None:
            pivot = cluster.baricenter
        self.pivot = pivot
        w = cluster.adc_values*weights
        wsum = numpy.sum(w)
        # Calculate the offsets with respect to the pivot.
        dx = (cluster.x - pivot.x())
        dy = (cluster.y - pivot.y())
        # Solve for the angle of the principal axis (note that at this point
        # phi is comprised between -pi/2 and pi/2 and might indicate either
        # the major or the minor axis of the tensor of inertia).
        A = numpy.sum(dx*dy*w)
        B = numpy.sum((dy**2. - dx**2.)*w)
        phi = -0.5*numpy.arctan2(2.*A, B)
        # Rotate by an angle phi and calculate the eigenvalues of the tensor
        # of inertia.
        xp = numpy.cos(phi)*dx + numpy.sin(phi)*dy
        yp = -numpy.sin(phi)*dx + numpy.cos(phi)*dy
        mom2_long = numpy.sum((xp**2.)*w)/wsum
        mom2_trans = numpy.sum((yp**2.)*w)/wsum
        # We want mom2_long to be the highest eigenvalue, so we need to
        # check wheteher we have to swap the eigenvalues, here. Note that
        # at this point phi is still comprised between -pi/2 and pi/2.
        if mom2_long < mom2_trans:
            mom2_long, mom2_trans = mom2_trans, mom2_long
            phi -= 0.5*numpy.pi*numpy.sign(phi)
        # Set the class members.
        self.phi = phi
        self.mom2_long = mom2_long
        self.mom2_trans = mom2_trans

    def draw(self, ellipse=True, semiaxes=False, color='black', linewidth=1.5):
        """Draw the output of the moments analysis.
        """
        import matplotlib.pyplot as plt
        self.pivot.draw(color=color)
        self.axis = xpeRay2d(self.pivot, self.phi)
        self.axis.draw(color=color, linewidth=linewidth)
        if ellipse:
            e = Ellipse(xy=self.pivot, width=2*numpy.sqrt(self.mom2_long),
                        height=2*numpy.sqrt(self.mom2_trans),
                        angle=numpy.degrees(self.phi), facecolor='none',
                        edgecolor=color, linewidth=linewidth)
            plt.gca().add_artist(e)
        if semiaxes:
            self.major_semiaxis = xpeRay2d(self.pivot, self.phi)
            self.major_semiaxis.draw(r=numpy.sqrt(self.mom2_long), color=color,
                                     linewidth=linewidth)
            self.minor_semiaxis = xpeRay2d(self.pivot, self.phi + 0.5*numpy.pi)
            self.minor_semiaxis.draw(r=numpy.sqrt(self.mom2_trans), color=color,
                                     linewidth=linewidth)

    def __str__(self):
        """Terminal formatting.
        """
        return 'Moments analysis: phi = %.3f rad, moments = %.3e/%.3e' %\
            (self.phi, self.mom2_long, self.mom2_trans)

    

class xpePixyRecon:

    """Container class for the original Pixy reconstruction.
    """
    CLUSTER_TOO_SMALL = 0x1
    CLUSTER_TOO_LARGE = 0x2
    FIRST_MOMENTS_ANALYSIS_FAILED = 0x4
    
    def __init__(self, cluster=None, min_cluster_size=6, max_cluster_size=400,
                 small_radius=1.5, large_radius=3.5, weight_scale=0.05):
        """Constructor.
        """
        self.error_summary = 0
        if cluster is not None:
            self.run(cluster, min_cluster_size, max_cluster_size, small_radius,
                     large_radius, weight_scale)

    def run(self, cluster, min_cluster_size, max_cluster_size, small_radius,
            large_radius, weight_scale):
        """Run the actual reconstruction, as implemented in Pixy.
        """
        # If the cluster is too small or too big, return.
        if cluster.num_pixels() < min_cluster_size:
            self.error_summary += self.CLUSTER_TOO_SMALL
            return
        if cluster.num_pixels() >= max_cluster_size:
            self.error_summary += self.CLUSTER_TOO_LARGE
            return
        # Run the first-pass moments analysis.
        self.ma0 = xpeMomentsAnalysis(cluster, cluster.baricenter, 1.)
        if not self.ma0.successful():
            self.error_summary += self.FIRST_MOMENTS_ANALYSIS_FAILED
            return
        self.phi0 = self.phi1 = self.ma0.phi
        self.conversion_point = cluster.baricenter
        self.conversion_baricenter = cluster.baricenter
        # Calculate the third moment along the principal axis of the charge
        # distribution.
        self.mom3_long = cluster.moment(3, cluster.baricenter, self.phi0)
        if self.mom3_long == 0:
            return
        # Calculate the distances from the baricenter in units of the
        # longitudinal rms of the charge distribution.
        dx = (cluster.x - cluster.baricenter.x())
        dy = (cluster.y - cluster.baricenter.y())
        d = numpy.sqrt(dx**2 + dy**2)/numpy.sqrt(self.ma0.mom2_long)
        # Calculate the projection of the pixels along the major axis.
        xp = numpy.cos(self.phi0)*dx + numpy.sin(self.phi0)*dy
        # Select all the pixels whose distance from the baricenter is
        # comprised withing the two radii (small and large) and are lying
        # on the "right" side (i.e., that indicated by the third moment).
        _mask = (d > small_radius)*(d < large_radius)*(xp/self.mom3_long > 0.)
        _adc = cluster.adc_values[_mask]
        _adc_sum = float(numpy.sum(_adc))
        # Calculate the center of mass of the selected pixels---this is the
        # reconstructed conversion point.
        if _adc_sum == 0:
            return 
        _x = numpy.sum(cluster.x[_mask]*_adc)/_adc_sum
        _y = numpy.sum(cluster.y[_mask]*_adc)/_adc_sum
        self.conversion_point = xpePoint2d(_x, _y)
        # And now we can assign a direction to the original axis, based on the
        # sign of the third moment. (Note that we could have done this right
        # at the beginning, but that would have changed the logic of the
        # conversion point calculation, and this implementation is adherent
        # to the original code in Pixy).
        if self.mom3_long > 0:
            self.phi0 -= numpy.pi*numpy.sign(self.phi0)
        # And now the second moments analysis (need to wrap all this into a
        # single call to the do_moments_analysis() method).
        dx = (cluster.x - self.conversion_point.x())
        dy = (cluster.y - self.conversion_point.y())
        d = numpy.sqrt(dx**2 + dy**2)
        w = numpy.exp(-d/weight_scale)
        _adc_sum = float(numpy.sum(cluster.adc_values*w))
        _x = numpy.sum(cluster.x*cluster.adc_values*w)/_adc_sum
        _y = numpy.sum(cluster.y*cluster.adc_values*w)/_adc_sum
        self.conversion_baricenter = xpePoint2d(_x, _y)
        self.ma1 = xpeMomentsAnalysis(cluster, self.conversion_baricenter, w)
        self.phi1 = self.ma1.phi
        if abs(self.phi1 - self.phi0) > 0.5*numpy.pi:
            self.phi1 -= numpy.pi*numpy.sign(self.phi1)



class xpeRecon(xpePixyRecon):

    """Subclass for playing around with the recon stuff.
    """

    pass



def run_pixy_recon(file_path, num_events=1000000000, zero_suppression=9,
                   coordinate_system='pixy', output_path=None,
                   min_cluster_size=6, max_cluster_size=400,
                   small_radius=1.5, large_radius=3.5, weight_scale=0.05):
    """Run the event reconstruction on a binary file.
    """
    assert(file_path.endswith('.mdat'))
    if output_path is None:
        output_path = file_path.replace('.mdat', '_pixy.root')
    logger.info('Opening output file %s...' % output_path)
    output_file = ROOT.TFile(output_path, 'RECREATE')
    output_tree = xpePixyTree()
    event_id = 0
    for event in xpeBinaryFileWindowed(file_path):
        cluster_list = hierarchical_clustering(event, zero_suppression,
                                               coordinate_system)
        cluster = cluster_list[0]
        recon = xpePixyRecon(cluster)
        if not recon.error_summary:
            output_tree.set_value('fRunId', -1)
            output_tree.set_value('fEventId', event_id)
            output_tree.set_value('fNClusters', len(cluster_list))
            output_tree.set_value('fTrigWindow', event.num_pixels())
            output_tree.set_value('fTimeTick', -1)
            output_tree.set_value('fTimeStamp', -1)
            output_tree.set_value('fBufferId', event.buffer_id)
            output_tree.set_value('fCluSize', cluster.num_pixels())
            output_tree.set_value('fPHeight', cluster.pulse_height)
            output_tree.set_value('fStoN', -1)
            output_tree.set_value('fTotNoise', -1)
            output_tree.set_value('fBaricenterX', cluster.baricenter.x())
            output_tree.set_value('fBaricenterY', cluster.baricenter.y())
            output_tree.set_value('fTheta0', recon.phi0)
            output_tree.set_value('fTheta1', recon.phi1)
            output_tree.set_value('fMomX', recon.ma0.mom2_long)
            output_tree.set_value('fMomY', recon.ma0.mom2_trans)
            output_tree.set_value('fMomThirdX', recon.mom3_long)
            output_tree.set_value('fImpactX', recon.conversion_point.x())
            output_tree.set_value('fImpactY', recon.conversion_point.y())
            output_tree.Fill()
        event_id += 1
        if event_id == num_events:
            break
    num_proc_events = output_tree.GetEntries()
    output_tree.Write()
    output_file.Close()
    logger.info('Done, %d event(s) written to the output file.' %\
                num_proc_events)
    return output_path
        

def run_xpe_recon(file_path, num_events=1000000000, zero_suppression=9,
                  coordinate_system='pixy', output_path=None,
                  min_cluster_size=6, max_cluster_size=400,
                  small_radius=1.5, large_radius=3.5, weight_scale=0.05):
    """Run the event reconstruction on a binary file.
    """
    assert(file_path.endswith('.mdat'))
    if output_path is None:
        output_path = file_path.replace('.mdat', '_xpe.root')
    logger.info('Opening output file %s...' % output_path)
    output_file = ROOT.TFile(output_path, 'RECREATE')
    output_tree = xpeReconTree()
    event_id = 0
    for event in xpeBinaryFileWindowed(file_path):
        cluster_list = hierarchical_clustering(event, zero_suppression,
                                               coordinate_system)
        cluster = cluster_list[0]
        recon = xpePixyRecon(cluster)
        if not recon.error_summary:
            output_tree.set_value('fRunId', -1)
            output_tree.set_value('fEventId', event_id)
            output_tree.set_value('fNClusters', len(cluster_list))
            output_tree.set_value('fTrigWindow', event.num_pixels())
            output_tree.set_value('fTimeTick', -1)
            output_tree.set_value('fTimeStamp', -1)
            output_tree.set_value('fBufferId', event.buffer_id)
            output_tree.set_value('fCluSize', cluster.num_pixels())
            output_tree.set_value('fPHeight', cluster.pulse_height)
            output_tree.set_value('fStoN', -1)
            output_tree.set_value('fTotNoise', -1)
            output_tree.set_value('fBaricenterX', cluster.baricenter.x())
            output_tree.set_value('fBaricenterY', cluster.baricenter.y())
            output_tree.set_value('fTheta0', recon.phi0)
            output_tree.set_value('fTheta1', recon.phi1)
            output_tree.set_value('fMomX', recon.ma0.mom2_long)
            output_tree.set_value('fMomY', recon.ma0.mom2_trans)
            output_tree.set_value('fMomThirdX', recon.mom3_long)
            output_tree.set_value('fImpactX', recon.conversion_baricenter.x())
            output_tree.set_value('fImpactY', recon.conversion_baricenter.y())
            output_tree.Fill()
        event_id += 1
        if event_id == num_events:
            break
    num_proc_events = output_tree.GetEntries()
    output_tree.Write()
    output_file.Close()
    logger.info('Done, %d event(s) written to the output file.' %\
                num_proc_events)
    return output_path


if __name__ == '__main__':
    file_path = '/data/work/xpe/xpedaq/data/test_fe_500evts.mdat'
    for event in xpeBinaryFileWindowed(file_path):
        cluster_list = hierarchical_clustering(event, 9, 'xpe')
        cluster = cluster_list[0]
        recon = xpePixyRecon(cluster)
        print cluster.pulse_height, cluster.num_pixels(),\
            cluster.baricenter.x(), cluster.baricenter.y(), recon.phi0,\
            recon.ma0.mom2_long, recon.ma0.mom2_trans
        raw_input()
