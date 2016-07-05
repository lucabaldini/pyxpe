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

from pyxpe.recon.binio import xpeBinaryFileWindowed
from pyxpe.recon.rootio import xpePixyTree
from pyxpe.recon.clustering import hierarchical_clustering
from pyxpe.utils.logging_ import logger


def run_pixy(file_path, num_events=1000000000, zero_suppression=9,
             coordinate_system='xpedaq', output_path=None,
             min_cluster_size=6, max_cluster_size=400):
    """Run the event reconstruction on a binary file.
    """
    assert(file_path.endswith('.mdat'))
    if output_path is None:
        output_path = file_path.replace('.mdat', '.root')
    logger.info('Opening output file %s...' % output_path)
    output_file = ROOT.TFile(output_path, 'RECREATE')
    output_tree = xpePixyTree()
    event_id = 0
    for event in xpeBinaryFileWindowed(file_path):
        cluster_list = hierarchical_clustering(event, zero_suppression,
                                               coordinate_system)
        cluster = cluster_list[0]
        if cluster.num_pixels() >= min_cluster_size and\
           cluster.num_pixels() < max_cluster_size:
            _data = {
                'fRunId': -1,
                'fEventId': event_id, 
                'fNClusters': len(cluster_list), 
                'fTrigWindow': event.num_pixels(), 
                'fTimeTick': -1, 
                'fTimeStamp': -1, 
                'fBufferId': event.buffer_id, 
                'fCluSize': cluster.num_pixels(),
                'fPHeight': cluster.pulse_height,
                'fStoN': -1,
                'fTotNoise': -1,
                'fBaricenterX': cluster.baricenter.x(),
                'fBaricenterY': cluster.baricenter.y(),
                'fTheta0': cluster.phi0,
                'fTheta1': cluster.phi1,
                'fMomX': cluster.mom2_long,
                'fMomY': cluster.mom2_trans,
                'fMomThirdX': cluster.mom3_long,
                'fImpactX': cluster.conversion_point.x(),
                'fImpactY': cluster.conversion_point.y()
            }
            output_tree.fill(_data)
        event_id += 1
        if event_id == num_events:
            break
    num_proc_events = output_tree.GetEntries()
    output_tree.Write()
    output_file.Close()
    logger.info('Done, %d event(s) written to the output file.' %\
                num_proc_events)
    return output_path
        
