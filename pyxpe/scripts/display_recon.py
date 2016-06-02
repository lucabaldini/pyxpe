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
from matplotlib.patches import Ellipse, Circle, Wedge

from pyxpe.binio import xpeBinaryFileWindowed
from pyxpe.clustering import hierarchical_clustering
from pyxpe.geometry import xpeRay2d
from pyxpe.logging_ import logger


def annotate(text, pos, text_pos, text_size=15, color='gray',
             bbstyle='roundtooth'):
    """Annotate the current figure.
    """
    ax = plt.gca()
    _bbox = dict(boxstyle=bbstyle, color=color, fc='w')
    _arrowprops = dict(arrowstyle="->", color=color, 
                       connectionstyle='angle3,angleA=0,angleB=90')
    ax.annotate(text, xy=pos, xycoords='data', xytext=text_pos,
                textcoords='figure fraction', size=text_size, va='center',
                ha='center', color=color, bbox=_bbox,
                arrowprops=_arrowprops)
    

def display_recon(event, zero_suppression=9, coordinate_system='pixy'):
    """
    """
    cluster = hierarchical_clustering(event, zero_suppression,
                                      coordinate_system)[0]
    logger.info(event)
    logger.info(cluster)

    event_fig = event.draw(show=False)

    cluster_fig = cluster.draw(coordinate_system, show=False)

    mom_fig = cluster.draw(coordinate_system, hexcol_padding=0.1, show=False)
    _color = 'blue'
    _lw = 1.5
    cluster.baricenter.draw(color=_color)
    annotate('Baricenter', cluster.baricenter, (0.8, 0.85))
    cluster.axis0.draw(color=_color, lw=_lw, ls='dashed')
    p = cluster.axis0.at(-0.4)
    annotate('Principal axis', p, (0.1, 0.9))
    major_axis = xpeRay2d(cluster.baricenter, cluster.phi0)
    major_axis.draw(r=numpy.sqrt(cluster.mom2_long), color=_color, lw=_lw)
    minor_axis = xpeRay2d(cluster.baricenter, cluster.phi0 + 0.5*numpy.pi)
    minor_axis.draw(r=numpy.sqrt(cluster.mom2_trans), color=_color, lw=_lw)
    e = Ellipse(xy=cluster.baricenter, width=2*numpy.sqrt(cluster.mom2_long),
                height=2*numpy.sqrt(cluster.mom2_trans),
                angle=numpy.degrees(cluster.phi0), facecolor='none',
                edgecolor=_color, lw=_lw)
    plt.gca().add_artist(e)
    p = minor_axis.at(-numpy.sqrt(cluster.mom2_trans))
    annotate('Ellipsoid of inertia', p, (0.5, 0.1))

    mom3_fig = plt.figure(facecolor='w')
    ax = plt.subplot(111)
    x = cluster.projection1d(cluster.baricenter, cluster.phi0)
    ax.bar(x, cluster.adc_values, width=0.002, color='black')
    plt.xlabel('Projection along the principal axis [mm]')
    plt.ylabel('Pulse height [ADC counts]')
    x3 = numpy.sign(cluster.mom3_long)*abs(cluster.mom3_long)**(1./3.)
    annotate('Baricenter', (0., 0.), (0.65, 0.8))
    annotate('$\sqrt[3]{M_3}$', (x3, 0.), (-0.75, 0.5))

    conv_fig = cluster.draw(coordinate_system, hexcol_padding=0.75, show=False)
    _color = 'blue'
    _lw = 1.5
    cluster.baricenter.draw(color=_color)
    annotate('Baricenter', cluster.baricenter, (0.8, 0.85))
    cluster.axis0.draw(color=_color, lw=_lw, ls='dashed')
    p = cluster.axis0.at(-0.8)
    annotate('Principal axis', p, (-0.9, 0.9))
    minor_axis = xpeRay2d(cluster.baricenter, cluster.phi0 + 0.5*numpy.pi)
    minor_axis.draw(color=_color, lw=_lw)
    r1 = 1.5*numpy.sqrt(cluster.mom2_long)
    c1 = Circle(xy=cluster.baricenter, radius=r1, facecolor='none',
                edgecolor=_color, lw=_lw, hatch='///')
    plt.gca().add_artist(c1)
    r2 = 3.5*numpy.sqrt(cluster.mom2_long)
    c2 = Circle(xy=cluster.baricenter, radius=r2, facecolor='none',
                edgecolor=_color, lw=_lw)
    plt.gca().add_artist(c2)
    _phi1 = numpy.degrees(cluster.phi0) - 90
    _phi2 = _phi1 + 180.
    w2 = Wedge(cluster.baricenter, r2, _phi1, _phi2, facecolor='none',
               edgecolor=_color, lw=_lw, hatch='///')
    plt.gca().add_artist(w2)
    cluster.conversion_point.draw()
    annotate('Conversion point', cluster.conversion_point, (0.1, 0.1))

    dir2_fig = cluster.draw(coordinate_system, hexcol_padding=0.1, show=False)
    
    
    plt.show()


if __name__ == '__main__':
    file_path = '/data/work/xpe/xpedaq/data/test_fe_500evts.mdat'
    event = xpeBinaryFileWindowed(file_path).next()
    display_recon(event)

