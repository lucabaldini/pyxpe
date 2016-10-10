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

from pyxpe.recon.binio import xpeBinaryFileWindowed
from pyxpe.recon.clustering import hierarchical_clustering
from pyxpe.recon.geometry import xpeRay2d, xpePoint2d
from pyxpe.recon.recon import xpePixyRecon
from pyxpe.utils.logging_ import logger


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

def display_recon(event, zero_suppression=9, coordinate_system='xpe'):
    """
    """
    cluster = hierarchical_clustering(event, zero_suppression,
                                      coordinate_system)[0]
    logger.info(event)
    logger.info(cluster)

    event_fig = event.draw(show=False, grids=False)
    plt.savefig('sample_evt_raw.pdf')

    cluster_fig = cluster.draw(coordinate_system, show=False)
    cluster.baricenter.draw(color='blue')
    annotate('Baricenter', cluster.baricenter, (0.15, 0.75))
    plt.savefig('sample_evt_cluster.pdf')

    recon = xpePixyRecon(cluster)
    ma0 = recon.ma0
    ma1 = recon.ma1

    mom_fig = cluster.draw(coordinate_system, hexcol_padding=0.1, show=False)
    ma0.draw(color='blue', semiaxes=True)
    annotate('Baricenter', cluster.baricenter, (0.15, 0.75))
    p = ma0.axis.at(0.4)
    annotate('Principal axis', p, (0.45, 0.95))
    p = ma0.major_semiaxis.at(-0.5*numpy.sqrt(ma0.mom2_long))
    annotate('$\\sqrt{M_2^{\\rm long}}$', p, (0.05, 0.5))
    p = ma0.minor_semiaxis.at(-0.5*numpy.sqrt(ma0.mom2_trans))
    annotate('$\\sqrt{M_2^{\\rm trans}}$', p, (0.92, 0.6))
    p = ma0.minor_semiaxis.at(-numpy.sqrt(ma0.mom2_trans))
    annotate('Ellipsoid of inertia', p, (0.8, 0.2))
    plt.savefig('sample_evt_mom_analysis.pdf')

    mom3_fig = plt.figure(facecolor='w')
    ax = plt.subplot(111)
    x = cluster.projection1d(cluster.baricenter, ma0.phi)
    ax.bar(x, cluster.adc_values, width=0.002, color='black')
    plt.xlabel('Projection along the principal axis [mm]')
    plt.ylabel('Pulse height [ADC counts]')
    x3 = numpy.sign(recon.mom3_long)*abs(recon.mom3_long)**(1./3.)
    annotate('Baricenter', (0., 0.), (0.35, 0.85))
    annotate('$\sqrt[3]{M_3}$', (x3, 0.), (0.75, 0.5))
    plt.savefig('sample_evt_mom3.pdf')

    conv_fig = cluster.draw(coordinate_system, hexcol_padding=0.75, show=False)
    _color = 'blue'
    _lw = 1.5
    cluster.baricenter.draw(color=_color)
    annotate('Baricenter', cluster.baricenter, (0.1, 0.65))
    ma0.axis.draw(color=_color, ls='dashed')
    p = ma0.axis.at(0.7)
    annotate('Principal axis', p, (0.45, 0.95))
    ma0.minor_semiaxis.draw(color=_color, lw=_lw)
    r1 = 1.5*numpy.sqrt(ma0.mom2_long)
    c1 = Circle(xy=cluster.baricenter, radius=r1, facecolor='none',
                edgecolor=_color, lw=_lw, hatch='///')
    plt.gca().add_artist(c1)
    p = cluster.baricenter + xpePoint2d(r1, 0)
    annotate('$r_1 = 1.5 \\times \\sqrt{M_2^{\\rm long}}$', p, (0.875, 0.55))
    r2 = 3.5*numpy.sqrt(ma0.mom2_long)
    c2 = Circle(xy=cluster.baricenter, radius=r2, facecolor='none',
                edgecolor=_color, lw=_lw)
    plt.gca().add_artist(c2)
    p = cluster.baricenter + xpePoint2d(r2*numpy.sin(2.5), -r2*numpy.cos(2.5))
    annotate('$r_2 = 3.5 \\times \\sqrt{M_2^{\\rm long}}$', p, (0.9, 0.9))
    _phi1 = numpy.degrees(ma0.phi) + 90
    _phi2 = _phi1 - 180.
    w2 = Wedge(cluster.baricenter, r2, _phi1, _phi2, facecolor='none',
               edgecolor=_color, lw=_lw, hatch='///')
    plt.gca().add_artist(w2)
    recon.conversion_point.draw(color='green')
    annotate('Conversion point', recon.conversion_point, (0.15, 0.9))
    plt.savefig('sample_evt_conv_point.pdf')


    dir2_fig = cluster.draw(coordinate_system, hexcol_padding=0.1, show=False)
    cluster.baricenter.draw(color='blue')
    annotate('Baricenter', cluster.baricenter, (0.1, 0.65))
    ma0.axis.draw(color='blue', lw=_lw, ls='dashed')
    p = ma0.axis.at(0.4)
    annotate('Principal axis', p, (0.45, 0.95))    
    ma1.draw(color='green', semiaxes=True)
    annotate('Conversion point', ma1.pivot, (0.15, 0.9))
    p = ma1.axis.at(0.5)
    annotate('Final direction', p, (0.9, 0.35))
    plt.savefig('sample_evt_phi2.pdf')
    plt.show()
    

if __name__ == '__main__':
    file_path = '/data/work/xpe/xpedaq/data/test_fe_500evts.mdat'
    event = xpeBinaryFileWindowed(file_path).next()
    display_recon(event)
