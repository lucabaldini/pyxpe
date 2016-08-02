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
from scipy.sparse import lil_matrix

from pyxpe.recon.binio import xpeBinaryFileWindowed
from pyxpe.utils.logging_ import logger
from pyxpe.utils.profile import xpeChrono
from pyxpe.recon.daq import load_run_info
from pyxpe.recon.event import xpeAnsiColors



def compute_mean_rms(file_path, plot=False):
    """Step through a charge-injection binary file and calculate the
    mean and standard deviation of the charge distribution in each pixel.

    Warning
    -------
    Note that, since the signal in ADC counts can be a relatively large
    number, we use Welford's algorithm for the running variance:
    https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    And yes, this is not a nuisance---the first attempt at using the
    dumb running variance ended up in flames (in case you're tempted to
    change this).
    """
    assert(file_path.endswith('.mdat'))
    n = numpy.zeros((300, 352), 'd')
    mean = numpy.zeros((300, 352), 'd')
    M2 = numpy.zeros((300, 352), 'd')
    for event in xpeBinaryFileWindowed(file_path):
        indices = (slice(event.xmin, event.xmax + 1),
                   slice(event.ymin, event.ymax + 1))
        n[indices] += 1
        delta = event.adc_values - mean[indices]
        mean[indices] += delta/n[indices]
        M2[indices] += delta*(event.adc_values - mean[indices])
    rms = numpy.sqrt(M2/(n - 1))
    if plot:
        fig, axes = plt.subplots(nrows=1, ncols=2)
        img0 = axes[0].matshow(mean)
        img1 = axes[1].matshow(rms)
        plt.show()
    return mean, rms


def analyze_run(folder_path):
    """Analyze a charge-injection run.
    """
    run_info = load_run_info(folder_path)
    config = run_info['config']
    x0 = config['pixel_address_x']
    y0 = config['pixel_address_y']
    logger.info('Charge injected at address (%d, %d)...' % (x0, y0))
    logger.info('Readout at %.2f MHz, clock shift %d ns...' %\
                (config.clock_freq_mhz(), config.clock_shift_ns()))
    mean, rms = compute_mean_rms(run_info['data_file_path'])
    padding = 5
    for x in range(x0 - padding, x0 + padding + 1):
        line = '(%3d, %3d) -> %.2f +/- %.2f' % (x, y0, mean[x, y0], rms[x, y0])
        if x == x0:
            line = '%s%s%s' % (xpeAnsiColors.RED, line, xpeAnsiColors.ENDC)
        logger.info(line)

    

if __name__ == '__main__':
    import argparse
    formatter = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=formatter)
    parser.add_argument('folder', type=str,
                        help='the input run folder')
    parser.add_argument('-n', '--num_events', type=int, default=10,
                        help = 'number of events to be processed')
    args = parser.parse_args()
    #mean_variance(args.binfile, args.num_events)
    analyze_run(args.folder)
