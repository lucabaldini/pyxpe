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

from pyxpe.recon.binio import xpeBinaryFileFullFrame
from pyxpe.utils.logging_ import logger
from pyxpe.utils.profile import xpeChrono
from pyxpe.recon.xpol import XPOL_NUM_PIXELS, XPOL_NUM_COLUMNS, XPOL_NUM_ROWS


def mean_rms(file_path, num_events=1000000, plot=False):
    """Step through a full-frame binary file and calculate the
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
    n = 0
    mean = numpy.zeros((XPOL_NUM_COLUMNS, XPOL_NUM_ROWS), 'd')
    M2 = numpy.zeros((XPOL_NUM_COLUMNS, XPOL_NUM_ROWS), 'd')
    for adc_values in xpeBinaryFileFullFrame(file_path):
        n += 1
        delta = adc_values - mean
        mean += delta/n
        M2 += delta*(adc_values - mean)
        if n == num_events:
            break
    rms = numpy.sqrt(M2/(n - 1))
    if plot:
        fig_mean = plt.figure('Pedestal mean')
        img_mean = plt.gca().matshow(mean, vmin=750, vmax=1250)
        cax_mean = fig_mean.add_axes([0.9, 0.1, 0.03, 0.8])
        fig_mean.colorbar(img_mean, cax=cax_mean)
        fig_rms = plt.figure('Pedestal RMS')
        img_rms = plt.gca().matshow(rms, vmin=0, vmax=20)
        cax_rms = fig_rms.add_axes([0.9, 0.1, 0.03, 0.8])
        fig_rms.colorbar(img_rms, cax=cax_rms)
        plt.figure()
        plt.hist(mean.flatten(), bins=numpy.arange(0, 1800, 9))
        plt.figure()
        plt.hist(rms.flatten(), bins=numpy.arange(0, 50, 1))
        plt.show()
    return mean, rms


def process_run(input_file_path, output_file_path=None, num_events=1000000):
    """Process a pedestal run in full-frame mode and write a FITS file
    with the pedestal and noise values for each pixel in the matrix.
    """
    assert(input_file_path.endswith('.mdat'))
    if output_file_path is None:
        output_file_path = input_file_path.replace('.mdat', '_peds.fits')
    assert(output_file_path.endswith('.fits'))
    mean, rms = mean_rms(input_file_path, num_events, False)
    from astropy.io import fits
    import time
    primary_hdu = fits.PrimaryHDU()
    primary_hdu.header.set('CREATOR', 'pedestals.py',
                           's/w task which wrote this pixmap')
    _date = time.strftime('%Y-%m-%dT%H:%M:%S', time.gmtime())
    primary_hdu.header.set('DATE', _date,
                           'file creation date (YYYY-MM-DDThh:mm:ss UT')
    primary_hdu.header.set('NCOLS', XPOL_NUM_COLUMNS,
                           'number of columns in the matrix')
    primary_hdu.header.set('NROWS', XPOL_NUM_ROWS,
                           'number of rows in the matrix')
    pedmean_hdu = fits.ImageHDU(mean.transpose(), name='PEDMEAN')
    pedrms_hdu = fits.ImageHDU(rms.transpose(), name='PEDRMS')
    hdu_list = fits.HDUList([primary_hdu, pedmean_hdu, pedrms_hdu])
    hdu_list.writeto(output_file_path, clobber=True)
    logger.info('Pedestals written to %s.' % output_file_path)
    


if __name__ == '__main__':
    import argparse
    formatter = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=formatter)
    parser.add_argument('binfile', type=str,
                        help='the input binary file')
    parser.add_argument('-n', '--num_events', type=int, default=1000000,
                        help = 'number of events to be processed')
    args = parser.parse_args()
    process_run(args.binfile, None, args.num_events)
