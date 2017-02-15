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


from pyxpe.utils.logging_ import logger

import struct
import numpy

import matplotlib
import matplotlib.pyplot as plt

from pyxpe.recon.xpol import xpeHexagonalMatrix, pixel2world

# python2/3 compatibility fix
try:
    xrange
except NameError:
    xrange = range


class xpeAnsiColors:
    
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


    
class xpeEventBase:

    """
    """

    pass



class xpeEventFullFrame(xpeEventBase):

    """Basic class representing an event aquired in full-frame mode.
    """

    def __init__(self, adc_values):
        """Constructor.
        """
        self.adc_values = adc_values
        
    def size(self):
        """Return the total number of bytes in the event.
        """
        return 2*self.num_pixels()
        
    def num_columns(self):
        """Return the number of columns.
        """
        return 300
        
    def num_rows(self):
        """Return the number of rows.
        """
        return 352

    def num_pixels(self):
        """Return the total number of pixels in the window.
        """
        return self.num_rows()*self.num_columns() 
        
    def adc_value(self, col, row):
        """Return the pulse height for a given pixel in the window.
        """
        return self.adc_values[col, row]

    def highest_pixel(self):
        """Return the coordinats of the pixel with the maximum value of
        ADC counts.
        """
        return numpy.unravel_index(numpy.argmax(self.adc_values),
                                   self.adc_values.shape)
    
    def highest_adc_value(self):
        """Return the maximum value of ADC counts for the pixels in the event.
        """
        return self.adc_values.max()

    def draw(self, show = True):
        """
        """        
        im = plt.imshow(self.adc_values.T, cmap='hot', interpolation='none')
        plt.colorbar(im, orientation='horizontal')
        if show:
            plt.show()    
    
    
class xpeEventWindowed(xpeEventBase):
    
    """Basic class representing an event aquired in windowed mode.
    """
    
    HEADER_MARKER = 65535
    HEADER_LENGTH = 20
    
    def __init__(self, xmin, xmax, ymin, ymax, buffer_id, t1, t2, s1, s2,
                 adc_values):
        """Constructor.
        """
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.buffer_id = buffer_id
        self.microseconds = (t1 + t2*65534)*0.8
        self.adc_values = adc_values

    def size(self):
        """Return the total number of bytes in the event.
        """
        return self.HEADER_LENGTH + 2*self.num_pixels()

    def num_columns(self):
        """Return the number of columns.
        """
        return (self.xmax - self.xmin + 1)

    def num_rows(self):
        """Return the number of rows.
        """
        return (self.ymax - self.ymin + 1)

    def num_pixels(self):
        """Return the total number of pixels in the window.
        """
        return self.num_rows()*self.num_columns()
        
    def adc_value(self, col, row):
        """Return the pulse height for a given pixel in the window.
        """
        return self.adc_values[col, row]

    def highest_pixel(self):
        """Return the coordinats of the pixel with the maximum value of
        ADC counts.
        """
        return numpy.unravel_index(numpy.argmax(self.adc_values),
                                   self.adc_values.shape)

    def highest_adc_value(self):
        """Return the maximum value of ADC counts for the pixels in the event.
        """
        return self.adc_values.max()

    def pulse_height(self, zero_suppression):
        """Return the total pulse height for the event, i.e., the raw sum of
        all the ADC values above the zero-suppression threshold.

        Args
        ----
        zero_suppression : int
            The zero-suppression threshold.
        """
        return self.adc_values[self.adc_values > zero_suppression].sum()

    def hit_data(self, zero_suppression, coordinate_system):
        """Return three arrays (of the same length) containing the x and y
        coordinates and the charge in ADC counts for all the pixels in the
        event above the zero-suppression threshold.

        Args
        ----
        zero_suppression : int
            The zero-suppression threshold.

        coordinate_system : str
            The coordinate system to be used.
        """
        _mask = self.adc_values > zero_suppression
        adc_values = self.adc_values[_mask]
        col, row = numpy.where(_mask)
        x, y = pixel2world(self.xmin + col, self.ymin + row, coordinate_system)
        return x, y, adc_values

    def ascii(self, zero_suppression=9, max_threshold=0.75, width=4,
              color=True):
        """Return a pretty-printed ASCII representation of the event.
        """
        _fmt = '%%%dd' % width
        _max = self.highest_adc_value()
        text = ''
        text += ' '*(2*width + 2)
        for col in xrange(self.num_columns()):
            text += _fmt % (col + self.xmin)
        text += '\n'
        text += ' '*(2*width + 2)
        for col in xrange(self.num_columns()):
            text += _fmt % col
        text += '\n'
        text += ' '*(2*width + 1) + '+' + '-'*(width*self.num_columns()) + '\n'
        for row in xrange(self.num_rows()):
            text += (_fmt % (row + self.ymin)) + ' ' + (_fmt % row) + '|'
            for col in xrange(self.num_columns()):
                adc = self.adc_value(col, row)
                pix = _fmt % adc
                if color and adc == _max:
                    pix = '%s%s%s' %\
                          (xpeAnsiColors.RED, pix, xpeAnsiColors.ENDC)
                elif color and adc >= max_threshold*_max:
                    pix = '%s%s%s' %\
                          (xpeAnsiColors.YELLOW, pix, xpeAnsiColors.ENDC)
                elif color and adc > zero_suppression:
                    pix = '%s%s%s' %\
                          (xpeAnsiColors.GREEN, pix, xpeAnsiColors.ENDC)
                text += pix
            text += '\n%s|\n' % (' '*(2*width + 1))
        return text

    def draw_ascii(self, zero_suppression=9):
        """Print the ASCII representation of the event.
        """
        print(self.ascii(zero_suppression))

    def window_matrix(self):
        """Return an xpeHexagonalMatrix object corresponding to the
        readout window.
        """
        return xpeHexagonalMatrix(self.num_columns(), self.num_rows(),
                                  self.xmin, self.ymin)

    def draw(self, zero_suppression=9, grids=True, show=True):
        """
        """
        matrix = self.window_matrix()
        matrix.draw(self.adc_values, zero_suppression, grids=grids, show=False)
        if show:
            plt.show()

    def __str__(self):
        """String representation.
        """
        text = 'buffer %5d, w(%3d, %3d)--(%3d, %3d), %d px, t = %d us' %\
               (self.buffer_id, self.xmin, self.ymin, self.xmax, self.ymax,
                self.num_pixels(), self.microseconds)
        return text






