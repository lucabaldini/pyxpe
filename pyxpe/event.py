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


import logging
logging.basicConfig(level=logging.DEBUG)

import struct
import numpy

import matplotlib
import matplotlib.pyplot as plt

import pyxpedaq.xpol as xpol



PIXELS_PER_BUFFER = 13200
NUM_BUFFERS = 8
NUM_PIXELS = PIXELS_PER_BUFFER*NUM_BUFFERS


class pAnsiColors:
    
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


    
class pXpeEventBase:

    """
    """

    pass



class pXpeEventFullFrame(pXpeEventBase):

    """
    """

    pass



class pXpeEventWindowed(pXpeEventBase):
    
    """
    """
    
    HEADER_MARKER = 65535
    HEADER_LENGTH = 20
    
    def __init__(self, xmin, xmax, ymin, ymax, buffer_id, t1, t2, s1, s2,
                 adc_counts):
        """Constructor.
        """
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.buffer_id = buffer_id
        self.t1 = t1
        self.t2 = t2
        self.s1 = s1
        self.s2 = s2
        self.adc_counts = adc_counts

    def size(self):
        """Return the total number of bytes in the event.
        """
        return self.HEADER_LENGTH + 2*self.num_pixels()

    def num_columns(self):
        """
        """
        return (self.xmax - self.xmin + 1)

    def num_rows(self):
        """
        """
        return (self.ymax - self.ymin + 1)

    def num_pixels(self):
        """Return the total number of pixels in the window.
        """
        return self.num_rows()*self.num_columns()

    def microseconds(self):
        """Return the timestamp from the FPGA.
        """
        return (self.t1 + self.t2*65534)*0.8e-6

    def start_seconds(self):
        """Return the seconds from the start run.
        """
        return (self.s2 + self.s1*65536)
        
    def adc_value(self, col, row):
        """Return the pulse height for a given pixel in the window.
        """
        return self.adc_counts[col, row]

    def highest_pixel(self):
        """Return the coordinats of the pixel with the maximum value of
        ADC counts.
        """
        return numpy.unravel_index(numpy.argmax(self.adc_counts),
                                   self.adc_counts.shape)

    def highest_adc_value(self):
        """Return the maximum value of ADC counts for the pixels in the event.
        """
        return self.adc_counts.max()

    def pulse_height(self, zero_suppression=10):
        """Return the total pulse height for the event, i.e., the raw sum of
        all the ADC values above the zero-suppression threshold.
        """
        return self.adc_counts[self.adc_counts >= zero_suppression].sum()

    def baricenter(self, zero_suppression=10):
        """Return the baricenter of the event.
        """
        return None

    def ascii(self, zero_suppression=10, max_threshold=0.75, width=4,
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
                    pix = '%s%s%s' % (pAnsiColors.RED, pix, pAnsiColors.ENDC)
                elif color and adc >= max_threshold*_max:
                    pix = '%s%s%s' % (pAnsiColors.YELLOW, pix, pAnsiColors.ENDC)
                elif color and adc >= zero_suppression:
                    pix = '%s%s%s' % (pAnsiColors.GREEN, pix, pAnsiColors.ENDC)
                text += pix
            text += '\n%s|\n' % (' '*(2*width + 1))
        return text

    def draw_ascii(self, threshold=10):
        """
        """
        print(self.ascii())

    def draw(self, zero_suppression=10):
        """
        """
        matrix = xpol.pHexagonalMatrix(self.num_columns(), self.num_rows(),
                                       self.xmin, self.ymin)
        _vals = self.adc_counts.flatten()
        _maxval = float(self.adc_counts.max())
        _vals[_vals < zero_suppression] = -1.
        _vals = _vals/_maxval
        cmap = matplotlib.cm.get_cmap('Reds')
        cmap.set_under('white')
        _colors = cmap(_vals)
        matrix.draw(colors=_colors, show=False)
        for (x, y), val in zip(matrix.grid(), self.adc_counts.flatten()):
            if val >= zero_suppression:
                if val < 0.5*_maxval:
                    col = 'black'
                else:
                    col = 'white'
                plt.text(x, y, '%s' % val, horizontalalignment='center',
                         verticalalignment='center', size=8, color=col)
        plt.show()

    def __str__(self):
        """
        """
        text = 'buffer %5d, w(%3d, %3d)--(%3d, %3d), %d px, t = %d + %.5f s' %\
               (self.buffer_id, self.xmin, self.ymin, self.xmax, self.ymax,
                self.num_pixels(), self.start_seconds(), self.microseconds())
        return text



class pXpeBinaryFileBase(file):

    """ Base class for a xpedaq binary file.
    """

    def __init__(self, filePath):
        """Constructor.
        """
        logging.info('Opening input binary file %s...' % filePath)
        file.__init__(self, filePath, 'rb')

    def read_word(self):
        """Read and byte-swap a single 2-bytes binary word from file.

        Note that struct.unpack returns a tuple even when we read a single
        number, and here we're returning the first (and only) element of the
        tuple.
        """
        return struct.unpack('H', self.read(2))[0]

    def read_words(self, num_words):
        """Read and byte-swap a fixed number of 2-bytes binary words from file.
        """
        return struct.unpack('%dH' % num_words, self.read(2*num_words))

    def __iter__(self):
        """Iterator implementation.
        """
        return self

    def next(self):
        """Do-nothing next() method to be reimplemented in the derived classes.
        """
        pass

        

class pXpeBinaryFileFullFrame(pXpeBinaryFileBase):
    
    """Binary file acquired in full-frame mode.
    """

    def next(self):
        """Read the next event in the file.
        """
        data = self.read_words(NUM_PIXELS)
        adc_counts = numpy.array(data).reshape(NUM_BUFFERS, PIXELS_PER_BUFFER)
        return adc_counts
            


class pXpeBinaryFileWindowed(pXpeBinaryFileBase):
    
    """Binary file acquired in windowed mode.
    """

    def next(self):
        """Read the next event in the file.
        """
        try:
            header = self.read_word()
        except Exception, e:
            print(e)
            raise StopIteration()
        if header != pXpeEventWindowed.HEADER_MARKER:
            logging.error('Event header mismatch (got %s).' % hex(header))
            raise StopIteration()
        xmin, xmax, ymin, ymax, buffer_id, t1, t2, s1, s2 = self.read_words(9)
        num_columns = (xmax - xmin + 1)
        num_rows = (ymax - ymin + 1)
        data = self.read_words(num_rows*num_columns)
        adc_counts = numpy.array(data).reshape((num_rows, num_columns)).T
        return pXpeEventWindowed(xmin, xmax, ymin, ymax,
                                 buffer_id, t1, t2, s1, s2, adc_counts)



def test_windowed(filePath, num_events):
    """
    """
    input_file = pXpeBinaryFileWindowed(filePath)
    for i in xrange(args.num_events):
        event = input_file.next()
        print event
        event.draw_ascii()
        event.draw()

        
if __name__ == '__main__':
    import argparse
    formatter = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=formatter)
    parser.add_argument('infile', type=str,
                        help='the input binary file')
    parser.add_argument('num_events', type=int,
                        help = 'number of events to be processed')                    
    args = parser.parse_args()
    test_windowed(args.infile, args.num_events)
