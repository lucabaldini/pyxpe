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

from pyxpe.recon.event import xpeEventWindowed, xpeEventFullFrame
from pyxpe.recon.xpol import XPOL_NUM_PIXELS, XPOL_NUM_COLUMNS, XPOL_NUM_ROWS
from pyxpe.recon.xpol import XPOL_PIXELS_PER_BUFFER, XPOL_NUM_BUFFERS



class xpeBinaryFileBase:

    """ Base class for a xpedaq binary file.
    """

    def __init__(self, filePath):
        """Constructor.
        """
        logger.info('Opening input binary file %s...' % filePath)
        self.__file = open(filePath, 'rb')

    def seek(self, offset):
        """ redefine seek
        """
        self.__file.seek(offset)

    def read(self, n):
        """ redefine read
        """
        return self.__file.read(n)

    def close(self):
        """ redefine
        """
        self.__file.close()

    def read_word(self):
        """Read and byte-swap a single 2-bytes binary word from file.

        Note that struct.unpack returns a tuple even when we read a single
        number, and here we're returning the first (and only) element of the
        tuple.
        """
        return struct.unpack('H', self.read(2))[0]

    def read_words(self, num_words):
        """Read and byte-swap a fixed number of 2-bytes binary words from file.
        
        Args
        ----
        num_words : int
            The number of words to be read from the input file.
        """
        return struct.unpack('%dH' % num_words, self.read(2*num_words))

    def __iter__(self):
        """Basic iterator implementation.
        """
        return self

    def next(self):
        """Do-nothing next() method to be reimplemented in the derived classes.
        """
        pass

        

class xpeBinaryFileFullFrame(xpeBinaryFileBase):
    
    """Binary file acquired in full-frame mode.
    """

    def next(self):
        """Read the next event in the file.

        Warning
        -------
        We should return an event object instead of a plain numpy array.
        """
        try:
            data = self.read_words(XPOL_NUM_PIXELS)
        except Exception:
            raise StopIteration()
        adc_counts = numpy.array(data, numpy.uint16)
        adc_counts = adc_counts.reshape(XPOL_PIXELS_PER_BUFFER,
                                        XPOL_NUM_BUFFERS)
        adc_counts = adc_counts.transpose()     
        adc_counts = adc_counts.flatten()
        adc_counts = adc_counts.reshape(XPOL_NUM_ROWS, XPOL_NUM_COLUMNS).T
        return xpeEventFullFrame(adc_counts)           


class xpeBinaryFileWindowed(xpeBinaryFileBase):
    
    """Binary file acquired in windowed mode.
    """

    def next(self):
        """Read the next event in the file.
        """
        try:
            header = self.read_word()
        except Exception:
            raise StopIteration()
        if header != xpeEventWindowed.HEADER_MARKER:
            msg = 'Event header mismatch at byte %d' % self.tell()
            msg += ' (expected %s, got %s).' %\
                   (hex(xpeEventWindowed.HEADER_MARKER), hex(header))
            logger.error(msg)
            logger.info('Moving ahead to the next event header...')
            while header != xpeEventWindowed.HEADER_MARKER:
                header = self.read_word()
            logger.info('Got back in synch at byte %d.' % self.tell())
        xmin, xmax, ymin, ymax, buf_id, t1, t2, s1, s2 = self.read_words(9)
        num_columns = (xmax - xmin + 1)
        num_rows = (ymax - ymin + 1)
        data = self.read_words(num_rows*num_columns)
        adc = numpy.array(data).reshape((num_rows, num_columns)).T
        return xpeEventWindowed(xmin, xmax, ymin, ymax, buf_id, t1, t2, s1, s2,
                                adc)
    

def open_binary_file(filePath):
    """Open an xpedaq binary file.

    This is meant to be a generic interface to all kind of files produce by the
    DAQ.
    """
    word = xpeBinaryFileBase(filePath).read_word()
    if word == xpeEventWindowed.HEADER_MARKER:
        # No-header file in windowed mode.
        return xpeBinaryFileWindowed(filePath)
    elif word < 0x0fff:
        # No-header file in full-frame mode.
        return xpeBinaryFileFullFrame(filePath)

    
def test_fullframe(filePath, num_events):
    """ scicazzi
    """
    input_file = xpeBinaryFileFullFrame(filePath)
    for i in xrange(args.num_events):
        event = input_file.next()
        print (event)
    return event



def test_windowed(filePath, num_events):
    """Read a windowed event file and display the events.
    """
    input_file = xpeBinaryFileWindowed(filePath)
    for i in xrange(args.num_events):
        event = input_file.next()
        print (event)
        event.draw_ascii()
        event.draw()

        
if __name__ == '__main__':
    import argparse
    formatter = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=formatter)
    parser.add_argument('infile', type=str,
                        help='the input binary file')
    parser.add_argument('-n', '--num_events', type=int, default=10,
                        help = 'number of events to be processed')
    args = parser.parse_args()
    #test_windowed(args.infile, args.num_events)
    evt = test_fullframe(args.infile, args.num_events)
