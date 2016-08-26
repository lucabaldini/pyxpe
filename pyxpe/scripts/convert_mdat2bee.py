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

"""
This script convert a .mdat file in the binary format expected from the BEE
"""

import argparse
from pyxpe.recon.binio import xpeBinaryFileWindowed
from pyxpe.recon.event import xpeEventWindowed
import struct


# Settings for output file
EVT_HEADER    = 0xffff # None to skip header
HIT_THRESHOLD = 3 # for zero suppression
TIME_OFFSET   = 400*24*3600 # arbitrary offset for event time, 400 days...

class xpeBinaryBEEOutput:
    def __init__(self, outfilePath, hit_threshold, header_flag, time_offset):
        """Constructor.
        """
        self.__outfilePath   = "tmp_test.bee" #outfilePath
        self.__hit_threshold = hit_threshold
        self.__header_flag   = header_flag
        self.__time_offset   = time_offset
        self.__outfile = file(self.__outfilePath, "wb+") # cleanup

    def __write2file(self, wordseq, nwords):
        # write a word at a time to avoid confusion with little/big endian
        # one word is unsigned short, 16 bits! forcing big endian with '>'
        for i in xrange(nwords-1, -1, -1):
            self.__outfile.write(struct.pack('>H', (wordseq>>(16*i))&0xffff))

    def addEvent(self, evt, Verbose=True):
        if Verbose:
            print "\n", evt
        
        if self.__header_flag:
            self.__write2file(EVT_HEADER, 1)
        # TIME 28 + 20 bits
        # evt.microseconds isn time from FPGA, ignoring start runs bytes!
        time_s  = int(evt.microseconds/1e6)+ self.__time_offset
        time_us = int(evt.microseconds%1e6)
        if Verbose:
            print "TIME FPGA:", evt.microseconds, "(without offset)"
            print "TIME sec:", time_s, " us: ", time_us
            print "TIME     (hex) sec:", hex(time_s), " us: ", hex(time_us)
            print "TIME     (bin) sec:", bin(time_s), " us: ", bin(time_us)
        time_out = time_s<<20 | time_us # THIS IS THE TRICK
        if Verbose:
            print "TIME OUT (hex) sec:", hex(time_out)
            print "TIME OUT (bin) sec:", bin(time_out)
        self.__write2file(time_out, 3)

        # ROI 9+9 + 6+6 bits
        xmin = evt.xmin
        ymin = evt.ymin
        deltax = evt.xmax - evt.xmin
        deltay = evt.ymax - evt.ymin
        assert (deltax<64 and deltay<64), "ROI length must be <2**6 "+\
            "but is (%d, %d)" % (deltax, deltay)
        if Verbose:
            print "ROI (dec):", xmin, ymin, deltax, deltay
            print "ROI     (hex):", hex(xmin), hex(ymin),\
                hex(deltax), hex(deltay)
            print "ROI     (bin):", bin(xmin), bin(ymin), \
                bin(deltax), bin(deltay)

        roi_out =  ((xmin&0x1ff)<<23) | \
                   ((ymin&0x1ff)<<14) | \
                   ((deltax&0x3f)<<8) | \
                   ((deltay&0x3f)<<2)
        if Verbose:
            print "ROI OUT (hex):", hex(roi_out)
            print "ROI OUT (bin):", bin(roi_out)
        self.__write2file(roi_out, 2)

        # HIT
        # scan...

        
        
if __name__ == '__main__':
    formatter = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=formatter)
    parser.add_argument('binfile', type=str,
                        help='the input binary file')
    parser.add_argument('-n', '--num-events', type=int, default=10,
                        help='number of events to be read')
    parser.add_argument('-z', '--zero_suppression', type=int,
                        default=HIT_THRESHOLD,
                        help='zero-suppression threshold')
    parser.add_argument('-t', '--time_offset', type=int, default=TIME_OFFSET,
                        help='set event time offset')
    parser.add_argument('-H', '--header', action="store_true", default=False,
                        help='add event header')
    args = parser.parse_args()

    
    input_file = xpeBinaryFileWindowed(args.binfile)
    out_file_name = args.binfile
    assert out_file_name.endswith('.mdat'), "Input file MUST be a .mdat file"
    out_file_name = out_file_name.replace('.mdat', '.bee')
    output_file = xpeBinaryBEEOutput(out_file_name, args.zero_suppression,
                                     args.header, args.time_offset)
    for i in xrange(args.num_events):
        event = input_file.next()
        output_file.addEvent(event)



    
