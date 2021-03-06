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
EVT_HEADER    = 0xff 
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

    def __write2file(self, wordseq, nbytes):
        # write a word at a time, to avoid confusion with little/big endian
        # one word is unsigned short, 16 bits! forcing big endian with '>'
        for i in xrange(nbytes-1, -1, -1):
            #self.__outfile.write(struct.pack('>H', (wordseq>>(16*i))&0xffff))
            self.__outfile.write(struct.pack('>B', (wordseq>>(8*i))&0xff))

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
        self.__write2file(time_out, 6)


        # ROI 9+9 + 9+9 bits
        xmin = evt.xmin
        ymin = evt.ymin
        deltax = evt.num_columns()
        deltay = evt.num_rows()
        #assert (deltax<64 and deltay<64), "ROI length must be <2**6 "+\
        #    "but is (%d, %d)" % (deltax, deltay)
        if Verbose:
            print "ROI (dec):", xmin, ymin, deltax, deltay
            print "ROI     (hex):", hex(xmin), hex(ymin),\
                hex(deltax), hex(deltay)
            print "ROI     (bin):", bin(xmin), bin(ymin), \
                bin(deltax), bin(deltay)

        roi_out =  ((xmin&0x1ff)<<31) | \
                   ((ymin&0x1ff)<<22) | \
                   ((deltax&0x1ff)<<13) | \
                   ((deltay&0x1ff)<<4)
        if Verbose:
            print "ROI OUT (hex):", hex(roi_out)
            print "ROI OUT (bin):", bin(roi_out)
        self.__write2file(roi_out, 5)


        # HIT 2+12/18 + 4/2 type E (for PH) or type C ()
        # Marker: C 0b10 E 0b00 C-end-of-package 0b11      
        if Verbose:
            print "HIT ROI Size: %d x %d (tot %d)" % \
                (evt.num_columns(), evt.num_rows(), evt.num_pixels()) 
            print "nPixels above thr (%d)= %d" % \
                (self.__hit_threshold,
                 len(evt.adc_values[evt.adc_values > self.__hit_threshold]))
            print evt.ascii(0)

        hit_pixel_count = 0
        hit_byte_count  = 0
        ContiguousHitFlag = 0
        hit_out = 0x0
        for yId in xrange(evt.num_rows()):
            for xId in xrange(evt.num_columns()):
                currentX = xId+xmin
                currentY = yId+ymin
                currentH = evt.adc_value(xId,yId)
                if currentH>self.__hit_threshold: # above thr
                    if ContiguousHitFlag == 0: # no hit before
                        hit_out = (0b10<<22) | \
                                  ((xId&0x1ff)<<13) | ((yId&0x1ff)<<4)
                        hit_out = (hit_out<<16) | ((currentH&0xfff)<<2)
                        ContiguousHitFlag = 1
                        self.__write2file(hit_out, 5)
                        hit_byte_count  +=5
                        hit_pixel_count +=1
                        if Verbose:
                            print "HIT NEW", xId, yId, currentH 
                            print "HIT OUT", hex(hit_out), bin(hit_out)
                    else:
                        hit_out = ((currentH&0xfff)<<2)
                        self.__write2file(hit_out, 2)
                        hit_byte_count  +=2
                        hit_pixel_count +=1
                        if Verbose:
                            print "HIT CONT", xId, yId, currentH
                            print "HIT OUT", hex(hit_out), bin(hit_out)
                else:
                    ContiguousHitFlag = 0
        # End of package
        hit_out = (0b11<<22) | ((xId&0x1ff)<<13) | ((yId&0x1ff)<<4)
        self.__write2file(hit_out, 3)
        hit_byte_count +=3
        if Verbose:
            print "HIT LAST", xId, yId
            print "HIT OUT", hex(hit_out), bin(hit_out)
        
        return (evt.num_pixels(), hit_byte_count, hit_pixel_count)
                        

        

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

    list_n_hits   = [] # with no zero suppression each hit is 1 word
    list_n_words  = []
    list_n_hits_ovrThr = []
    print "Compressing events with thr %d " % args.zero_suppression
    for i in xrange(args.num_events):
        event = input_file.next()
        (nHit, nWords, nHitOverThr) = output_file.addEvent(event, False)
        list_n_hits.append(nHit)
        list_n_words.append(nWords)
        list_n_hits_ovrThr.append(nHitOverThr)
        #print "Evt %d, n hits = %d, n word in out file = %d" % (i, nHit, nWords)

    if args.num_events>1:
        import numpy
        list_n_hits  = numpy.array(list_n_hits)
        list_n_words = numpy.array(list_n_words)
        list_n_hits_ovrThr  = numpy.array(list_n_hits_ovrThr)
        len_n_hits = len(list_n_hits)
        ave_n_hits = float(sum(list_n_hits))/len_n_hits
        rms_n_hits = numpy.sqrt(sum((list_n_hits-ave_n_hits)**2)/(len_n_hits -1))
        len_n_words = len(list_n_words)
        ave_n_words = float(sum(list_n_words))/len_n_words
        rms_n_words = numpy.sqrt(sum((list_n_words-ave_n_words)**2)/(len_n_words -1))
        len_n_hits_ovrThr = len(list_n_hits_ovrThr)
        ave_n_hits_ovrThr = float(sum(list_n_hits_ovrThr))/len_n_hits_ovrThr
        rms_n_hits_ovrThr = numpy.sqrt(sum((list_n_hits_ovrThr-ave_n_hits_ovrThr)**2)/(len_n_hits_ovrThr -1))
        
        print "Number of hits stats: [min, max, average, rms] = [%d, %d, %f, %f]" %\
            (min(list_n_hits), max(list_n_hits), ave_n_hits, rms_n_hits )
        
        print "Number of bytes stats: [min, max, average, rms] = [%d, %d, %f, %f]"%\
            (min(list_n_words), max(list_n_words), ave_n_words, rms_n_words)

        print "Number of hits_ovrThr stats: [min, max, average, rms] = [%d, %d, %f, %f]" %\
            (min(list_n_hits_ovrThr), max(list_n_hits_ovrThr), ave_n_hits_ovrThr, rms_n_hits_ovrThr)
    
        n_address = (list_n_words -2.*list_n_hits_ovrThr)/3 -1
        print "Averave Number of Address = ", sum(n_address)/len(n_address)

        n_byte_per_hit = 1.*list_n_words/list_n_hits_ovrThr
        print "Averave Number of byte/hit = ", sum(n_byte_per_hit)/len(n_byte_per_hit)

        print "SUMMARY: thr, <nHit>, <nByte>, sigma_nByte, <byte/nHit>", \
            args.zero_suppression, ave_n_hits_ovrThr, ave_n_words, rms_n_words, \
            sum(n_byte_per_hit)/len(n_byte_per_hit)

        import matplotlib.pyplot as plt
        plt.grid(color='gray')
        plt.title("Number of bytes in zero-suppressed output file")
        plt.ylabel("Entries/bin")
        plt.xlabel("# Bytes")
        plt.hist(list_n_words, bins=200, range=(0,1500))
        plt.show()
