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


__description__ = 'Run the XPE reconstruction'


import os
from pyxpe.utils.logging_ import startmsg
from pyxpe.recon.xpol import XPOL_COORDINATE_SYSTEMS
from pyxpe.recon.recon import run_pixy


"""Command-line switches.
"""
import argparse
import ast

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('binfile', type=str,
                    help='path to the input xpedaq binary file')
PARSER.add_argument('-o', '--output-path', type=str, default=None,
                    help='path to the output ROOT file')
PARSER.add_argument('-n', '--num-events', type=int, default=1000000000,
                    help = 'the maximum number of events to be processed')
PARSER.add_argument('-z', '--zero-suppression', type=int, default=9,
                    help = 'the zero-suppression threshold')
PARSER.add_argument('-c', '--coordinate-system', type=str, default='pixy',
                    choices=XPOL_COORDINATE_SYSTEMS,
                    help = 'the coordinate system for the clustering')
#PARSER.add_argument('--clobber', type=ast.literal_eval, choices=[True, False],
#                    default=True,
#                    help='overwrite or do not overwrite existing output files')



def xpepixy(file_path, **kwargs):
    """
    """    
    return run_pixy(file_path, **kwargs)



if __name__=='__main__':
    args = PARSER.parse_args()
    startmsg()
    file_path = args.__dict__.pop('binfile')
    xpepixy(file_path, **args.__dict__)
