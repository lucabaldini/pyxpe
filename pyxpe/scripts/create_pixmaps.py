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


import os

from pyxpe.recon.xpol import xpeXpolMatrix, XPOL_COORDINATE_SYSTEMS
from pyxpe import PYXPE_DATA


OUTPUT_FOLDER = os.path.join(PYXPE_DATA, 'pixmaps')

matrix = xpeXpolMatrix()
for coordinate_system in XPOL_COORDINATE_SYSTEMS:
    filePath = os.path.join(OUTPUT_FOLDER, 'pixmap_%s.dat' % coordinate_system)
    matrix.write_pixmap_ascii(filePath, coordinate_system)
    filePath = os.path.join(OUTPUT_FOLDER, 'pixmap_%s.fits' % coordinate_system)
    matrix.write_pixmap_fits(filePath, coordinate_system)
