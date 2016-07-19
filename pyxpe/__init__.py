#!/usr/bin/env python
#
# Copyright (C) 2015, the X-ray Polarimetry Explorer (XPE) team.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU GengReral Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


import os

PACKAGE_NAME = 'pyxpe'

"""Basic folder structure of the package.
"""
PYXPE_ROOT = os.path.abspath(os.path.dirname(__file__))
PYXPE_BASE = os.path.join(PYXPE_ROOT, os.pardir)
PYXPE_BIN = os.path.join(PYXPE_BASE, 'bin')
PYXPE_RECON = os.path.join(PYXPE_ROOT, 'recon')
PYXPE_SCRIPTS = os.path.join(PYXPE_ROOT, 'scripts')
PYXPE_SIMULATION = os.path.join(PYXPE_ROOT, 'simulation')
PYXPE_UTILS = os.path.join(PYXPE_ROOT, 'utils')
