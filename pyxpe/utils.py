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


import time


class xpeChrono:

    """Small chronometer class.

    A chronometer essentially measures the elapsed time since it has been
    started and is equipped to print itself to the standard output. (Note the
    chronometer is reset unpon the instantiation of a class object.)

    Examples
    --------
    >>> from ximpol.utils.profile import xChrono
    >>> c = xChrono()
    >>> # ... do something.
    >>> print(c)
    """

    def __init__(self):
        """Constructor.
        """
        self.reset()

    def reset(self):
        """Reset the chronometer.
        """
        self.start_time = time.time()

    def __call__(self):
        """Return the elapsed time.
        """
        return time.time() - self.start_time

    def __str__(self):
        """ String formatting.
        """
        return '[t0 + %.3f s]' % self()


def main():
    c = xpeChrono()
    print(c)
    time.sleep(2)
    print(c)
    c.reset()
    print(c)


if __name__ == '__main__':
    main()
