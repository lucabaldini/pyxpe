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


# Wrapper lib for random number generation
# to change engine (python, root, heprep), by changing only this file

# default now is python pseudo-random generator
# https://docs.python.org/2.7/library/random.html
import random

class xperandom:
    """Wrapper random class
    """

    def __init__(self):
        self.engine = random

    def get_engine(self):
        """ Return the engine, in order to use it directly
        """
        return self.engine
    
    def set_seed(self, seed=None):
        """ Set the seed of the pseudo-random generator
        """
        self.engine.seed(seed)

    def random(self):
        """ Main engine call
        Return the next random floating point number in the range [0.0, 1.0).
        """
        return self.engine.random()

    def get_state(self):
        """
        """
        return self.engine.getstate()

    def set_state(self, state):
        """
        """
        return self.engine.setstate(state)
    

        
if __name__ == '__main__':
     r = xperandom()
     
     r.get_engine().seed(2)
     print("Direct engine call:", \
         r.get_engine().random(), \
         r.get_engine().random())
     
     r.set_seed(2)
     print("Wrapped call:", r.random(), r.random())

     currentstate = r.get_state()
     #print ("Get state", currentstate)
     print("NextRandom", r.random())
     r.set_state(currentstate)
     print("PrevRandom", r.random())
     
