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


import ROOT
import numpy

from pyxpe.utils.logging_ import logger


ROOT_TO_NUMPY_DICT = {'C': '',
                      'B' : 'int8',
                      'b' : 'uint8',
                      'S' : 'int16',
                      's' : 'uint16',
                      'I' : 'int32',
                      'i' : 'uint32',
                      'F' : 'float32',
                      'D' : 'float64',
                      'L' : 'int64',
                      'l' : 'uint64',
                      'O' : 'bool_'
                  }


def root2numpy(rootType):
    """ Convert a ROOT branch type into the corresponding numpy type.
    """
    return ROOT_TO_NUMPY_DICT[rootType]




class xpeRootTreeBase(ROOT.TTree):
    
    """ Small wrapper around the ROOT.Tree class.

    Compared to the base class, this class provides some additional facility to
    create and manipulate branches and it's mainly used to write ROOT tree
    to file.

    The list of supported types, straight from the ROOT documentation, is:
    - C : a character string terminated by the 0 character
    - B : an 8 bit signed integer (Char_t)
    - b : an 8 bit unsigned integer (UChar_t)
    - S : a 16 bit signed integer (Short_t)
    - s : a 16 bit unsigned integer (UShort_t)
    - I : a 32 bit signed integer (Int_t)
    - i : a 32 bit unsigned integer (UInt_t)
    - F : a 32 bit floating point (Float_t)
    - D : a 64 bit floating point (Double_t)
    - L : a 64 bit signed integer (Long64_t)
    - l : a 64 bit unsigned integer (ULong64_t)
    - O : [the letter 'o', not a zero] a boolean (Bool_t)
    """

    NAME = 'Tree'
    TITLE = 'Tree'
    BRANCHES = []
    ALIAS_DICT = {}

    def __init__(self):
        """ Constructor.
        """
        logger.info('Creating ROOT tree %s...' % self.NAME)
        ROOT.TTree.__init__(self, self.NAME, self.TITLE)
        self.__array_dict = {}
        for branch_name, branch_type in self.BRANCHES:
            self.add_branch(branch_name, branch_type)
        self.branch_list = [item[0] for item in self.BRANCHES]
        for key, value in self.ALIAS_DICT.items():
            logger.info('Setting alias "%s" -> "%s"...' % (key, value))
            self.SetAlias(key, value)

    def add_branch(self, branch_name, branch_type):
        """ Add a branch to the output tree.
        """
        branchTitle = '%s/%s' % (branch_name, branch_type)
        logger.info('Adding branch %s to %s...' % (branchTitle, self.GetName()))
        a = numpy.array([0], dtype = root2numpy(branch_type))
        self.__array_dict[branch_name] = a
        self.Branch(branch_name, a, branchTitle)

    def fill(self, entry_dict):
        """ Fill a row of the tree.

        This really set the value for all the branches of the tree and calls the
        ROOT.TTree.Fill() method at the end.

        The argument row is essentially a dictionary which is supposed to
        contain all the branch names as keys (note that the loop is actually
        done over the branch names, so that the dictionary can contain a
        superset of the branch names as its keys.)
        """
        for branch_name in self.branch_list:
            self.set_value(branch_name, entry_dict[branch_name])
        self.Fill()

    def set_value(self, branch_name, value):
        """ Set the value of a specific array.
        """
        self.__array_dict[branch_name][0] = value

    def value(self, branch_name, entry=None):
        """ Return the value of a specific array.
        """
        if entry is not None:
            self.GetEntry(entry)
        return self.__array_dict[branch_name][0]



class xpePixyTree(xpeRootTreeBase):

    NAME = 'tree'
    TITLE = 'Analized Data Tree'
    BRANCHES = [
        ('fRunId', 'i'),
        ('fEventId', 'I'),
        ('fNClusters', 'I'),
        ('fTrigWindow', 'I'),
        ('fTimeTick', 'l'),
        ('fTimeStamp', 'D'),
        ('fBufferId', 'I'),
        ('fCluSize', 'I'),
        ('fPHeight', 'F'),
        ('fStoN', 'F'),
        ('fTotNoise', 'F'),
        ('fBaricenterX', 'F'),
        ('fBaricenterY', 'F'),
        ('fTheta0', 'F'),
        ('fTheta1', 'F'),
        ('fMomX', 'F'),
        ('fMomY', 'F'),
        ('fMomThirdX', 'F'),
        ('fImpactX', 'F'),
        ('fImpactY', 'F')
    ]



def test():
    """ Test code.
    """
    t = xpeRootTree()
    t.add_branch('Var1', 'I')
    t.add_branch('Var2', 'F')
    for i in range(10):
        t.set_value('Var1', i)
        t.set_value('Var2', i**0.5)
        t.Fill()
    t.Draw('Var2')
    raw_input('Press enter to exit')



if __name__ == '__main__':
    test()
