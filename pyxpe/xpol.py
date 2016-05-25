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


import numpy

import matplotlib.pyplot as plt
from matplotlib import collections, transforms


class pHexagonalMatrix():

    """
    """

    COLUMN_PITCH = 0.0500
    ROW_PITCH = 0.0433

    def __init__(self, num_columns, num_rows, start_column=0, start_row=0):
        """
        """
        self.num_columns = num_columns
        self.num_rows = num_rows
        self.start_column = start_column
        self.start_row = start_row
        self.__grid = []
        for col in xrange(start_column, start_column + num_columns):
            for row in xrange(start_row, start_row + num_rows):
                self.__grid.append(self.pixel2world(col, row))

    def grid(self):
        """
        """
        return self.__grid

    def pixel2world(self, col, row):
        """Convert from pixel coordinates to world coordinates.
        """
        return self.pixel2world_asic(col, row)

    def pixel2world_asic(self, col, row):
        """Convert from pixel coordinates to world coordinates.

        This is using the DAQ coordinate system, where the origin is in the
        top-left corner of the array, the x coordinate is spanning the columns
        from left to right and the y coordinate is spanning the rows from top
        to bottom.
        """
        x = (col + 0.5*(row % 2))*self.COLUMN_PITCH
        y = -row*self.ROW_PITCH
        return (x, y)

    def pixel2world_recon(self, col, row):
        """Convert from pixel coordinates to world coordinates.

        This is using the reconstruction coordinate system.
        """
        x = (row - 0.5*(self.num_rows - 1))*self.ROW_PITCH
        y = (col - 0.5*(self.num_columns - 0.5 + row % 2))*self.COLUMN_PITCH
        return (x, y)

    def asic2recon(self, x, y):
        """Convert from ASIC coordinates to recon coordinates.
        """
        _x = -(y + 0.5*(self.num_rows - 1)*self.ROW_PITCH)
        _y = x - 0.5*(self.num_columns - 0.5)*self.COLUMN_PITCH
        return (_x, _y)

    def recon2asic(self, x, y):
        """Convert from recon coordiates to ASIC coordinates.
        """
        _x = y + 0.5*(self.num_columns - 0.5)*self.COLUMN_PITCH
        _y = -(x + 0.5*(self.num_rows - 1)*self.ROW_PITCH)
        return (_x, _y)

    def border(self, col, row):
        """Return true if the pixel at the specified position is on the
        border of the array.
        """
        return (col == 0) or (col == self.num_columns - 1) or\
            (row == 0) or (row == self.num_rows - 1)

    def write_pixmap(self, filePath):
        """Write a pixmap.dat-like file containing the pixel hash table used
        by the reconstruction.

        352
        300
        x y u v ch mask border
        -7.59915 -7.4875 0 0 0 1 1
        -7.59915 -7.4375 0 1 1 1 1
        ...
        7.59915 7.3875 351 298 105598 1 1
        7.59915 7.4375 351 299 105599 1 1
        """
        chan = 0
        mask = 1
        f = open(filePath, 'w')
        f.write('%d\n%d\n' % (self.num_rows, self.num_columns))
        f.write('x y u v ch mask border\n')
        for row in xrange(self.start_row, self.start_row + self.num_rows):
            for col in xrange(self.start_column, self.start_column +
                              self.num_columns):
                x, y = self.pixel2world_recon(col, row)
                bord = self.border(col, row)
                line = '%.5f %.4f %d %d %d %d %d\n' %\
                       (x, y, row, col, chan, mask, bord)
                f.write(line)
                chan += 1
        f.close()

    def world2pixel_asic(self, x, y):
        """
        """
        row = int(-y/self.ROW_PITCH + 0.5)
        dx = x % self.COLUMN_PITCH
        row += int(abs(dx - 0.5*self.COLUMN_PITCH)/(0.5*self.COLUMN_PITCH))
        col = int(x/self.COLUMN_PITCH + 0.5) - (row % 2)
        return col, row

    def frame(self, padding=0.1):
        """Return a (xmin, ymin, xmax, ymax) containing the entire matrix.
        """
        xmin, ymin = self.pixel2world(self.start_column,
                                      self.start_row + self.num_rows - 1)
        xmax, ymax = self.pixel2world(self.start_column + self.num_columns - 1,
                                      self.start_row)
        # Make sure we include the pixel edges.
        xmin -= 0.5*self.COLUMN_PITCH
        xmax += 0.5*self.COLUMN_PITCH
        ymin -= 0.5*self.ROW_PITCH
        ymin += 0.5*self.ROW_PITCH
        # Add the padding.
        dx = padding*(xmax - xmin)
        dy = padding*(ymax - ymin)
        xmin -= dx
        xmax += dx
        ymin -= dy
        ymax += dy
        return xmin, ymin, xmax, ymax

    def draw(self, colors='white', show=True):
        """
        """
        xmin, ymin, xmax, ymax = self.frame()
        fig = plt.figure(figsize=(10, 10), dpi=80, facecolor='w')
        ax = plt.subplot(111, aspect='equal', adjustable='box-forced')
        ax.set_xlim([xmin, xmax])
        ax.set_ylim([ymin, ymax])
        # Overall average canvas dimensions in pixels.
        dim = (ax.transData.transform((xmin, ymin)) -
               ax.transData.transform((xmax, ymax))).mean()
        # Calculate a something proportional to the hexagon area in px**2.
        scale = max((ymax - ymin), (xmax - xmin))
        dim = (dim/scale*0.8*self.COLUMN_PITCH)**2
        # Create the hexagin collection.
        poly = collections.RegularPolyCollection(6,
                                                 offsets=self.__grid,
                                                 sizes=(dim,),
                                                 transOffset=ax.transData,
                                                 facecolors=colors,
                                                 edgecolors='gray')
        ax.add_collection(poly, autolim=True)
        plt.grid()

        def _htxt(x, y, s, **kwargs):
            """
            """
            return plt.text(x, y + 0.05, '%s' % s, horizontalalignment='center',
                            verticalalignment='center', **kwargs)

        def _vtxt(x, y, s, **kwargs):
            """
            """
            return plt.text(x - 0.03, y, '%s' % s, horizontalalignment='right',
                            verticalalignment='center', **kwargs)  

        x1, y1 = self.__grid[0]
        x2, y2 = self.__grid[-2]
        x3, y3 = self.__grid[-1]
        ht1 = _htxt(x1, y1, self.start_column)
        ht2 = _htxt(x2, y1, self.start_column + self.num_columns - 1)
        ht3 = _htxt(0.5*(x1 + x2), y1, '--- Readout column ---')        
        vt1 = _vtxt(x1, y1, self.start_row)
        vt2 = _vtxt(x1, y3, self.start_row + self.num_rows - 1)
        vt3 = _vtxt(x1, 0.5*(y1 + y3), '--- Readout row ---', rotation=90.)
        plt.xlabel('ASIC reference frame x [mm]')
        plt.ylabel('ASIC reference frame y [mm]')
        if show:
            plt.show()
        return fig



if __name__ == '__main__':
    matrix = pHexagonalMatrix(10, 12, 0, 0)
    #matrix.write_pixmap('pixmap.dat')
    import random
    x, y = random.uniform(0, 0.45), random.uniform(-0.45, 0)
    col, row = matrix.world2pixel_asic(x, y)
    print col, row
    matrix.draw(show=False)
    plt.plot([x], [y], 'o')
    plt.show()

    
