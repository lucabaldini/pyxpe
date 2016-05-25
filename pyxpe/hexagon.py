# lib for Hex grid coordinate systems
# from http://www.redblobgames.com/grids/hexagons/#pixel-to-hex

from math import sqrt

class Cube:
    def __init__(self, x,y,z):
        # Cube coordinate object, 3 coord.
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        return "Cube(%d, %d, %d)" % (self.x, self.y, self.z)

class Hex:
    def __init__(self, q, r):
        # Axial coordinate object, 2 coord.
        self.q = q
        self.r = r

    def __str__(self):
        return "Hex(%d, %d)" % (self.q, self.r)
        
class Offset:
    def __init__(self, col, row):
        # Axial coordinate object, 2 coord.
        self.row = row
        self.col = col

    def __str__(self):
        return "Offset(%d, %d)" % (self.col, self.row)

    def getNeighbors(self):
        print "not implemented yet"
        return None
        
        
# converting 'cube' (Cube) <--> 'axial' (Hex)

def cube_to_hex(h): # axial
    q = h.x
    r = h.z
    return Hex(q, r)

def hex_to_cube(h): 
    x = h.q
    z = h.r
    y = -x-z
    return Cube(x, y, z)


# converting 'offset' ("even-r" horizontal layout) <--> 'cube'

# convert cube to even-r offset
def cube_to_offset(h):
    col = h.x + (h.z + (h.z&1)) / 2
    row = h.z
    return Offset(col,row)

# convert even-r offset to cube
def offset_to_cube(h):
    x = h.col - (h.row + (h.row&1)) / 2
    z = h.row
    y = -x-z
    return Cube(x,y,z)

# converting 'offset' ("even-r" horizontal layout) <--> 'axial'
    
# convert axial to even-r offset
def hex_to_offset(h):
    col = h.q + (h.r + (h.r&1)) / 2
    row = h.r
    return Offset(col,row)

# convert even-r offset to hex
def offset_to_hex(h):
    q = h.col - (h.row + (h.row&1)) / 2
    r = h.row
    return Hex(q, r)
    
# Take a location and convert it into a hex grid coordinate
# NOTE:
# Y location goes down (positive as 'row' increases)
# with (0,0) at the center of Pixel 0,0
# Need to set dimention, all in mm
Pitch = 0.050 #mm
Size  = Pitch/sqrt(3)

def location_to_hex(x, y):
    q = (x * sqrt(3)/3 - y / 3.)/Size
    r = y * 2./3 / Size
    return hex_round(Hex(q, r))

def cube_round(h):
    rx = round(h.x)
    ry = round(h.y)
    rz = round(h.z)

    x_diff = abs(rx - h.x)
    y_diff = abs(ry - h.y)
    z_diff = abs(rz - h.z)

    if x_diff > y_diff and x_diff > z_diff:
        rx = -ry-rz
    elif y_diff > z_diff:
        ry = -rx-rz
    else:
        rz = -rx-ry
    if int(rx) != round(rx) or int(ry) != round(ry) or int(rz) != round(rz):
        print "Error, this should never happen:", rx, ry, rz
    return Cube(int(rx), int(ry), int(rz))

def hex_round(h):
    return cube_to_hex(cube_round(hex_to_cube(h)))

if __name__ == '__main__':
    print "test offset - axial"
    for (col, row) in [(0,0), (0,1), (0,2), (0,3),\
                       (1,0), (2,0), \
                       (1,1), (2,2), \
                       (299,0), (0, 351), (299,351)]:
        o = Offset(col, row)
        h  = offset_to_hex(o)
        o1 = hex_to_offset(h)
        print o, h, o1

    print "\n test location to hex"
    for (X,Y) in [(0.,0.), (0.01, 0), (0.03, 0)]:
        h = location_to_hex(X, Y)
        o = hex_to_offset(h)
        print "X,Y", X,Y, "->", h,o

    import numpy
    import matplotlib.pyplot as plt
    myX_00 = []; myY_00 = [];
    myX_01 = []; myY_01 = [];
    myX_10 = []; myY_10 = [];
    myX_11 = []; myY_11 = [];

    for iX in range(400):
        X = iX*0.0004 -0.06
        for iY in range(400):
            Y = iY*0.0004 -0.06
            h = location_to_hex(X, Y)
            o = hex_to_offset(h)
            if o.col == 0 and o.row == 0:
                myX_00.append(X)
                myY_00.append(-Y)
            if o.col == 1 and o.row == 0:
                myX_10.append(X)
                myY_10.append(-Y)
            if o.col == 0 and o.row == 1:
                myX_01.append(X)
                myY_01.append(-Y)
            if o.col == 1 and o.row == 1:
                myX_11.append(X)
                myY_11.append(-Y)
                #print X,Y

    plt.ion()
    plt.plot(numpy.array(myX_00),numpy.array(myY_00),'.')
    plt.plot(numpy.array(myX_10),numpy.array(myY_10),'.')
    plt.plot(numpy.array(myX_01),numpy.array(myY_01),'.')
    plt.plot(numpy.array(myX_11),numpy.array(myY_11),'.')
    plt.grid(color = 'gray')
    plt.xlabel('X [mm]')
    plt.ylabel('-Y [mm]')
