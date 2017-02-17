import numpy as np
import matplotlib.pyplot as plt
from pyxpe.recon.xpol import XPOL_NUM_COLUMNS, XPOL_NUM_ROWS
from pyxpe.recon.binio import xpeBinaryFileWindowed
#from pyxpe.recon.event import xpeEventWindowed
#from pyxpe.recon.xpol import XPOL_PIXELS_PER_BUFFER, XPOL_NUM_BUFFERS

N_MAX     = 1000
ZERO_SUPP = 5

def pixels_overthr(event, thr=9):
    """ Get all pixels above threshold
    """
    _mask = event.adc_values > thr
    adc_values = event.adc_values[_mask]
    col, row = np.where(_mask)
    return (event.xmin + col, event.ymin + row)

def guess_trg_pxl(event):
    """ Function to eval pixels that likely triggered the event.
    TBD
    """
    xtrg = -1
    ytrg = -1
    if event.ymax-event.ymin == 21 and event.xmax-event.xmin == 17:
        # ~center of the asic
        submatrix = event.adc_values[8:10, 10:12]
        x, y = np.unravel_index(np.argmax(submatrix),(2,2))
        return (x+8+event.xmin, y+10+event.ymin)

    return (xtrg, ytrg)

def highest_occupancy_pxl(pxl_matrix, npxl=10):
    """ Get the list (col, row, occupancy) of the
    npxl pixels with the highest occpancy
    """
    i = 0
    outlist = []
    while i<npxl:
        row_max, col_max = np.unravel_index(pxl_matrix.argmax(),
                                            pxl_matrix.shape)
        outlist.append((col_max, row_max, pxl_matrix[row_max][col_max]))
        pxl_matrix[row_max][col_max] = 0
        i+=1
    return outlist


file_path = '/data/xpedata/001/001_0000669/001_0000669_data.mdat'

#file_path = '/data/xpedata/001/001_0000670/001_0000670_data.mdat' # ci(100,100)
#file_path = '/data/xpedata/001/001_0000671/001_0000671_data.mdat' # ci(9,100)
#file_path = '/data/xpedata/001/001_0000672/001_0000672_data.mdat' # ci(10,100)
#file_path = '/data/xpedata/001/001_0000676/001_0000676_data.mdat' # ci(11,100)
#file_path = '/data/xpedata/001/001_0000677/001_0000677_data.mdat' # ci(101,100)
#file_path = '/data/xpedata/001/001_0000678/001_0000678_data.mdat' # ci(101,101)
input_file = xpeBinaryFileWindowed(file_path)

# matrix for pixel occupancy
# DON"T FORGET: np.zeros([row,col])
pxl_occupancy = np.zeros([XPOL_NUM_ROWS,XPOL_NUM_COLUMNS])

for i in range(N_MAX):
    event = input_file.next()
    #event.draw_ascii(ZERO_SUPP)
    #print (event)
    #print event.xmin+8, event.xmax-9,"-",event.ymin+10, event.ymax-11,\
    #    "---", guess_trg_pxl(event)
    col, row = pixels_overthr(event, ZERO_SUPP)
    #for (c,r) in zip(col, row):
    #    pxl_occupancy[r][c] +=event.adc_value(c-event.xmin, r-event.ymin)
    #    pxl_occupancy[r][c] += 1
    col, row = guess_trg_pxl(event)
    pxl_occupancy[row][col] += 1
    #event.draw(ZERO_SUPP)

from copy import copy
cp_occ = copy(pxl_occupancy)

print highest_occupancy_pxl(pxl_occupancy, 10)

# OK, here my debug in ROOT
import ROOT
h = ROOT.TH2F("h","", 300,0,300,352,0,352)
for j in xrange(300):
    for i in xrange(352):
        h.SetBinContent(j+1,i+1,cp_occ[i][j])

h.Draw("colz")
    

