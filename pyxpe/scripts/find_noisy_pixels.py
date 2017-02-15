import numpy as np
import matplotlib.pyplot as plt
from pyxpe.recon.xpol import XPOL_NUM_COLUMNS, XPOL_NUM_ROWS
from pyxpe.recon.binio import xpeBinaryFileWindowed
#from pyxpe.recon.event import xpeEventWindowed
#from pyxpe.recon.xpol import XPOL_PIXELS_PER_BUFFER, XPOL_NUM_BUFFERS

N_MAX     = 1000
ZERO_SUPP = 9

def pixels_overthr(event, thr=9):
    """ Get all pixels above threshold
    """
    _mask = event.adc_values > thr
    adc_values = event.adc_values[_mask]
    col, row = np.where(_mask)
    return (event.xmin + col, event.ymin + row)

def eval_trg_pxl():
    """ Function to eval pixels that likely triggered the event.
    TBD
    """
    pass

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
input_file = xpeBinaryFileWindowed(file_path)

# matrix for pixel occupancy
# DON"T FORGET: np.zeros([row,col])
pxl_occupancy = np.zeros([XPOL_NUM_ROWS,XPOL_NUM_COLUMNS])

for i in range(N_MAX):
    event = input_file.next()
    #event.draw_ascii(ZERO_SUPP)
    col, row = pixels_overthr(event, ZERO_SUPP)
    for (c,r) in zip(col, row):
        pxl_occupancy[r][c] +=event.adc_value(c-event.xmin, r-event.ymin)
        pxl_occupancy[r][c] += 1
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
    

