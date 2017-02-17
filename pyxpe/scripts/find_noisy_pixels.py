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
    Finding the triggering minicluster and 
    selecting its pixel with highest pulse height 
    """
    if event.xmax-event.xmin > 17 or event.ymax-event.ymin > 21:
        # window is too large and triggering pixels is too ambiguous
        return (-1, -1)
    # triggering minicluster ids
    if event.xmin>0:
        x0 = 8
    else:
        x0 = event.xmax-9
    if event.ymin>0:
        y0 = 10
    else:
        y0 = event.ymax-11
    submatrix = event.adc_values[x0:x0+2, y0:y0+2]
    x, y = np.unravel_index(np.argmax(submatrix),(2,2))
    return (x+x0+event.xmin, y+y0+event.ymin)
    

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


#runId = 670 # ci(100,100)
#runId = 679 # ci(9,101)
#runId = 680 # ci(8,101)
#runId = 681 # ci(7,101)
#runId = 682 # ci(2,101)
#runId = 683 # ci(3,101)
#runId = 684 # ci(0,101)
#runId = 685 # ci(298,10)
#runId = 686 # ci(298,3)
#runId = 687 # ci(298,350)
#runId = 688 # ci(0,0)
#file_path = '/data/xpedata/001/001_%07d/001_%07d_data.mdat'% (runId,runId)
input_file = xpeBinaryFileWindowed(file_path)

# matrix for pixel occupancy
# DON"T FORGET: np.zeros([row,col])
pxl_occupancy = np.zeros([XPOL_NUM_ROWS,XPOL_NUM_COLUMNS])

for i in range(N_MAX):
    event = input_file.next()
    #event.draw_ascii(ZERO_SUPP)
    #print (event)
    #print  guess_trg_pxl(event)
    #col, row = pixels_overthr(event, ZERO_SUPP)
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
    

