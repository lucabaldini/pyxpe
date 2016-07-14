# script designed for plotting result of long run with Fe55 in Rome
# for flat field, gain map study etc.
import numpy as np
import ROOT
import os, sys
ROOT.gStyle.SetOptFit(1)

CUT = "(fNClusters==1)"# "(1)"#

TH = 9
t = ROOT.TChain("tree")
t.Add('/data/xpe/xpedata/101_0002936/101_0002936_data_TH%d.root' % TH)
t.Add('/data/xpe/xpedata/101_0002937/101_0002937_data_TH%d.root' % TH)
t.Add('/data/xpe/xpedata/101_0002938/101_0002938_data_TH%d.root' % TH)
print "Getting Recon with TH %d for total number of events %d, with cut:\n %s" %\
    (TH, t.GetEntries(CUT), CUT)

out_filepath = "./tmpGainMap.root" # None to avoid it
#"GainMap_TH9_clu1_baricenter.root" #
if out_filepath !=None:
    out_file = ROOT.TFile(out_filepath, "RECREATE")

def getGainRes(histo, nsigfit= 1.8 , debug = False):
    """ Function to get gain and resolution
    from a fPHeight distribution
    """
    #m = histo.GetBinCenter(histo.GetMaximumBin())
    m = histo.GetMean()
    r = histo.GetRMS()
    if debug:
        print "First gaus fit around %.1f +- %.1f" % (m,r)
    g = ROOT.TF1("g", "gaus", m-r, m+1.2*r)
    histo.Fit("g", "RNQ")
    g1 = g.GetParameter(1)
    g2 = g.GetParameter(2)
    g = ROOT.TF1("g", "gaus", g1 -nsigfit*g2, g1 +nsigfit*g2)
    if debug:
        print "Second gaus fit in range [%.1f ; %.1f]" %\
            (g1 -nsigfit*g2, g1 +nsigfit*g2)
    histo.Fit("g", "RQ")
    gPeak  = g.GetParameter(1)
    gSigma = g.GetParameter(2)
    if debug:
        ctmp = ROOT.TCanvas()
        histo.Draw()
        print "Gain %.1f Res sigma %.3f %%" % (gPeak, gSigma/gPeak)
        raw_input("E2Go")
    return (gPeak, gSigma/gPeak)

# Setting binning for gain scan
#gain_binx = np.append(np.append([-7.6], np.linspace(-7, 7,35+1)), [7.6]) 
#gain_biny = np.append(np.append([-7.5], np.linspace(-7, 7,35+1)), [7.5])
#bin_pad = 0
#gain_binx = np.append(np.append([-7.6], np.linspace(-7, 7,4*35+1)), [7.6]) 
#gain_biny = np.append(np.append([-7.5], np.linspace(-7, 7,4*35+1)), [7.5])
#bin_pad = 2
gain_binx = np.append(np.append([-7.6], np.linspace(-6, 6,12+1)), [7.6]) 
gain_biny = np.append(np.append([-7.5], np.linspace(-6, 6,12+1)), [7.5])
bin_pad = 0
gain_binx_n = len(gain_binx) -1
gain_biny_n = len(gain_biny) -1
gain_binph = np.linspace(0, 10000, 200+1)


cGainTmp = ROOT.TCanvas()
hGainMapPH = ROOT.TH3F("hGainMapPH", "X Y PH",\
                       gain_binx_n, gain_binx,\
                       gain_biny_n, gain_biny,\
                       len(gain_binph)-1, gain_binph)
hGainMap   = ROOT.TH2F("hGainMap", "Gain (adc @ 5.9 keV);X [mm];Y [mm]",\
                       gain_binx_n, gain_binx,\
                       gain_biny_n, gain_biny)
hResMap    = ROOT.TH2F("hResMap", "Resolution (sigma);X [mm];Y [mm]",\
                       gain_binx_n, gain_binx,\
                       gain_biny_n, gain_biny)
hNumEvtMap = ROOT.TH2F("hNumEvtMap", "NumEvtMap;X [mm];Y [mm]",\
                       gain_binx_n, gain_binx,\
                       gain_biny_n, gain_biny)
gGainResCorr = ROOT.TGraph((gain_binx_n-2)*(gain_binx_n-2))
gn = 0
#t.Project("hGainMapPH", "fPHeight[0]:fImpactY[0]:fImpactX[0]", CUT)
t.Project("hGainMapPH", "fPHeight[0]:fBaricenterY[0]:fBaricenterX[0]", CUT)
for j in xrange(gain_biny_n):
    print "j=", j
    for i in xrange(gain_binx_n):
        htmp = hGainMapPH.ProjectionZ("htmp_%d_%d" %(i,j),
                                      max(1, (i+1)-bin_pad),
                                      min((i+1)+bin_pad, gain_binx_n),
                                      max(1, (j+1)-bin_pad),
                                      min((j+1)+bin_pad, gain_biny_n)
                                      )
        #print i, j,
        #print "\tx:", max(1, (i+1)-bin_pad), "-", i+1,"-",min((i+1)+bin_pad, gain_binx_n),
        #print "\ty:", max(1, (j+1)-bin_pad),  "-", j+1,"-",min((j+1)+bin_pad, gain_biny_n)
        if out_filepath !=None:
            htmp.Write()
        n = htmp.GetEntries()
        hNumEvtMap.SetBinContent(i+1, j+1, n)
        (G, R) = getGainRes(htmp, 1.8, False)
        hGainMap.SetBinContent(i+1, j+1, G)
        hResMap.SetBinContent(i+1, j+1, R)
        if i>0 and i<(gain_binx_n-1) and j>0 and j<(gain_biny_n-1):
            gGainResCorr.SetPoint(gn, G, R)
            gn+=1

cGainResMap = ROOT.TCanvas("cGainResMap", "cGainResMap", 1200, 600)
cGainResMap.Divide(2,1)
cGainResMap.cd(1)
hGainMap.SetStats(0)
hGainMap.Draw("colz")
hGainMap.GetZaxis().SetRangeUser(2800, 4600)
hGainMap.GetZaxis().SetLabelSize(0.03)
cGainResMap.cd(2)
hResMap.SetStats(0)
hResMap.Draw("colz")
hResMap.GetZaxis().SetRangeUser(0.05, 0.25)
hResMap.GetZaxis().SetLabelSize(0.03)

cGainResMap1 = ROOT.TCanvas("cGainResMap1", "cGainResMap1", 1200, 600)
cGainResMap1.Divide(2,1)
cGainResMap1.cd(1)
hNumEvtMap.SetStats(0)
hNumEvtMap.Draw("colz")
cGainResMap1.cd(2)
gGainResCorr.SetMarkerStyle(7)
gGainResCorr.Draw("AP")
gGainResCorr.SetTitle("")
gGainResCorr.GetXaxis().SetTitle("Gain (adc)")
gGainResCorr.GetYaxis().SetTitle("Resolution (sigma)")
ROOT.gPad.SetGridx(True)
ROOT.gPad.SetGridy(True)

ROOT.gPad.Update()

if out_filepath !=None:
    hNumEvtMap.Write()
    hGainMap.Write()
    hResMap.Write()
    gGainResCorr.Write()
    out_file.Write()
    #out_file.Close()

