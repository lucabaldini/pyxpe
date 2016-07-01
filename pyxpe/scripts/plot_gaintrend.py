# standalone scipt to plot gain and resolution vs time.
# Designed for the specific purpose of GPD018 trend after refill

import ROOT
import numpy as np
ROOT.gStyle.SetPadGridX(True)
ROOT.gStyle.SetPadGridY(True)
G2FWHM = 2.3548200450309493

tt = ROOT.TChain("tree")
# SELECT RUNS to be analized
for runid in [393, 394, 395, 396, 397, 398]:
    tt.Add("/data0/xpe/xpedata/002_%07d/002_%07d_data_TH5.root" %(runid, runid))

# SELECT OPTIONS
Label      = 'GPD018_fill1'
Nbins      = 20
nsigfit    = 1.5
CUT        = "(fNClusters==1)"
dh         = ROOT.TDatime(2016,6,30,10,56,21) # started on Jun 30 10:56:21 2016
TimeOffset = dh.Convert()
TimeHours  = "((fTimeStamp-%.1f)/3600.)" % TimeOffset
useFWHM      = False
useSmallSpot = True
SmallSpotCut = "((fImpactY[0]-0)**2 + (fImpactX[0]-0.5)**2)<(1**2)"

minTime   = (tt.GetMinimum("fTimeStamp")-TimeOffset)/3600.
maxTime   = (tt.GetMaximum("fTimeStamp")-TimeOffset)/3600.
print "Eval trend from %f hr to %f hr after start" % (minTime, maxTime)

timeBins = np.linspace(0,maxTime, Nbins+1)
timeVal = 0.5*(timeBins[1:] +timeBins[:-1])
timeErr = 0.5*(timeBins[1:] -timeBins[:-1])
gainVal = np.zeros(Nbins)
gainErr = np.zeros(Nbins)
resVal  = np.zeros(Nbins)
resErr  = np.zeros(Nbins)

ctmp = ROOT.TCanvas()
htmp = ROOT.TH1F("htmp", "htmp", 200, 0, 10000)
g0 = ROOT.TF1("g0", "gaus", 3000, 4000)
for i in xrange(Nbins):
    TimeCut = "(%s>=%f && %s <%f)" % \
              (TimeHours, timeBins[i], TimeHours, timeBins[i+1])
    thecut =  CUT + "&&"+ TimeCut
    if useSmallSpot:
        thecut +=  " && " + SmallSpotCut
    print i, thecut
    tt.Project("htmp", "fPHeight[0]", thecut)
    htmp.Draw()
    htmp.Fit("g0", "RNQ")
    g = ROOT.TF1("g", "gaus", \
                 g0.GetParameter(1) - nsigfit*g0.GetParameter(2),\
                 g0.GetParameter(1) + nsigfit*g0.GetParameter(2))
    htmp.Fit("g", "R")
    gPeak  = g.GetParameter(1)
    gSigma = g.GetParameter(2)
    print gPeak, gSigma, "Eres=", gSigma/gPeak, "FWHM", 2.355*(gSigma/gPeak)
    gainVal[i] = gPeak
    gainErr[i] = g.GetParError(1)
    resVal[i]  = 100.*gSigma/gPeak
    #resErr[i]  = resVal[i]*np.sqrt((g.GetParError(1)/gPeak)**2 +\
    #                               (g.GetParError(2)/gSigma )**2)
    resErr[i]  = resVal[i]*(g.GetParError(2)/gSigma ) # ~ok for now
    if useFWHM:
        resVal[i] = resVal[i]*G2FWHM
        resErr[i] = resErr[i]*G2FWHM

    ROOT.gPad.Update()
    #raw_input("inspect and pres enter")

cTrend = ROOT.TCanvas("gaintrend_%s" %Label,"gaintrend_%s" %Label,1000,700)
cTrend.Divide(1,2)
cTrend.cd(1)
gGain = ROOT.TGraphErrors(Nbins, timeVal*3600, gainVal, timeErr*3600, gainErr)
gGain.SetMarkerStyle(20)
gGain.SetTitle("Gain (gaussian peak) @ 5.9 keV (Fe55)")
gGain.GetYaxis().SetTitle("Peak")
gGain.GetXaxis().SetTitle("")
gGain.GetXaxis().SetTimeDisplay(1);
gGain.GetXaxis().SetTimeFormat("#splitline{%d/%m/%y}{%H:%M:%S}");
gGain.GetXaxis().SetTimeOffset(TimeOffset)
gGain.GetXaxis().SetLabelOffset(0.03)
gGain.Draw("ap")
cTrend.cd(2)
gRes = ROOT.TGraphErrors(Nbins, timeVal, resVal, timeErr, resErr)
gRes.SetMarkerStyle(20)
gRes.SetTitle("Resolution (gaussian sigma/peak) @ 5.9 keV (Fe55)")
if useFWHM:
    gRes.SetTitle("Resolution (gaussian sigma FWHM) @ 5.9 keV (Fe55)")
gRes.GetXaxis().SetTitle("Elapsed Time (hours)")
gRes.GetYaxis().SetTitle("#Delta E/E (%)")
gRes.Draw("ap")

