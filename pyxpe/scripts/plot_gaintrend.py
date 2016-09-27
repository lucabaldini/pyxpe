# standalone scipt to plot gain and resolution vs time.
# Designed for the specific purpose of GPD018 trend after refill

import ROOT
import numpy as np
ROOT.gStyle.SetOptFit(1)
ROOT.gStyle.SetPadGridX(True)
ROOT.gStyle.SetPadGridY(True)
G2FWHM = 2.3548200450309493

tt = ROOT.TChain("tree") 
# SELECT RUNS to be analized
# GPD18
for runid in [393, 394, 395, 396, 397]:#, 398, 399, 400,
              #401, 402, 403, 406, 407, 408, 409, 410, 
              #411, 412, 413, 414, 415, 416,
              #417, 418, 420, 421, 422, 423, 424, 425, 426, 427, 428,
              #432, 433, 434, 435, 436, 437, 438, 439, 440]:
    #tt.Add("/data1/xpe/xpedata/002_%07d/002_%07d_data_TH5.root"%(runid, runid))
    pass
# Note:
# 1 night gap between run 400 (timeout at 20:11 1/7) and 402 (start 9:37 2/7)
# 402 ends with timeout, 403 started after ~1 hr -- 403 is very small
# 1 night gap between run 416 (timeout at 21:15 4/7) and 417 (start 9:13 5/7)
# 418 interrupted by mistake (low voltage turnoff), skip 419.
# CHANGE SOURCE POSITION since run 432 (fTimeStamp threshold 1467898100.)
# 432 ends with timeout,  about 1.5hr gap (16:22 - 17:50) before 434 starts

# SELECT RUNS to be analized
# GPD19
for runid in [2724, 2725, 2726, 2727, 2728, 2729, 2730]:
    tt.Add("/data1/xpe/xpedata/002_%07d/002_%07d_data_TH5.root"%(runid, runid))
    pass
# SELECT OPTIONS
#Label      = 'GPD018_fill1_reference'
Label      = 'GPD019_fill0'
TimeBin    = 1. #hours
nsigfit    = 1.5
#dh         = ROOT.TDatime(2016,6,30,10,56,21) # started on Jun 30 10:56:21 2016
dh         = ROOT.TDatime(2016,9,26,10,59,0) # started on Mon Sep 26 10:59:13 2016
TimeOffset = dh.Convert()
#CUT        = '(fTrigWindow > 600*(1 + (fTimeStamp - %f)*0.5e-03/4000.))' % TimeOffset
CUT = '(1)'
TimeHours  = "((fTimeStamp-%.1f)/3600.)" % TimeOffset
minEvtInHist = 500
useFWHM      = True
useSmallSpot = False
SpotCut_tthr = 1467898100
SpotCut_0    = "((fImpactY[0]-0)**2 + (fImpactX[0]-0.5)**2)<(1.**2)" # max 2.5 to avoid bad chan area
SpotCut_1    = "((fImpactY[0]+4.)**2 + (fImpactX[0]+4.)**2)<(1.**2)" #
SmallSpotCut = "((fTimeStamp<=%d)*(%s) + (fTimeStamp>%d)*(%s))" %\
               (SpotCut_tthr, SpotCut_0, SpotCut_tthr, SpotCut_1)

"""
# time bins logic has to take into account time threshold
if tt.GetMaximum("fTimeStamp")>=SpotCut_tthr:
    split_time_flag = True
    minTime    = (tt.GetMinimum("fTimeStamp")-TimeOffset)/3600.
    maxTime    = (SpotCut_tthr-TimeOffset)/3600.
    Nbins      = int((maxTime-minTime)/TimeBin)
    timeBins_0 = np.linspace(0,maxTime, Nbins+1)
    minTime    = (SpotCut_tthr-TimeOffset)/3600.
    maxTime    = (tt.GetMaximum("fTimeStamp")-TimeOffset)/3600.
    Nbins      = max(round((maxTime-minTime)/TimeBin),1)
    timeBins_1 = np.linspace(minTime, maxTime, Nbins+1)
    timeBins   = np.append(timeBins_0, timeBins_1[1:])
    print (SpotCut_tthr-TimeOffset)/3600.
    print (tt.GetMaximum("fTimeStamp")-TimeOffset)/3600.
    print (maxTime-minTime)/TimeBin
    #raise
else:
    split_time_flag = False
    minTime    = (tt.GetMinimum("fTimeStamp")-TimeOffset)/3600.
    maxTime    = (tt.GetMaximum("fTimeStamp")-TimeOffset)/3600.
    Nbins      = int((maxTime-minTime)/TimeBin)
    timeBins   = np.linspace(0,maxTime, Nbins+1)

print SpotCut_tthr, (SpotCut_tthr-TimeOffset)/3600.
print timeBins
"""
split_time_flag = False
minTime    = (tt.GetMinimum("fTimeStamp")-TimeOffset)/3600.
maxTime    = (tt.GetMaximum("fTimeStamp")-TimeOffset)/3600.
Nbins      = max(int((maxTime-minTime)/TimeBin), 1)
timeBins   = np.linspace(0,maxTime, Nbins+1)
print "Eval trend from %f hr to %f hr after start in %d bins" %\
    (min(timeBins), max(timeBins), len(timeBins))

Nbins    = len(timeBins)-1
timeVal  = 0.5*(timeBins[1:] +timeBins[:-1])
timeErr  = 0.5*(timeBins[1:] -timeBins[:-1])
timeX    = []
timeXErr = []
gainVal  = []
gainErr  = []
resVal   = []
resErr   = []
chbins   = np.linspace(0,10000,201)

# project entire tree in a 2D histo
hAll = ROOT.TH2F("hAll", "hAll", len(chbins)-1, chbins, \
                 len(timeBins)-1, timeBins)
#hAllNP = ROOT.TH2F("hAllNP", "hAllNP", len(chbins)-1, chbins, \
#                   len(timeBins)-1, timeBins)

thecut =  CUT
if useSmallSpot:
    thecut +=  " && " + SmallSpotCut
print "Projection tree with cut:\n%s" % thecut
tt.Project("hAll", "%s:fPHeight[0]" % TimeHours, thecut)
#tt.Project("hAllNP", "%s:fPHeight[0]" % TimeHours, SmallSpotCut)

ctmp = ROOT.TCanvas()
ctmp.Print("tmpAllHisto.ps[");
hAll.Draw("colz")
ctmp.Print("tmpAllHisto.ps");
g0 = ROOT.TF1("g0", "gaus", 3000, 5500)
for i in xrange(Nbins):
    # get projection in time (y coord in hAll)
    htmp = hAll.ProjectionX("htmp", i+1,i+1)
    #htmpnp = hAllNP.ProjectionX("htmpnp", i+1,i+1)
    htmp.Draw()
    if htmp.GetEntries() >= minEvtInHist:
        #print '---> ', htmpnp.GetEntries(), htmp.GetEntries(),\
        #    htmp.GetEntries()/htmpnp.GetEntries()
        
        htmp.Fit("g0", "RNQ")
        g = ROOT.TF1("g", "gaus", \
                     g0.GetParameter(1) - nsigfit*g0.GetParameter(2),\
                     g0.GetParameter(1) + nsigfit*g0.GetParameter(2))
        htmp.Fit("g", "RQ")
        timeX.append(timeVal[i])
        timeXErr.append(timeErr[i])
        gPeak  = g.GetParameter(1)
        gSigma = g.GetParameter(2)
        print  "\n",i, gPeak, gSigma, "Eres=", gSigma/gPeak, \
            "FWHM", 2.355*(gSigma/gPeak)
        gainVal.append(gPeak)
        gainErr.append(g.GetParError(1))
        resVal.append(100.*gSigma/gPeak)
        resErr.append(resVal[-1]*(g.GetParError(2)/gSigma )) # ~ok for now
        if useFWHM:
            resVal[-1] = resVal[-1]*G2FWHM
            resErr[-1] = resErr[-1]*G2FWHM

    ROOT.gPad.Update()
    ctmp.Print("tmpAllHisto.ps");
    #raw_input("inspect and pres enter")

ctmp.Print("tmpAllHisto.ps]");
timeX    = np.array(timeX)
timeXErr = np.array(timeXErr)
gainVal  = np.array(gainVal)
gainErr  = np.array(gainErr)
resVal   = np.array(resVal)
resErr   = np.array(resErr)
N        = len(timeX)
cTrend = ROOT.TCanvas("gaintrend_%s" %Label+"_smallspot"*useSmallSpot,\
                      "gaintrend_%s" %Label+"_smallspot"*useSmallSpot,1000,700)
cTrend.Divide(1,2)
cTrend.cd(1)
gGain = ROOT.TGraphErrors(N, timeX*3600, gainVal, timeXErr*3600, gainErr)
gGain.SetMarkerStyle(20)
gGain.SetTitle("Gain (gaussian peak) @ 5.9 keV (Fe55)")
gGain.GetYaxis().SetTitle("Peak")
gGain.GetXaxis().SetTitle("")
gGain.GetXaxis().SetTimeDisplay(1);
gGain.GetXaxis().SetTimeFormat("#splitline{%d/%m/%y}{%H:%M:%S}");
gGain.GetXaxis().SetTimeOffset(TimeOffset)
gGain.GetXaxis().SetLabelOffset(0.03)
gGain.Draw("ap")
#fGain = ROOT.TF1('fGain', 'pol1', 0, (SpotCut_tthr-TimeOffset))
fGain = ROOT.TF1('fGain', 'pol1', 0, maxTime)
fGain.SetParameters(10, 0.002)
#gGain.Fit("fGain", "R")
ROOT.gPad.Update()
#psGain = gGain.GetListOfFunctions().FindObject("stats");
#psGain.SetX1NDC(0.65);
#psGain.SetY1NDC(0.18);
#psGain.SetX2NDC(0.87);
#psGain.SetY2NDC(0.34);
ROOT.gPad.Update()
if split_time_flag:
    time_bound_line = ROOT.TLine()
    time_bound_line.SetLineColor(ROOT.kBlue)
    time_bound_line.SetLineStyle(2)
    time_bound_line.DrawLine((SpotCut_tthr-TimeOffset),
                             gGain.GetYaxis().GetXmin(),
                             (SpotCut_tthr-TimeOffset),
                             gGain.GetYaxis().GetXmax())
    ROOT.gPad.Update()
                             

cTrend.cd(2)
gRes = ROOT.TGraphErrors(N, timeX, resVal, timeXErr, resErr)
gRes.SetMarkerStyle(20)
gRes.SetTitle("Resolution (gaussian sigma/peak) @ 5.9 keV (Fe55)")
if useFWHM:
    gRes.SetTitle("Resolution (gaussian sigma FWHM) @ 5.9 keV (Fe55)")
gRes.GetXaxis().SetTitle("Elapsed Time (hours)")
gRes.GetYaxis().SetTitle("#Delta E/E (%)")
gRes.Draw("ap")
fRes = ROOT.TF1('fRes', 'pol1', 0, (SpotCut_tthr-TimeOffset)/3600.)
fRes.SetParameters(0.1, 0.001)
#gRes.Fit("fRes", "R")
ROOT.gPad.Update()
#psRes = gRes.GetListOfFunctions().FindObject("stats");
#psRes.SetX1NDC(0.33);
#psRes.SetY1NDC(0.64);
#psRes.SetX2NDC(0.54);
#psRes.SetY2NDC(0.80);
cTrend.Update()
if split_time_flag:
    time_bound_line.DrawLine((SpotCut_tthr-TimeOffset)/3600.,
                             gRes.GetYaxis().GetXmin(),
                             (SpotCut_tthr-TimeOffset)/3600.,
                             gRes.GetYaxis().GetXmax())
    ROOT.gPad.Update()
