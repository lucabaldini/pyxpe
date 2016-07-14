# standalone scipt to plot gain vs position for attach. coeff.
# Designed for the specific purpose of GPD018 trend after refill

import ROOT
import numpy as np
ROOT.gStyle.SetOptFit(1)
ROOT.gStyle.SetPadGridX(True)
ROOT.gStyle.SetPadGridY(True)

tt = ROOT.TChain("tree") 
for runid in [2915, 2917]:
    tt.Add("/data/xpe/xpedata/101_%07d/101_%07d_data_TH9.root" %(runid, runid))
    print "adding run %d, entries %d" % (runid, tt.GetEntries())

thv = ROOT.TChain("tree")
thv.Add("/data/xpe/xpedata/101_0002945/101_0002945_data_TH9.root")

GainNormList = [2921, 2922, 2923, 2924, 2925, 2926, 2927, 2928, 2930, 2931]
#    (2921, 0.), (2922, 0.5), (2923, 1.), (2924, 1.5),\
#                (2925, 2.), (2926, 2.5), (2927, 3.), (2928, 3.5),\
#                (2930, 4.), (2931, 4.5) ]
    
Label = 'GPD018_fill1_attach'
Xmin  = -4.75#-6.25
Xmax  = 0.25#1.25
Xbins = 10 #15
CUT   = "(1)"
minEvtInHist = 2000
nsigfit      = 1.

c = ROOT.TCanvas()
c.Divide(2,2)
c.cd(1)
tt.Draw("fImpactY[0]", CUT)
thv.Draw("fImpactY[0]", CUT, "sames")
c.cd(2)
tt.Draw("fImpactY[0]:fImpactX[0]", CUT)
c.cd(4)
tt.Draw("fImpactX[0]", CUT)
thv.Draw("fImpactX[0]", CUT, "sames")

cc = ROOT.TCanvas("cGvsX", "cGvsX", 800, 600)
cc.Divide(1,3)
cc.cd(1)
hGvsX = ROOT.TH2F("hGvsX", "PH vs X;Impact position X [mm];PHeight [adc]",
                  Xbins,Xmin,  Xmax, 200,0, 10000)
tt.Project("hGvsX", "fPHeight[0]:fImpactX[0]", CUT)
hGvsX.Draw("colz")
cc.cd(2)
hGvsXhv = ROOT.TH2F("hGvsXhv", "PH vs X (HV);Impact position X [mm];PHeight [adc]",
                  Xbins,Xmin,  Xmax, 200,0, 10000)
thv.Project("hGvsXhv", "fPHeight[0]:fImpactX[0]", CUT)
hGvsXhv.Draw("colz")

cc.cd(3)
hGvsX_peak = ROOT.TH1F("hGvsX_peak", "PH peak vs X;Impact position X [mm];PH peak [adc]",Xbins,Xmin,  Xmax)
hGvsX_peak.SetLineColor(ROOT.kRed)
hGvsX_peak.SetLineWidth(2)
hGvsXhv_peak = ROOT.TH1F("hGvsXhv_peak", "PH peak vs X (HV);Impact position X [mm];PH peak [adc]",Xbins,Xmin,  Xmax)
hGvsXhv_peak.SetLineColor(ROOT.kGreen)
hGvsXhv_peak.SetLineWidth(2)
g0 = ROOT.TF1("g0", "gaus", 2500, 3500)
for i in xrange(Xbins):
    # get projection in time (y coord in hAll)
    htmp = hGvsX.ProjectionY("htmp", i+1,i+1)
    htmp.Draw()
    print i, htmp.GetEntries()
    if htmp.GetEntries() >= minEvtInHist:
        #print '---> ', htmpnp.GetEntries(), htmp.GetEntries(),\
        #    htmp.GetEntries()/htmpnp.GetEntries()
        
        htmp.Fit("g0", "RNQ")
        g = ROOT.TF1("g", "gaus", \
                     g0.GetParameter(1) - nsigfit*g0.GetParameter(2),\
                     g0.GetParameter(1) + nsigfit*g0.GetParameter(2))
        htmp.Fit("g", "RQ")
        gPeak  = g.GetParameter(1)
        gSigma = g.GetParameter(2)
        print i, gPeak, gSigma, "Eres=", gSigma/gPeak
        hGvsX_peak.SetBinContent(i+1, gPeak)
        hGvsX_peak.SetBinError(i+1, g.GetParError(1))
    else:
        hGvsX_peak.SetBinContent(i+1, 0)
        hGvsX_peak.SetBinError(i+1, 0)

    # do it again for test with HV DRIFT 3100
    htmp = hGvsXhv.ProjectionY("htmp", i+1,i+1)
    htmp.Draw()
    print "HV:", i, htmp.GetEntries()
    if htmp.GetEntries() >= minEvtInHist:
        #print '---> ', htmpnp.GetEntries(), htmp.GetEntries(),\
        #    htmp.GetEntries()/htmpnp.GetEntries()
        
        htmp.Fit("g0", "RNQ")
        g = ROOT.TF1("g", "gaus", \
                     g0.GetParameter(1) - nsigfit*g0.GetParameter(2),\
                     g0.GetParameter(1) + nsigfit*g0.GetParameter(2))
        htmp.Fit("g", "RQ")
        gPeak  = g.GetParameter(1)
        gSigma = g.GetParameter(2)
        print i, gPeak, gSigma, "Eres=", gSigma/gPeak
        hGvsXhv_peak.SetBinContent(i+1, gPeak)
        hGvsXhv_peak.SetBinError(i+1, g.GetParError(1))
    else:
        hGvsXhv_peak.SetBinContent(i+1, 0)
        hGvsXhv_peak.SetBinError(i+1, 0)
    
    ROOT.gPad.Update()
    #raw_input("inspect and pres enter")

# get Gain map from single vertical runs
hphtmp = ROOT.TH1F("hphtmp", "hphtmp", 200,0,10000)
gGainRef = ROOT.TGraph(len(GainNormList))
gGainRef.SetMarkerStyle(22)
hGvsX_ref = ROOT.TH1F("hGvsX_ref", "PH ref vs X;Impact position X [mm];PH peak [adc]", 10, -4.75, 0.25) # non TOCCARE!!!
for (i, runid) in enumerate(GainNormList):
    #f = ROOT.TFile("/data/xpe/xpedata/101_%07d/101_%07d_data_TH9.root" %(runid, runid))
    #myt = f.Get("tree")
    myt = ROOT.TChain("tree")
    myt.Add("/data/xpe/xpedata/101_%07d/101_%07d_data_TH9.root" %(runid, runid))
    myt.Project("hphtmp", "fPHeight[0]", CUT)
    hphtmp.Draw()
    hphtmp.Fit("g0", "RNQ")
    g = ROOT.TF1("g", "gaus", \
                 g0.GetParameter(1) - nsigfit*g0.GetParameter(2),\
                 g0.GetParameter(1) + nsigfit*g0.GetParameter(2))
    hphtmp.Fit("g", "RQ")
    gPeak  = g.GetParameter(1)
    gSigma = g.GetParameter(2)
    print "REF GAIN:", i, gPeak, gSigma, "Eres=", gSigma/gPeak
    myt.Draw("fImpactX[0]>>hpostmp", CUT)
    x = ROOT.hpostmp.GetMean()
    gGainRef.SetPoint(i, x, gPeak)
    hGvsX_ref.Fill(x, gPeak)
    ROOT.gPad.Update()
    #raw_input("e2g")
        
#hGvsX_pf = hGvsX.ProfileX()
#hGvsX_pf.SetLineColor(ROOT.kBlack)
#hGvsX_pf.SetLineWidth(2)
#hGvsX_pf.Draw()
hGvsX_peak.Draw("HE")
hGvsXhv_peak.Draw("HE,sames")
gGainRef.Draw("P,same")
hGvsX_ref.Draw("sames")

f1 = ROOT.TF1("f1", "pol1", -3.5, 0.5)
#hGvsX_pf.Fit("f1", "R")


# una merda!!! va fatto bene!!!
gGainNorm = ROOT.TGraph(Xbins)
gGainNorm.SetMarkerStyle(20)
gGainNormhv = ROOT.TGraph(Xbins)
gGainNormhv.SetMarkerStyle(22)
gGainNormhv.SetMarkerColor(hGvsXhv_peak.GetLineColor())
if Xbins != len(GainNormList):
    print "AHHHHH!!!!"
for i in range(Xbins):
    mygain   = hGvsX_ref.GetBinContent(i+1)
    myattach = hGvsX_peak.GetBinContent(i+1)
    myx      = hGvsX_peak.GetBinCenter(i+1)
    print "NORMGAIN", i, myx, myattach, mygain, myattach/mygain
    gGainNorm.SetPoint(i, myx, myattach/mygain)
    gGainNormhv.SetPoint(i, myx, hGvsXhv_peak.GetBinContent(i+1)/mygain)

#cc.cd(3)
cGainVsPosition = ROOT.TCanvas("cGainVsPosition", "cGainVsPosition", 800, 600)
cGainVsPosition.Divide(1,2)
cGainVsPosition.cd(1)
hGvsX_peak.SetStats(0)
hGvsX_ref.SetStats(0)
hGvsX_peak.Draw()
hGvsX_ref.Draw("sames")
hGvsXhv_peak.Draw("sames")
hGvsX_peak.SetTitle("Pulse Height vs impact point position")
hGvsX_peak.GetYaxis().SetTitle("Pulse Height Peak [adc]")
hGvsX_peak.GetYaxis().SetRangeUser(2500, 3200)
leg = ROOT.TLegend(0.65, 0.7, 0.95, 0.95)
leg.AddEntry(hGvsX_peak, "Gain for 30^{#circ} beam", "l")
leg.AddEntry(hGvsXhv_peak, "Gain for 30^{#circ} beam HVDrift 3100", "l")
leg.AddEntry(hGvsX_ref, "Gain for vertical beam (reference)", "l")
leg.Draw()
cGainVsPosition.cd(2)
gGainNorm.SetTitle("Relative Gain")
gGainNorm.GetXaxis().SetTitle("Impact position X [mm]")
gGainNorm.GetYaxis().SetTitle("Normalized gain")
gGainNorm.Draw("AP")
gGainNorm.Fit('pol1')
ROOT.gPad.Update()
gGainNormhv.Draw("P,same")
gGainNormhv.Fit('pol1')

c1 = ROOT.TCanvas()
hGvsXhv_ratio = hGvsXhv_peak.Clone("hGvsXhv_ratio")
hGvsXhv_ratio.Divide(hGvsX_peak)
hGvsXhv_ratio.Draw()
