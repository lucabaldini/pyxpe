import ROOT
import os, sys
ROOT.gStyle.SetOptFit(1)

fpath = sys.argv[1]
print "open file", fpath
flabel = os.path.basename(fpath).replace(".root", "_Theta")
print flabel

t = ROOT.TChain('tree')
t.Add(fpath)

cut = "fNClusters==1"
#&& ((fImpactY[0]-0.5)**2 + (fImpactX[0]+3.)**2)<(1**2)  "
#&& fPHeight[0]>2000"#" && (fMomX[0]/fMomY[0])>1.5"#"(1)"#
nsigfit = 1.8

print "CUT = ", cut
nTot    = t.GetEntries()
nEvtCut = t.GetEntries(cut)
print "Total Entries %d - after cut %d (eff %.2f %%)" %\
    (nTot, nEvtCut, 100.*nEvtCut/nTot)

print "AvgRate %.3f Hz" % (nTot/(t.GetMaximum('fTimeStamp')-t.GetMinimum('fTimeStamp')))
                           
kPi = ROOT.TMath.Pi()
h0 = ROOT.TH1F("h0", "Theta0", 100, -kPi, kPi)
h1 = ROOT.TH1F("h1", "Theta1", 100, -kPi, kPi)

Phi0 = ROOT.TF1("Phi0", "[0] + [1]*cos(x-[2])*cos(x-[2])", -kPi, kPi) 
Phi1 = ROOT.TF1("Phi1", "[0] + [1]*cos(x-[2])*cos(x-[2])",  -kPi, kPi)

c= ROOT.TCanvas(flabel,flabel,700,700)
c.Divide(2,2)
c.cd(1)
t.Draw("fTheta0[0]>>h0", cut, "HE")
Phi0.SetParameters(h0.GetMaximum(), h0.GetMaximum(), 0.5);
h0.Fit("Phi0", "MR");
A0 = Phi0.GetParameter(0);
B0 = Phi0.GetParameter(1);
ModulationFactor0 = B0/(B0+2*A0)*100.0;
h0.SetTitle("Theta0 modulation factor %.4f %%" % ModulationFactor0)
c.cd(2)
t.Draw("fTheta1[0]>>h1", cut, "HE")
Phi1.SetParameters(h1.GetMaximum(), h1.GetMaximum(), 0.5);
h1.Fit("Phi1", "MR");
A1 = Phi1.GetParameter(0);
B1 = Phi1.GetParameter(1);
ModulationFactor1 = B1/(B1+2*A1)*100.0;
h1.SetTitle("Theta1 modulation factor %.4f %%" % ModulationFactor1)
h1.GetYaxis().SetRangeUser(00, 400)
c.cd(3)
hPH = ROOT.TH1F("hPH", "Pulse Height;[adc]", 200, 0, 10000)
t.Draw("fPHeight[0]>>hPH", cut, "")
m = hPH.GetBinCenter(hPH.GetMaximumBin())
r = hPH.GetRMS()
g = ROOT.TF1("g", "gaus", m-r, m+r)
hPH.Fit("g", "RNQ")
g1 = g.GetParameter(1)
g2 = g.GetParameter(2)
g = ROOT.TF1("g", "gaus", g1 -nsigfit*g2, g1 +nsigfit*g2)
hPH.Fit("g", "RQ")
gPeak  = g.GetParameter(1)
gSigma = g.GetParameter(2)
print gPeak, gSigma, "Eres=", gSigma/gPeak, "FWHM", 2.355*(gSigma/gPeak)
hPH.SetTitle(hPH.GetTitle()+ " {#Delta E/E = %.1f %%}" % (100.*gSigma/gPeak))
c.cd(4)
t.Draw("fCluSize[0]>>hCluSize", cut, "")

ROOT.gPad.Update()
