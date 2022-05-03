import sys
import ROOT
import ROOT.TH1D as TH1D
import ROOT.TH2D as TH2D
import ROOT.TGraph2D as TGraph2D
from ROOT import TCanvas, TColor, TGaxis, TH1F, TPad, TF1
from ROOT import kBlack, kGreen, kRed
from ROOT import TFile, gStyle, TRatioPlot, gPad, TGraphErrors, TGraphAsymmErrors
import math
import array
import glob
import subprocess
import os
import numpy as np
import pandas as pd

from plottingTools import *
from scipy.special import zeta, polygamma, factorial

sys.path.insert(0, './rundec-python/build/lib.linux-x86_64-2.7')
import rundec

debug = True

ROOT.gROOT.SetBatch(True)

pdfList=["NNPDF31_nlo_as_0118"]

massRange = array.array('d',[169.5,170.0,170.5,171.0,171.5,172.0,172.5,173.0,173.5,174.0,174.5,175.0,175.5])
crd = rundec.CRunDec()
#test pole mass
#runningMass=163
#scale=163.2
runningMass=163.2
scale=442.08
nloops = 2
alphaSatScale = crd.AlphasExact(0.118, 91.1876, scale, 5, nloops)
if debug:
    
    print "alpha S at Scale ",scale ," = " , alphaS(0.118,scale)
    for runningMass in np.arange(163.0,164.0,0.1):
    #for runningMass in np.arange(160.0,165.0,0.1):
        print "Scale mu= ", runningMass," , m(mu) = ", runningMass, " | Pole Mass,  m_t = ", poleMass(runningMass,runningMass)
        print "Scale mu= ", runningMass," , m(mu) = ", runningMass, " | Pole Mass,  m_t = ", crd.mMS2mOS(runningMass, None, crd.AlphasExact(0.118, 91.1876, runningMass, 5, nloops), runningMass, 5, nloops)
        print "_____________________________________________________________________________"
        #print "Scale mu= ", (2*runningMass+125.0)/4.," , m(mu) = ", runningMass, " | Pole Mass,  m_t = ", poleMass(runningMass,(2*runningMass+125.0)/4.)





print "m(m) = ", runningMass, " | Pole Mass = ", poleMass(runningMass,runningMass)

print "RUNDEC alphaS  = " , crd.AlphasExact(0.118, 91.1876, scale, 5, nloops), " MSbar Top Mass = ", crd.mMS2mOS(runningMass, None, alphaSatScale, scale, 5, nloops)
print "Running mass at scale ", scale , " = ", crd.mMS2mMS(runningMass, runningMass, scale, 5,nloops)
print "Running mass at scale ", scale , " = ", crd.mOS2mMS(172.5, None, alphaSatScale, scale,5, 2);
#massRange=array.array('d',[161.56])
#runningMass=massRange[0]

sys.exit()
#variableInMadgraph = "pt H"
#variableInMadgraph = "y top"
#pdfList=["CT10_nlo","NNPDF31_nlo_as_0118","NNPDF30_nlo_as_0118","MMHT2014nlo68clas118","ABMP16_5_nlo","ABMP16free_5_nlo"]
pdfList=["MMHT2014nlo68clas118"]

#variableInMadgraphList = ["pt H","y H" ,"inv mass","pt anti","y anti","tt inv","pt top","y top"]


variableInMadgraphList = ["y top"]


massRangePole = array.array('d',[172.5,172.5,172.5])
scaleArrayPole = array.array('d',[470.,470.,470.])
poleVariations = array.array('d',[0.5,0.25,-1]) # -1 is nominal 

#massRangeRunning = array.array('d',[158.54,161.56,155.85])
#scaleArrayRunning = array.array('d',[442.08,224.06,873.41])
massRangeRunning = array.array('d',[161.56,164.95,158.54])
scaleArrayRunning = array.array('d',[224.06,113.72,442.08])

xBins = array.array('d',np.arange(149.75, 175, 0.5))
yBins = array.array('d',np.arange(0., 300., 30.))
yBins = np.append(yBins,array.array('d',[360.,420.,480.,540.,640]))
#yBins=np.append(yBins,np.arange(250., 450., 50.))
#yBins=np.append(yBins,np.arange(550., 750., 100.))
#yBins=np.append(yBins,np.arange(400., 600., 50.))
#yBinsHigh = array.array('d',np.arange(350., 750., 50))
#yBins=np.append(yBins,array.array('d',[380.,480.,580.,740,750]))
#yBins=np.append(yBins,array.array('d',[350.,400.,450.,500.,520.]))

#yBins = array.array('d',np.arange(300, 600, 10))

#xBins=[]
#yBins=[]
crossCheck= False
#Shift entry values in a bin against each other, so they are better readable
shiftValues=True
# Range for ratio shape comparison
minRatio = 0.85
maxRatio = 1.15
if "inv" in variableInMadgraphList[0]:
    minRatio = 0.75
    maxRatio = 1.35

#define colors
poleMassColor=ROOT.kBlue+4
runningMassColor=ROOT.kRed

outputFolder = "RunningMass-DiffXSWithUncertainties_Shifted"

for variableInMadgraph in variableInMadgraphList:
    title, xLabel, yLabel= titleAndLabels(variableInMadgraph)
    #yBins = getYBins(variableInMadgraph)
    yBins = []
    for pdf in pdfList:
        runningMassXS = []
        poleMassHistos = []
        outputFileName = variableInMadgraph+"_"+pdf
        for iterator in range(len(scaleArrayRunning)):
        
            poleMassHisto = make1DHistoList("ttH_NLO_PoleMass_"+pdf+"_M_", "/Events/run_01/MADatNLO.top", variableInMadgraph,[massRangePole[iterator]],scaleArrayPole[iterator],poleVariations[iterator])
            nloXSHisto = make1DHistoList("ttH_runningMass_NLO_"+pdf+"_M_", "/Events/run_01/MADatNLO.top", variableInMadgraph,[massRangeRunning[iterator]],scaleArrayRunning[iterator])

            if len(xBins) and len(yBins):
                poleMassHisto[0]=poleMassHisto[0].Rebin(len(yBins)-1,"Rebinned poleMassHisto",yBins)
                nloXSHisto[0]=nloXSHisto[0].Rebin(len(yBins)-1,"Rebinned nloXSHisto",yBins)
                runningMassXS.append(calcuateDiffXSRunningMass(nloXSHisto, pdf, massRangeRunning[iterator],scaleArrayRunning[iterator], variableInMadgraph,xBins,yBins))
            else:
                runningMassXS.append(calcuateDiffXSRunningMass(nloXSHisto, pdf, massRangeRunning[iterator],scaleArrayRunning[iterator], variableInMadgraph))

            if len(xBins) and len(yBins):
#                poleMassHisto[0]=poleMassHisto[0].Rebin(len(yBins)-1,"Rebinned poleMassHisto",yBins)
#                runningMassXS[iterator]=runningMassXS[iterator].Rebin(len(yBins)-1,"Rebinned runningMassXS",yBins)
                for xbin in range(poleMassHisto[0].GetNbinsX()+2):
                    poleMassHisto[0].SetBinContent(xbin,poleMassHisto[0].GetBinContent(xbin)/poleMassHisto[0].GetBinWidth(xbin))
                    poleMassHisto[0].SetBinError(xbin,poleMassHisto[0].GetBinError(xbin)/poleMassHisto[0].GetBinWidth(xbin))

                    runningMassXS[iterator].SetBinContent(xbin,runningMassXS[iterator].GetBinContent(xbin)/runningMassXS[iterator].GetBinWidth(xbin))
                    runningMassXS[iterator].SetBinError(xbin,runningMassXS[iterator].GetBinError(xbin)/runningMassXS[iterator].GetBinWidth(xbin))


            poleMassHistos.append(poleMassHisto[0])
            #nloXSHisto = make1DHistoList("ttH_runningMass_NLO_"+pdf+"_M", "/Events/run_01/MADatNLO.top", variableInMadgraph,[runningMass])
            #runningMassXS.append(calcuateDiffXSRunningMass(nloXSHisto, pdf, runningMass, scale, variableInMadgraph))


        yBinsPT = getYBins(variableInMadgraph)
        #yBinsPT =[]
        for iterHisto in range(len(poleMassHistos)):
            runningMassXS[iterHisto].Scale(1000.)
            poleMassHistos[iterHisto].Scale(1000.)
            if len(yBinsPT):
                print yBinsPT
                poleMassHistos[iterHisto] = poleMassHistos[iterHisto].Rebin(len(yBinsPT)-1,"Rebinned nloXSHisto",yBinsPT)
                runningMassXS[iterHisto]  = runningMassXS[iterHisto].Rebin(len(yBinsPT)-1,"Rebinned Running Mass",yBinsPT)
                poleMassHistos[iterHisto].Print("all")
                for iBin in range(poleMassHistos[iterHisto].GetNbinsX()+1):
                    poleMassHistos[iterHisto].SetBinContent(iBin+1, poleMassHistos[iterHisto].GetBinContent(iBin+1)/poleMassHistos[iterHisto].GetBinWidth(iBin+1))
                    poleMassHistos[iterHisto].SetBinError(iBin+1, poleMassHistos[iterHisto].GetBinError(iBin+1)/poleMassHistos[iterHisto].GetBinWidth(iBin+1))
                    
                    runningMassXS[iterHisto].SetBinContent(iBin+1, runningMassXS[iterHisto].GetBinContent(iBin+1)/runningMassXS[iterHisto].GetBinWidth(iBin+1))
                    runningMassXS[iterHisto].SetBinError(iBin+1, runningMassXS[iterHisto].GetBinError(iBin+1)/runningMassXS[iterHisto].GetBinWidth(iBin+1))


        #poleMassHistos[2].Print("all")
        #sys.exit()
        gStyle.SetPadLeftMargin(0.16)
        canvas, pad1, pad2, pad3 = createCanvas3Pads()
        gStyle.SetPadLeftMargin(0.16)
        gStyle.SetTitleFontSize(0.07)
        canvas.cd()

        pad1.cd()

        #pad1.SetLogy()
        #canvas = TCanvas("","")
        if crossCheck:
            poleMassHistos[0].Draw()
            poleMassHistos[1].SetLineColor(ROOT.kGreen)
            poleMassHistos[1].Draw("same")
            poleMassHistos[2].SetLineColor(ROOT.kBlack)
            poleMassHistos[2].Draw("same")

            runningMassXS[0].SetLineColor(ROOT.kViolet)
            runningMassXS[0].Draw("same")
            runningMassXS[1].SetLineColor(ROOT.kViolet+1)
            runningMassXS[1].Draw("same")
            runningMassXS[2].SetLineColor(ROOT.kViolet+2)
            runningMassXS[2].Draw("same")
            canvas.Update()
            canvas.Print(outputFolder+"/test_SystematicBand.pdf")
            sys.exit()

        
        
        graphPoleWithUncertainties = makeDiffXSGraphError(poleMassHistos)
        graphRunningWithUncertainties = makeDiffXSGraphError(runningMassXS,False,False)
        graphPoleWithUncertainties.Print("all")
        graphPoleWithUncertainties.SetTitle(title)
        graphRunningWithUncertainties.SetTitle(title)

        setGraphStyle(graphPoleWithUncertainties, 1, poleMassColor, 1, -1, poleMassColor, -1,1154,poleMassColor)
        setGraphStyle(graphRunningWithUncertainties, 1, runningMassColor, 1, -1, runningMassColor, -1,3345,runningMassColor)
        yAxis = graphPoleWithUncertainties.GetYaxis()
        yAxis.SetTitle(yLabel)
        yAxis.SetNdivisions(505)
        yAxis.SetTitleSize(40)
        yAxis.SetTitleFont(43)
        yAxis.SetTitleOffset(1.4)
        yAxis.SetLabelFont(43)
        yAxis.SetLabelSize(25)
        graphPoleWithUncertainties.GetHistogram().Draw()
        graphPoleWithUncertainties.Draw("E3")
        graphPoleWithUncertainties.Draw("P")
        graphRunningWithUncertainties.Draw("E3")
        graphRunningWithUncertainties.Draw("P")

        leg = ROOT.TLegend(.65,.65,.9,0.87)
        leg.SetBorderSize(0)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.06)
        #leg.AddEntry(h2,"MMHT NNLO","L")
        #leg.AddEntry(h1,"MMHT NLO","L")
        leg.AddEntry(graphPoleWithUncertainties,"Pole Mass","L")
        leg.AddEntry(graphRunningWithUncertainties,"#bar{MS} Mass","L")
        leg.Draw("same")
#        DrawErrorBand(graphPoleWithUncertainties)
        #graphPoleWithUncertainties.Draw("same")
        
#        DrawErrorBand(graphRunningWithUncertainties)
        #graphRunningWithUncertainties.Draw("same")


        #h3 = createRatio( runningMassXS[0], poleMassHistos[0], variableInMadgraph,"#frac{#sigma(#mu)}{#sigma_{m_{Top}}}")
        #h3.Draw()
#        drawRatioPad(canvas, 0.8+0.01, 1.2-0.01, "#frac{Running Mass}{Pole Mass}")
#        gPad.SetGrid(0, 1)
#        gPad.RedrawAxis("g")
        #print "----------------------------NOMINAL-----------------------------------------"
        #poleMassHistos[0].Print("all")
        #print "----------------------------DOWN VARIATION-----------------------------------------"
        #poleMassHistos[1].Print("all")
        #print "----------------------------UP VARIATION-----------------------------------------"
        #poleMassHistos[2].Print("all")
        #sys.exit()

        

        graphPoleRatio = makeDiffXSGraphError(poleMassHistos, False)
        graphRunningRatio = makeDiffXSGraphError(runningMassXS, False)

        print "============= Printing CROSS SECTIONS ============="
        printCrossSections(poleMassHistos,runningMassXS)
        print "==================================================="

        graphPoleRatio.SetTitle("")
        graphRunningRatio.SetTitle("")

        ratioGraphPoleErrorRelative = makeRatioGraph(graphPoleRatio,graphPoleRatio,1)
        ratioGraphScaleErrorRelative = makeRatioGraph(graphRunningRatio,graphPoleRatio,1,shiftValues)

        ratioGraphPoleErrorRelative.SetTitle("")
        ratioGraphScaleErrorRelative.SetTitle("")

        ratioGraphPoleError = makeRatioGraph(graphPoleRatio,graphPoleRatio,1)
        ratioGraphScaleError = makeRatioGraph(graphRunningRatio,graphRunningRatio,1,shiftValues)

        ratioGraphPoleError.SetTitle("")
        ratioGraphScaleError.SetTitle("")

        setGraphStyle(ratioGraphPoleErrorRelative, graphPoleWithUncertainties.GetLineStyle(), 
            graphPoleWithUncertainties.GetLineColor(), graphPoleWithUncertainties.GetLineWidth(), 
            graphPoleWithUncertainties.GetMarkerStyle(), graphPoleWithUncertainties.GetMarkerColor(),
            graphPoleWithUncertainties.GetMarkerSize(), graphPoleWithUncertainties.GetFillStyle(),
            graphPoleWithUncertainties.GetFillColor())

        setGraphStyle(ratioGraphScaleErrorRelative, graphRunningWithUncertainties.GetLineStyle(), 
            graphRunningWithUncertainties.GetLineColor(), graphRunningWithUncertainties.GetLineWidth(), 
            graphRunningWithUncertainties.GetMarkerStyle(), graphRunningWithUncertainties.GetMarkerColor(),
            graphRunningWithUncertainties.GetMarkerSize(), graphRunningWithUncertainties.GetFillStyle(),
            graphRunningWithUncertainties.GetFillColor())

        
        setGraphStyle(ratioGraphPoleError, graphPoleWithUncertainties.GetLineStyle(), 
            graphPoleWithUncertainties.GetLineColor(), graphPoleWithUncertainties.GetLineWidth(), 
            graphPoleWithUncertainties.GetMarkerStyle(), graphPoleWithUncertainties.GetMarkerColor(),
            graphPoleWithUncertainties.GetMarkerSize(), graphPoleWithUncertainties.GetFillStyle(),
            graphPoleWithUncertainties.GetFillColor())

        setGraphStyle(ratioGraphScaleError, graphRunningWithUncertainties.GetLineStyle(), 
            graphRunningWithUncertainties.GetLineColor(), graphRunningWithUncertainties.GetLineWidth(), 
            graphRunningWithUncertainties.GetMarkerStyle(), graphRunningWithUncertainties.GetMarkerColor(),
            graphRunningWithUncertainties.GetMarkerSize(), graphRunningWithUncertainties.GetFillStyle(),
            graphRunningWithUncertainties.GetFillColor())


        pad2.cd()
        ratioGraphPoleErrorRelative = setHistoProperties(ratioGraphPoleErrorRelative, "", "#frac{#sigma(m(#mu))}{#sigma(m_{Pole})}", minRatio, maxRatio)
        ratioGraphPoleErrorRelative.GetHistogram().Draw()
        ratioGraphPoleErrorRelative.Draw("LP")
        ratioGraphScaleErrorRelative.Draw("LP")

        pad3.cd()
        ratioGraphPoleError = setHistoProperties(ratioGraphPoleError, xLabel, "#frac{#frac{#Delta d#sigma}{dX}}{#frac{d#sigma}{dX}}")
        ratioGraphPoleError.GetHistogram().Draw()

        ratioGraphPoleError.Draw("C3")
        ratioGraphPoleError.Draw("P")
        ratioGraphScaleError.Draw("C3")
        ratioGraphScaleError.Draw("P")

#        DrawErrorBand(ratioGraphPoleError)
#        DrawErrorBand(ratioGraphScaleError)
        #ratioGraphScaleError.Draw("same")
        

        gPad.RedrawAxis()
        gPad.Update()
        gPad.Modified()

        canvas.Update()

        canvas.SaveAs(outputFolder+"/"+variableInMadgraph.replace(" ","")+"_SystematicBand.pdf")
        sys.exit()
        
        # Draw grid
        
        

        

        # Draw horizontal line at y=1
        axisHisto = getPadAxisHisto(canvas)
        xmin = axisHisto.GetXaxis().GetXmin()
        xmax = axisHisto.GetXaxis().GetXmax()
        height = TString("")
        height += 1
        line1 = TF1("line1", height, xmin, xmax)
        line1.SetLineStyle(1)
        line1.SetLineWidth(1)
        line1.SetLineColor(kBlack)
        line1.Draw("l,same")

        gPad.RedrawAxis()
        gPad.Update()
        gPad.Modified()

        #
        

        c, pad1, pad2 = createCanvasPads()        

        pad1.cd()

        histoWithUncertainties.SetTitle("M(mu) = "+str(runningMass)+" Scale"+str(scale)+" , "+pdf)
        histoWithUncertainties.SetFillColor(38)
        histoWithUncertainties.SetFillStyle(3001)
        #histoWithUncertainties.SetLineColor(kBlack)
        histoWithUncertainties.Draw()
        #runningMassXS.SetLineColor(kRed-1)
        #histoWithUncertainties.Draw("same")

        leg = ROOT.TLegend(.65,.65,.9,0.87)
        leg.SetBorderSize(0)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.04)
        #leg.AddEntry(h2,"MMHT NNLO","L")
        #leg.AddEntry(h1,"MMHT NLO","L")
        leg.AddEntry(histoWithUncertainties,"XS, Scale Variations","L")
       # leg.AddEntry(runningMassXS,"RunningMass","L")
        leg.Draw("same")
        # to avoid clipping the bottom zero, redraw a small axis
        #h1.GetYaxis().SetLabelSize(0.0)
        axis = TGaxis(-5, 20, -5, 220, 20, 220, 510, "")
        axis.SetLabelFont(43)
        axis.SetLabelSize(15)
        axis.Draw()

        c.SaveAs(outputFolder+"/"+outputFileName+".pdf")


sys.exit()

OutputFileName = "Massdependence/"
#variableInMadgraph = "pt H"
variableInMadgraph = "y H"
#variableInMadgraph = "inv mass"
#variableInMadgraph = "pt anti"
#variableInMadgraph = "y anti"
#variableInMadgraph = "tt inv"
#variableInMadgraph = "pt top"
#variableInMadgraph = "y top"
#VariableName= "ptH_*.dat"

histoTitle = "p_{T} of Higgs Boson, MG5aMC@NLO"
#histoTitle = "p_{T} of Higgs Boson, NNLL"
#histoTitle = "M_{ttH}, MG5aMC@NLO"
#YTitle = "#frac{d#sigma}{dy}"
YTitle = ""
XTitle = "Top Mass [GeV]"
ZTitle = "#frac{d#sigma}{dX}"
#variableInMadgraph = "inv mass"
VariableName= "mass_*.dat"
#histoTitle = "Rapidity of Higgs Boson, NNLL"
massRange = array.array('d',[169.5,170.0,170.5,171.0,171.5,172.0,172.5,173.0,173.5,174.0,174.5,175.0,175.5])

#binCenters,binEntries,binUncertainties = makeHistogramFromMG5File("ttH_BORN_NNPDFLO30_M175.5/Events/run_01_LO/MADatNLO.top", variableInMadgraph)
#histo = fillAndReturnHistogram(histoTitle,binCenters,binEntries,binUncertainties)
#binCenters,binEntries,binUncertainties = makeHistogramFromMG5File("ttH_BORN_NNPDFLO30_M175.5/Events/run_01_LO/MADatNLO.top", variableInMadgraph)
#histo1 = fillAndReturnHistogram(histoTitle,binCenters,binEntries,binUncertainties)
#yAxisLowEdges = array.array('d',binCenters)
#print histo1.GetNbinsX(), make2DHistfromListOf1DHisto([histo1],massRange)
##################################3
#1dHistoList=make1DHistoList("ttH_BORN_NNPDFLO30_M", "/Events/run_01_LO/MADatNLO.top", variableInMadgraph,massRange)

#binLowEdges=makeBinLowEdges(massRange)
#allhistos=makeListOfMassDependenceHistos(make1DHistoList("ttH_BORN_CT14lo_M", "/Events/run_01_LO/MADatNLO.top", variableInMadgraph,massRange),binLowEdges)
#allhistos=fitListofHistos(allhistos)


c = TCanvas("c", "canvas", 800, 800)
for histoID, histo in enumerate(allhistos):
    histo.Draw()
    if histoID == 0:
        c.Print("CT14lo_Massdependence_ptH.pdf(","pdf")
    elif histoID == len(allhistos)-1:
        c.Print("CT14lo_Massdependence_ptH.pdf)","pdf")
    else:
        c.Print("CT14lo_Massdependence_ptH.pdf","pdf")

sys.exit()
##################################3
histo2D, histoTest=make2DHistfromListOf1DHisto(make1DHistoList("ttH_BORN_NNPDFLO30_M", "/Events/run_01_LO/MADatNLO.top", variableInMadgraph,massRange),massRange)
gx,gy=calculateGradient(histoTest)
c1=TCanvas()
c1.cd()
histogram1D = histoTest.ProjectionX("Projection in Y",0,-1)
histogram1D.Draw()
input()
sys.exit()
#histo2D.GetXaxis().SetTitle(XTitle)
#histo2D.GetYaxis().SetTitle(YTitle)
#histo2D.GetZaxis().SetTitle(ZTitle)

#histo2D.Draw("colz")
if ' ' in variableInMadgraph:
    variableInMadgraph=variableInMadgraph.replace('', '_')
c1.SaveAs(OutputFileName+variableInMadgraph+".root")

print gx, len(gy)


sys.exit()

# draw everything
pad1.cd()
#h1.Draw()
#h2.Draw("same")
histo.SetStats(0)
histo.SetLineColor(kBlack)
histo.Draw()
histo1.Draw("same")

leg = ROOT.TLegend(.65,.65,.9,0.87)
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetTextFont(42)
leg.SetTextSize(0.04)
#leg.AddEntry(h2,"MMHT NNLO","L")
#leg.AddEntry(h1,"MMHT NLO","L")
leg.AddEntry(histo,"MMHT NNLO","L")
leg.AddEntry(histo1,"MMHT NLO","L")
leg.Draw("same")
# to avoid clipping the bottom zero, redraw a small axis
#h1.GetYaxis().SetLabelSize(0.0)
axis = TGaxis(-5, 20, -5, 220, 20, 220, 510, "")
axis.SetLabelFont(43)
axis.SetLabelSize(15)
axis.Draw()
pad2.cd()
h3.Draw("ep")

c.SaveAs(OutputFileName+".pdf")

