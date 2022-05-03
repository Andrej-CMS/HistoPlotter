import sys
import ROOT
import ROOT.TH1D as TH1D
import ROOT.TH2D as TH2D
import ROOT.TGraph2D as TGraph2D
from ROOT import TCanvas, TColor, TGaxis, TH1F, TPad, TF1, TString
from ROOT import kBlack, kGreen, kRed
from ROOT import TFile, gStyle, TRatioPlot, gPad, TGraphErrors, TGraphAsymmErrors
import math
import array
import glob
import subprocess
import os
import numpy as np
import pandas as pd

from plottingUtils import *

ROOT.gROOT.SetBatch(True)


outputPath = "ComparisonPlots"
variable = "jet_pt"
normalize = True


# Set properties of plots based on dictionary of variables
variableFixedOrder = getListOfVariableProperties(variable)[0]
xTitle = getListOfVariableProperties(variable)[1]
yTitle = getListOfVariableProperties(variable)[2]
numberOfBins = getListOfVariableProperties(variable)[3]
xBinsLow = getListOfVariableProperties(variable)[4][0]
xBinsHigh = getListOfVariableProperties(variable)[4][1]
xBins = array.array('d',np.linspace(xBinsLow, xBinsHigh, numberOfBins+1))
print xBins

ratioMin = getListOfVariableProperties(variable)[5][0]
ratioMax = getListOfVariableProperties(variable)[5][1]
histoTitle = getListOfVariableProperties(variable)[6]
totalXS_FixedOrder_Name = getListOfVariableProperties(variable)[7]
treeName = getListOfVariableProperties(variable)[8]


yTitleRatio = "ratio"
yTitleAddition = ""
if normalize:
    yTitleAddition = "#frac{1}{#sigma}"

histoList = []
histo     = ROOT.TH1D()

legendNameList = []
# run over all input files and plot them into one histogram
for argNumber in range(len(sys.argv)-1):
    
    
    fileName = sys.argv[argNumber+1]
    print " Reading from file", fileName 
    inFile = ROOT.TFile.Open( fileName ," READ ")

    legendNameList.append(fileName[:-5])

    tree  = inFile.Get (treeName)

    histo = ROOT.TH1D( fileName, histoTitle+str(argNumber+1) , numberOfBins , xBins )
    # SetDirectory(0) to keep the histogram instead of being deleted after file close
    histo.SetDirectory(0)
    histo.SetMinimum(0.001)

    sumOfWeights = 0.
    factor_GeV = 1.
    isRhoVariable = True


    if not "_rho_" in variable:
        factor_GeV = 1000.
        isRhoVariable = False

    for entryNum in range (0 , tree.GetEntries()):

        tree.GetEntry( entryNum )

        if variable == "jet_HT":
            eventValue      = getattr( tree ,"jet_pt")
        else:
            eventValue      = getattr( tree ,variable)

        eventWeight     = getattr( tree , "weight_mc")
        sumOfWeights    += eventWeight
        jet_HT = 0.

        if isRhoVariable:
            jetPt = getattr( tree ,"MC_j_afterFSR_pt")
            if jetPt < 30.:
                continue
        if variable == "jet_pt":
            for value in eventValue:
                histo.Fill(value/factor_GeV, eventWeight)
        elif variable == "jet_HT":
            for value in eventValue:
                jet_HT+=value
            histo.Fill(jet_HT/factor_GeV, eventWeight)
        else:
            histo.Fill(eventValue/factor_GeV, eventWeight)

    crossSectionFactor = sumOfWeights/tree.GetEntries()

    if normalize:
        histo.Scale(1./histo.Integral())
    else:
        #need to normalize to expected events at 1 inverse femtobarn
        histo.Scale(1./crossSectionFactor)

    print "crossSectionFactor, histo   ", crossSectionFactor, " ", histo

    histoList.append(histo)


print histoList

# set style for canvas
gStyle.SetPadLeftMargin(0.16)
gStyle.SetTitleFontSize(0.05)

canvas, pad1, pad2 = createCanvasPads()
canvas.cd()

#create legend
legend = ROOT.TLegend(.65,.65,.82,0.87)
legend.SetBorderSize(0)
legend.SetFillColor(0)
legend.SetFillStyle(0)
legend.SetTextFont(42)
legend.SetTextSize(0.05)


listOfLegendNames = []
# loop over all histos in histolist and draw them into the canvas and ratios
for histoNumber in range(len(histoList)):
    
    histoList[histoNumber].GetYaxis().SetTitle(yTitleAddition+yTitle)
    histoList[histoNumber].GetYaxis().SetTitleSize(0.05)
    histoList[histoNumber].GetXaxis().SetTitle(xTitle)
    histoList[histoNumber].SetMaximum(1.2*histoList[argNumber].GetMaximum())
    histoList[histoNumber].SetTitle(histoTitle)

    histoList[histoNumber].SetStats(0)
    histoList[histoNumber].SetLineWidth(2)
    histoList[histoNumber].SetLineColor(histoNumber+1)

    legendName_tmp  = legendNameList[histoNumber].split("/")
    legendName      = legendName_tmp[len(legendName_tmp)-1]

    legend.AddEntry(histoList[histoNumber], getLegendNames(legendName),"L")
    listOfLegendNames.append(legendName)
    pad1.cd()
    if histoNumber == 0:
        histoList[histoNumber].Draw()
    else:
        histoList[histoNumber].Draw("same")

    if histoNumber == len(histoList)-1:
        legend.Draw("same")

pad2.cd()
ratioList=[]
for histoNumber in range(len(histoList)):
    ratio = createRatio(histoList[histoNumber], histoList[histoNumber], "Top quark p_{T}", "#frac{Pythia 8}{Herwig 7}")
    ratio.SetDirectory(0)
    ratio = setHistoProperties(ratio, xTitle, yTitleRatio , ratioMin, ratioMax)

    ratioList.append(ratio)


for ratioNumber in range(len(ratioList)):
    #ratioList[ratioNumber].Print("all")
    printUncertaintiesToLatex(ratioList[ratioNumber],ratioList[0],listOfLegendNames[ratioNumber])
    if ratioNumber == 0:
        ratioList[ratioNumber].Draw()
    else:
        ratioList[ratioNumber].Draw("same")


if not os.path.exists(outputPath):
  
  # Create a new directory because it does not exist 
  os.makedirs(outputPath)
  print("The new directory is created!")

canvas.Update()
string_normalize = ""
if normalize:
    string_normalize = "_normalized"

canvas.SaveAs(outputPath+"/"+variable+string_normalize+".pdf")


# pythiaColor = ROOT.kViolet-1
# herwigColor = ROOT.kGreen+4
# mg5Color = ROOT.kBlue
# fixedOrderColor = ROOT.kBlack





# sumOfWeights = 0.
# factor_GeV = 1.
# isRhoVariable = True
# if not "_rho_" in variable:
#     factor_GeV = 1000.
#     isRhoVariable = False
# for entryNum in range (0 , treePythia.GetEntries()):
#         treePythia.GetEntry( entryNum )
#         topPt_afterFSR_Pythia = getattr( treePythia ,variable)
#         eventWeight           = getattr( treePythia , "weight_mc")
#         sumOfWeights += eventWeight

#         if isRhoVariable:
#             jetPt = getattr( treePythia ,"MC_j_afterFSR_pt")
#             if jetPt < 30.:
#                 continue

#         if topPt_afterFSR_Pythia > 0.:
#             histo_topPt_afterFSR_Pythia.Fill(topPt_afterFSR_Pythia/factor_GeV, eventWeight)

# crossSectionFactorPythia8 = sumOfWeights/treePythia.GetEntries()
# sumOfWeights = 0.
# for entryNum in range (0 , treeHerwig.GetEntries()):
#         treeHerwig.GetEntry( entryNum )
#         topPt_afterFSR_Herwig = getattr( treeHerwig ,variable)
#         eventWeight           = getattr( treeHerwig , "weight_mc")
#         sumOfWeights += eventWeight

#         if isRhoVariable:
#             jetPt = getattr( treeHerwig ,"MC_j_afterFSR_pt")
#             if jetPt < 30.:
#                 continue

#         if topPt_afterFSR_Herwig > 0.:
#             histo_topPt_afterFSR_Herwig.Fill(topPt_afterFSR_Herwig/factor_GeV, eventWeight)

# crossSectionFactorHerwig = sumOfWeights/treeHerwig.GetEntries()
# sumOfWeights = 0.
# for entryNum in range (0 , treeMG5.GetEntries()):
#         treeMG5.GetEntry( entryNum )
#         topPt_afterFSR_MG5 = getattr( treeMG5 ,variable)
#         eventWeight        = getattr( treeMG5 , "weight_mc")
#         sumOfWeights += eventWeight
        
#         if isRhoVariable:
#             jetPt = getattr( treeMG5 ,"MC_j_afterFSR_pt")
#             if jetPt < 30.:
#                 continue

#         if topPt_afterFSR_MG5 > 0.:
#             histo_topPt_afterFSR_MG5.Fill(topPt_afterFSR_MG5/factor_GeV, eventWeight)

# crossSectionFactorMG5 = sumOfWeights/treeMG5.GetEntries()

# if normalize:
#     histo_topPt_afterFSR_Pythia.Scale(1./histo_topPt_afterFSR_Pythia.Integral())
#     histo_topPt_afterFSR_Herwig.Scale(1./histo_topPt_afterFSR_Herwig.Integral())
#     histo_topPt_afterFSR_FixedOrder.Scale(1./histo_topPt_afterFSR_FixedOrder.Integral())
#     histo_topPt_afterFSR_MG5.Scale(1./histo_topPt_afterFSR_MG5.Integral())
# else:
#     #need to normalize to expected events at 1 inverse femtobarn
#     histo_topPt_afterFSR_Pythia.Scale(1./crossSectionFactorPythia8)
#     histo_topPt_afterFSR_Herwig.Scale(1./crossSectionFactorHerwig)
#     # print "Norm Fixed Order = ", inFileFixedOrder.Get(totalXS_FixedOrder_Name).GetBinContent(1)
#     histo_topPt_afterFSR_FixedOrder.Scale(inFileFixedOrder.Get(totalXS_FixedOrder_Name).GetBinContent(1))
#     histo_topPt_afterFSR_MG5.Scale(1./crossSectionFactorMG5)


# # canvas, pad1, pad2 = createCanvasPads()

# # canvas = TCanvas("c", "canvas", 1024, 1024)
# print "Integral FO = ", histo_topPt_afterFSR_FixedOrder.Integral()
# print "Integral H7 = ", histo_topPt_afterFSR_Herwig.Integral()
# # canvas.cd()

# # histo_topPt_afterFSR_Pythia.Draw()



# # histo_topPt_afterFSR_Herwig.Draw("same")

# histo_topPt_afterFSR_Pythia.GetYaxis().SetTitle(yTitleAddition+yTitle)
# histo_topPt_afterFSR_Pythia.GetYaxis().SetTitleSize(0.04)
# histo_topPt_afterFSR_Pythia.GetXaxis().SetTitle(xTitle)
# histo_topPt_afterFSR_Pythia.SetMaximum(1.2*histo_topPt_afterFSR_Pythia.GetMaximum())
# histo_topPt_afterFSR_Pythia.SetTitle(histoTitle)

# histo_topPt_afterFSR_Pythia.SetStats(0)

# histo_topPt_afterFSR_FixedOrder.GetYaxis().SetTitle(yTitleAddition+yTitle)
# histo_topPt_afterFSR_FixedOrder.GetYaxis().SetTitleSize(0.06)
# histo_topPt_afterFSR_FixedOrder.GetXaxis().SetTitle(xTitle)

# histo_topPt_afterFSR_FixedOrder.SetStats(0)
# histo_topPt_afterFSR_FixedOrder.SetMinimum(0.001) #avoid showing 0 on axis

# gStyle.SetPadLeftMargin(0.16)
# gStyle.SetTitleFontSize(0.05)

# canvas, pad1, pad2 = createCanvasPads()

# canvas.cd()

# pad1.cd()

# histo_topPt_afterFSR_Pythia.SetLineColor(pythiaColor)
# histo_topPt_afterFSR_Herwig.SetLineColor(herwigColor)
# histo_topPt_afterFSR_FixedOrder.SetLineColor(fixedOrderColor)
# histo_topPt_afterFSR_MG5.SetLineColor(mg5Color)

# histo_topPt_afterFSR_Pythia.SetLineWidth(2)
# histo_topPt_afterFSR_Herwig.SetLineWidth(2)
# histo_topPt_afterFSR_FixedOrder.SetLineWidth(2)
# histo_topPt_afterFSR_MG5.SetLineWidth(2)


# histo_topPt_afterFSR_FixedOrder.Draw()
# histo_topPt_afterFSR_Pythia.Draw("same")
# histo_topPt_afterFSR_Herwig.Draw("same")
# histo_topPt_afterFSR_MG5.Draw("same")

# leg = ROOT.TLegend(.65,.65,.82,0.87)
# leg.SetBorderSize(0)
# leg.SetFillColor(0)
# leg.SetFillStyle(0)
# leg.SetTextFont(42)
# leg.SetTextSize(0.05)
# leg.AddEntry(histo_topPt_afterFSR_FixedOrder,"t#bar{t}+jet FO","L")
# leg.AddEntry(histo_topPt_afterFSR_Pythia,"t#bar{t}+Pwhg+P8","L")
# leg.AddEntry(histo_topPt_afterFSR_Herwig,"t#bar{t}+Pwhg+H7","L")
# leg.AddEntry(histo_topPt_afterFSR_MG5,"t#bar{t} MG5+P8","L")
# leg.Draw("same")

# pad2.cd()

# ratioHerwig = createRatio(histo_topPt_afterFSR_Herwig, histo_topPt_afterFSR_FixedOrder, "Top quark p_{T}", "#frac{Pythia 8}{Herwig 7}")
# ratioPythia = createRatio(histo_topPt_afterFSR_Pythia, histo_topPt_afterFSR_FixedOrder, "Top quark p_{T}", "#frac{Pythia 8}{Herwig 7}")
# ratioMG5    = createRatio(histo_topPt_afterFSR_MG5, histo_topPt_afterFSR_FixedOrder, "Top quark p_{T}", "#frac{Pythia 8}{Herwig 7}")



# ratioHerwig = setHistoProperties(ratioHerwig, xTitle, yTitleRatio , ratioMin, ratioMax)
# ratioPythia = setHistoProperties(ratioPythia, xTitle, yTitleRatio , ratioMin, ratioMax)
# ratioMG5 = setHistoProperties(ratioMG5, xTitle, yTitleRatio , ratioMin, ratioMax)

# ratioHerwig.GetYaxis().SetRangeUser(ratioMin,ratioMax)
# ratioPythia.GetYaxis().SetRangeUser(ratioMin,ratioMax)
# ratioMG5.GetYaxis().SetRangeUser(ratioMin,ratioMax)

# ratioHerwig.Draw("")
# ratioPythia.Draw("same")
# ratioMG5.Draw("same")
# # draw a horizontal line on the pad
# xmin = ratioHerwig.GetXaxis().GetXmin()
# xmax = ratioHerwig.GetXaxis().GetXmax()
# # height = TString("") 
# # height += 1
# f = TF1("f", "1", xmin, xmax)
# f.SetLineStyle(1)
# f.SetLineWidth(1)
# f.SetLineColor(kBlack)
# f.Draw("L same")



# # axis = TGaxis(-5, 20, -5, 220, 20, 220, 510, "")
# # axis.SetLabelFont(43)
# # axis.SetLabelSize(15)
# # axis.Draw()

# if not os.path.exists(outputPath):
  
#   # Create a new directory because it does not exist 
#   os.makedirs(outputPath)
#   print("The new directory is created!")

# canvas.Update()
# string_normalize = ""
# if normalize:
#     string_normalize = "normalized"

# canvas.SaveAs(outputPath+"/"+variable+"_PythiaHerwig_"+string_normalize+".pdf")

