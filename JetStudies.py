#import ROOT
#import os
#os.environ["X509_USER_PROXY"] = "/eos/home-t/ttedesch/SWAN_projects/RDataFrame_porting/proxy"
#!xrdcp root://cms-xrd-global.cern.ch//store/mc/RunIIFall17NanoAODv7/VBS_SSWW_LL_polarization_TuneCP5_13TeV-madgraph-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/100000/AA101ED1-3427-E44D-9D6E-82569A5A589C.root .
#!xrdcp root://cms-xrd-global.cern.ch//store/mc/RunIIFall17NanoAODv7/QCD_HT1000to1500_BGenFilter_TuneCP5_13TeV-madgraph-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/100000/14659586-1A5E-8E40-921D-05D85DBCA1DC.root .

import ROOT
import os

ROOT.ROOT.EnableImplicitMT()

ROOT.gInterpreter.Declare('#include "postselection.h"')

df = ROOT.RDataFrame("Events", "AA101ED1-3427-E44D-9D6E-82569A5A589C.root ")
df_atleast2Jets = df.Filter("nJet>2", "At least two jets")
df_GoodJets = df_atleast2Jets.Define("GoodJets_idx", "GoodJets(Jet_jetId, Jet_eta, Jet_pt, Jet_puId)")
df_atleast2GoodJets = df_GoodJets.Filter("atleast2GoodJets(GoodJets_idx)", "At least two good jets")
df_VBSjets = df_atleast2GoodJets.Define("VBSJet_idx", "SelectVBSJets_invmass(Jet_pt, Jet_eta, Jet_phi, Jet_puId, GoodJets_idx)")
df_2VBSjets = df_VBSjets.Filter("VBSJet_idx[0] != VBSJet_idx[1]", "2 VBS jets")


df_VBSleadingjets = df_2VBSjets.Define("Leadingjet_pt", "GetLeading(Jet_pt, VBSJet_idx)")\
                              .Define("Leadingjet_eta", "GetLeading(Jet_eta, VBSJet_idx)")\
                              .Define("Leadingjet_phi", "GetLeading(Jet_phi, VBSJet_idx)")\
                              .Define("Leadingjet_mass", "GetLeading(Jet_mass, VBSJet_idx)")\
                              .Define("Leadingjet_area", "GetLeading(Jet_area, VBSJet_idx)")\
                              .Define("Leadingjet_nConstituents", "GetLeading(Jet_nConstituents, VBSJet_idx)")\
                              .Define("Leadingjet_nElectrons", "GetLeading(Jet_nElectrons, VBSJet_idx)")\
                              .Define("Leadingjet_nMuons", "GetLeading(Jet_nMuons, VBSJet_idx)")\
                              .Define("Leadingjet_neEmEF", "GetLeading(Jet_neEmEF, VBSJet_idx)")\
                              .Define("Leadingjet_neHEF", "GetLeading(Jet_neHEF, VBSJet_idx)")\
                              .Define("Leadingjet_qgl", "GetLeading(Jet_qgl, VBSJet_idx)")\
                              .Define("SubLeadingjet_pt", "GetSubLeading(Jet_pt, VBSJet_idx)")\
                              .Define("SubLeadingjet_eta", "GetSubLeading(Jet_eta, VBSJet_idx)")\
                              .Define("SubLeadingjet_phi", "GetSubLeading(Jet_phi, VBSJet_idx)")\
                              .Define("SubLeadingjet_mass", "GetSubLeading(Jet_mass, VBSJet_idx)")\
                              .Define("SubLeadingjet_area", "GetSubLeading(Jet_area, VBSJet_idx)")\
                              .Define("SubLeadingjet_nConstituents", "GetSubLeading(Jet_nConstituents, VBSJet_idx)")\
                              .Define("SubLeadingjet_nElectrons", "GetSubLeading(Jet_nElectrons, VBSJet_idx)")\
                              .Define("SubLeadingjet_nMuons", "GetSubLeading(Jet_nMuons, VBSJet_idx)")\
                              .Define("SubLeadingjet_neEmEF", "GetSubLeading(Jet_neEmEF, VBSJet_idx)")\
                              .Define("SubLeadingjet_neHEF", "GetSubLeading(Jet_neHEF, VBSJet_idx)")\
                              .Define("SubLeadingjet_qgl", "GetSubLeading(Jet_qgl, VBSJet_idx)")\
                              .Define("mjj", "GetInvMass(Jet_pt, Jet_eta, Jet_phi, Jet_mass, VBSJet_idx)")\

h_VBS = {}

#h = df_atleast2GoodJets.Histo1D("nJet")
h_VBS["Leadingjet_pt"] = df_VBSleadingjets.Histo1D(("Leadingjet_pt", "Leadingjet_pt", 100, 0, 1000),"Leadingjet_pt")
h_VBS["Leadingjet_eta"] = df_VBSleadingjets.Histo1D(("Leadingjet_eta", "Leadingjet_eta", 100, -5, 5),"Leadingjet_eta")
h_VBS["Leadingjet_phi"] = df_VBSleadingjets.Histo1D(("Leadingjet_phi", "Leadingjet_phi", 100, -3.5, 3.5),"Leadingjet_phi")
h_VBS["Leadingjet_mass"] = df_VBSleadingjets.Histo1D(("Leadingjet_mass", "Leadingjet_mass", 100, 0, 200),"Leadingjet_mass")
h_VBS["Leadingjet_area"] = df_VBSleadingjets.Histo1D(("Leadingjet_area", "Leadingjet_area", 100, 0, 0.8),"Leadingjet_area")
h_VBS["Leadingjet_nConstituents"] = df_VBSleadingjets.Histo1D(("Leadingjet_nConstituents", "Leadingjet_nConstituents", 100, 0, 150),"Leadingjet_nConstituents")
h_VBS["Leadingjet_nElectrons"] = df_VBSleadingjets.Histo1D(("Leadingjet_nElectrons", "Leadingjet_nElectrons", 10, 0, 10),"Leadingjet_nElectrons")
h_VBS["Leadingjet_nMuons"] = df_VBSleadingjets.Histo1D(("Leadingjet_nMuons", "Leadingjet_nMuons", 10, 0, 10),"Leadingjet_nMuons")
h_VBS["Leadingjet_neEmEF"] = df_VBSleadingjets.Histo1D(("Leadingjet_neEmEF", "Leadingjet_neEmEF", 100, 0, 1),"Leadingjet_neEmEF")
h_VBS["Leadingjet_neHEF"] = df_VBSleadingjets.Histo1D(("Leadingjet_neHEF", "Leadingjet_neHEF", 100, 0, 1),"Leadingjet_neHEF")
h_VBS["Leadingjet_qgl"] = df_VBSleadingjets.Histo1D(("Leadingjet_qgl", "Leadingjet_qgl", 100, 0, 1),"Leadingjet_qgl")
h_VBS["SubLeadingjet_pt"] = df_VBSleadingjets.Histo1D(("SubLeadingjet_pt", "SubLeadingjet_pt", 100, 0, 1000),"SubLeadingjet_pt")
h_VBS["SubLeadingjet_eta"] = df_VBSleadingjets.Histo1D(("SubLeadingjet_eta", "SubLeadingjet_eta", 100, -5, 5),"SubLeadingjet_eta")
h_VBS["SubLeadingjet_phi"] = df_VBSleadingjets.Histo1D(("SubLeadingjet_phi", "SubLeadingjet_phi", 100, -3.5, 3.5),"SubLeadingjet_phi")
h_VBS["SubLeadingjet_mass"] = df_VBSleadingjets.Histo1D(("SubLeadingjet_mass", "SubLeadingjet_mass", 100, 0, 200),"SubLeadingjet_mass")
h_VBS["SubLeadingjet_area"] = df_VBSleadingjets.Histo1D(("SubLeadingjet_area", "SubLeadingjet_area", 100, 0, 0.8),"SubLeadingjet_area")
h_VBS["SubLeadingjet_nConstituents"] = df_VBSleadingjets.Histo1D(("SubLeadingjet_nConstituents", "SubLeadingjet_nConstituents", 100, 0, 150),"SubLeadingjet_nConstituents")
h_VBS["SubLeadingjet_nElectrons"] = df_VBSleadingjets.Histo1D(("SubLeadingjet_nElectrons", "SubLeadingjet_nElectrons", 10, 0, 10),"SubLeadingjet_nElectrons")
h_VBS["SubLeadingjet_nMuons"] = df_VBSleadingjets.Histo1D(("SubLeadingjet_nMuons", "SubLeadingjet_nMuons", 10, 0, 10),"SubLeadingjet_nMuons")
h_VBS["SubLeadingjet_neEmEF"] = df_VBSleadingjets.Histo1D(("SubLeadingjet_neEmEF", "SubLeadingjet_neEmEF", 100, 0, 1),"SubLeadingjet_neEmEF")
h_VBS["SubLeadingjet_neHEF"] = df_VBSleadingjets.Histo1D(("SubLeadingjet_neHEF", "SubLeadingjet_neHEF", 100, 0, 1),"SubLeadingjet_neHEF")
h_VBS["SubLeadingjet_qgl"] = df_VBSleadingjets.Histo1D(("SubLeadingjet_qgl", "SubLeadingjet_qgl", 100, 0, 1),"SubLeadingjet_qgl")
h_VBS["mjj"] = df_VBSleadingjets.Histo1D(("mjj", "mjj", 100, 0, 10000),"mjj")

allCutsReport = df.Report()
allCutsReport.Print()

df = ROOT.RDataFrame("Events", "14659586-1A5E-8E40-921D-05D85DBCA1DC.root")
df_atleast2Jets = df.Filter("nJet>2", "At least two jets")
df_GoodJets = df_atleast2Jets.Define("GoodJets_idx", "GoodJets(Jet_jetId, Jet_eta, Jet_pt, Jet_puId)")
df_atleast2GoodJets = df_GoodJets.Filter("atleast2GoodJets(GoodJets_idx)", "At least two good jets")
df_VBSjets = df_atleast2GoodJets.Define("VBSJet_idx", "SelectVBSJets_invmass(Jet_pt, Jet_eta, Jet_phi, Jet_mass, GoodJets_idx)")
df_2VBSjets = df_VBSjets.Filter("VBSJet_idx[0] != VBSJet_idx[1]", "2 VBS jets")

df_VBSleadingjets = df_2VBSjets.Define("Leadingjet_pt", "GetLeading(Jet_pt, VBSJet_idx)")\
                               .Define("Leadingjet_eta", "GetLeading(Jet_eta, VBSJet_idx)")\
                               .Define("Leadingjet_phi", "GetLeading(Jet_phi, VBSJet_idx)")\
                               .Define("Leadingjet_mass", "GetLeading(Jet_mass, VBSJet_idx)")\
                               .Define("Leadingjet_area", "GetLeading(Jet_area, VBSJet_idx)")\
                               .Define("Leadingjet_nConstituents", "GetLeading(Jet_nConstituents, VBSJet_idx)")\
                               .Define("Leadingjet_nElectrons", "GetLeading(Jet_nElectrons, VBSJet_idx)")\
                               .Define("Leadingjet_nMuons", "GetLeading(Jet_nMuons, VBSJet_idx)")\
                               .Define("Leadingjet_neEmEF", "GetLeading(Jet_neEmEF, VBSJet_idx)")\
                               .Define("Leadingjet_neHEF", "GetLeading(Jet_neHEF, VBSJet_idx)")\
                               .Define("Leadingjet_qgl", "GetLeading(Jet_qgl, VBSJet_idx)")\
                               .Define("SubLeadingjet_pt", "GetSubLeading(Jet_pt, VBSJet_idx)")\
                               .Define("SubLeadingjet_eta", "GetSubLeading(Jet_eta, VBSJet_idx)")\
                               .Define("SubLeadingjet_phi", "GetSubLeading(Jet_phi, VBSJet_idx)")\
                               .Define("SubLeadingjet_mass", "GetSubLeading(Jet_mass, VBSJet_idx)")\
                               .Define("SubLeadingjet_area", "GetSubLeading(Jet_area, VBSJet_idx)")\
                               .Define("SubLeadingjet_nConstituents", "GetSubLeading(Jet_nConstituents, VBSJet_idx)")\
                               .Define("SubLeadingjet_nElectrons", "GetSubLeading(Jet_nElectrons, VBSJet_idx)")\
                               .Define("SubLeadingjet_nMuons", "GetSubLeading(Jet_nMuons, VBSJet_idx)")\
                               .Define("SubLeadingjet_neEmEF", "GetSubLeading(Jet_neEmEF, VBSJet_idx)")\
                               .Define("SubLeadingjet_neHEF", "GetSubLeading(Jet_neHEF, VBSJet_idx)")\
                               .Define("SubLeadingjet_qgl", "GetSubLeading(Jet_qgl, VBSJet_idx)")\
                               .Define("mjj", "GetInvMass(Jet_pt, Jet_eta, Jet_phi, Jet_mass, VBSJet_idx)")\


h_QCD = {}

#h = df_atleast2GoodJets.Histo1D("nJet")
h_QCD["Leadingjet_pt"] = df_VBSleadingjets.Histo1D(("Leadingjet_pt", "Leadingjet_pt", 100, 0, 1000),"Leadingjet_pt")
h_QCD["Leadingjet_eta"] = df_VBSleadingjets.Histo1D(("Leadingjet_eta", "Leadingjet_eta", 100, -5, 5),"Leadingjet_eta")
h_QCD["Leadingjet_phi"] = df_VBSleadingjets.Histo1D(("Leadingjet_phi", "Leadingjet_phi", 100, -3.5, 3.5),"Leadingjet_phi")
h_QCD["Leadingjet_mass"] = df_VBSleadingjets.Histo1D(("Leadingjet_mass", "Leadingjet_mass", 100, 0, 200),"Leadingjet_mass")
h_QCD["Leadingjet_area"] = df_VBSleadingjets.Histo1D(("Leadingjet_area", "Leadingjet_area", 100, 0, 0.8),"Leadingjet_area")
h_QCD["Leadingjet_nConstituents"] = df_VBSleadingjets.Histo1D(("Leadingjet_nConstituents", "Leadingjet_nConstituents", 100, 0, 150),"Leadingjet_nConstituents")
h_QCD["Leadingjet_nElectrons"] = df_VBSleadingjets.Histo1D(("Leadingjet_nElectrons", "Leadingjet_nElectrons", 10, 0, 10),"Leadingjet_nElectrons")
h_QCD["Leadingjet_nMuons"] = df_VBSleadingjets.Histo1D(("Leadingjet_nMuons", "Leadingjet_nMuons", 10, 0, 10),"Leadingjet_nMuons")
h_QCD["Leadingjet_neEmEF"] = df_VBSleadingjets.Histo1D(("Leadingjet_neEmEF", "Leadingjet_neEmEF", 100, 0, 1),"Leadingjet_neEmEF")
h_QCD["Leadingjet_neHEF"] = df_VBSleadingjets.Histo1D(("Leadingjet_neHEF", "Leadingjet_neHEF", 100, 0, 1),"Leadingjet_neHEF")
h_QCD["Leadingjet_qgl"] = df_VBSleadingjets.Histo1D(("Leadingjet_qgl", "Leadingjet_qgl", 100, 0, 1),"Leadingjet_qgl")
h_QCD["SubLeadingjet_pt"] = df_VBSleadingjets.Histo1D(("SubLeadingjet_pt", "SubLeadingjet_pt", 100, 0, 1000),"SubLeadingjet_pt")
h_QCD["SubLeadingjet_eta"] = df_VBSleadingjets.Histo1D(("SubLeadingjet_eta", "SubLeadingjet_eta", 100, -5, 5),"SubLeadingjet_eta")
h_QCD["SubLeadingjet_phi"] = df_VBSleadingjets.Histo1D(("SubLeadingjet_phi", "SubLeadingjet_phi", 100, -3.5, 3.5),"SubLeadingjet_phi")
h_QCD["SubLeadingjet_mass"] = df_VBSleadingjets.Histo1D(("SubLeadingjet_mass", "SubLeadingjet_mass", 100, 0, 200),"SubLeadingjet_mass")
h_QCD["SubLeadingjet_area"] = df_VBSleadingjets.Histo1D(("SubLeadingjet_area", "SubLeadingjet_area", 100, 0, 0.8),"SubLeadingjet_area")
h_QCD["SubLeadingjet_nConstituents"] = df_VBSleadingjets.Histo1D(("SubLeadingjet_nConstituents", "SubLeadingjet_nConstituents", 100, 0, 150),"SubLeadingjet_nConstituents")
h_QCD["SubLeadingjet_nElectrons"] = df_VBSleadingjets.Histo1D(("SubLeadingjet_nElectrons", "SubLeadingjet_nElectrons", 10, 0, 10),"SubLeadingjet_nElectrons")
h_QCD["SubLeadingjet_nMuons"] = df_VBSleadingjets.Histo1D(("SubLeadingjet_nMuons", "SubLeadingjet_nMuons", 10, 0, 10),"SubLeadingjet_nMuons")
h_QCD["SubLeadingjet_neEmEF"] = df_VBSleadingjets.Histo1D(("SubLeadingjet_neEmEF", "SubLeadingjet_neEmEF", 100, 0, 1),"SubLeadingjet_neEmEF")
h_QCD["SubLeadingjet_neHEF"] = df_VBSleadingjets.Histo1D(("SubLeadingjet_neHEF", "SubLeadingjet_neHEF", 100, 0, 1),"SubLeadingjet_neHEF")
h_QCD["SubLeadingjet_qgl"] = df_VBSleadingjets.Histo1D(("SubLeadingjet_qgl", "SubLeadingjet_qgl", 100, 0, 1),"SubLeadingjet_qgl")
h_QCD["mjj"] = df_VBSleadingjets.Histo1D(("mjj", "mjj", 100, 0, 10000),"mjj")

allCutsReport = df.Report()
allCutsReport.Print()

variables = ["Leadingjet_pt",
"Leadingjet_eta",
"Leadingjet_phi",
"Leadingjet_mass",
"Leadingjet_area",
"Leadingjet_nConstituents",
"Leadingjet_nElectrons",
"Leadingjet_nMuons",
"Leadingjet_neEmEF",
"Leadingjet_neHEF",
"Leadingjet_qgl",
"SubLeadingjet_pt",
"SubLeadingjet_eta",
"SubLeadingjet_phi",
"SubLeadingjet_mass",
"SubLeadingjet_area",
"SubLeadingjet_nConstituents",
"SubLeadingjet_nElectrons",
"SubLeadingjet_nMuons",
"SubLeadingjet_neEmEF",
"SubLeadingjet_neHEF",
"SubLeadingjet_qgl",
"mjj",]

canvases = []
legends = []

for i, v in enumerate(variables):
    canvases.append(ROOT.TCanvas())
    variable = variables[i]
    h_qcd = h_QCD[variable]
    h_vbs = h_VBS[variable]
    h_qcd.Scale( 1./h_qcd.Integral(),"WIDTH")
    h_vbs.Scale( 1./h_vbs.Integral(),"WIDTH")
    maximum = max(h_qcd.GetMaximum(), h_vbs.GetMaximum())
    h_qcd.SetLineColor(3)
    h_qcd.GetYaxis().SetRangeUser(0,maximum)
    h_qcd.SetStats(False)
    h_vbs.SetStats(False)
    h_vbs.GetYaxis().SetRangeUser(0,maximum)
    h1 = h_vbs.DrawCopy('HIST')
    h2 = h_qcd.DrawCopy('HIST SAME')
    legends.append(ROOT.TLegend(.73,.32,.97,.53))
    legends[i].SetBorderSize(0)
    legends[i].SetFillColor(0)
    legends[i].SetFillStyle(0)
    legends[i].SetTextFont(42)
    legends[i].SetTextSize(0.035)
    legends[i].AddEntry(h1, "VBS", "L")
    legends[i].AddEntry(h2, "QCD", "L")
    legends[i].Draw()
    canvases[i].Draw()
    #c.Show()