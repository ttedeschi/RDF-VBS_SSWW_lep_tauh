import ROOT
import os
os.environ["X509_USER_PROXY"] = "/eos/home-t/ttedesch/SWAN_projects/RDataFrame_porting/proxy"
path = "root://cms-xrd-global.cern.ch//store/user/apiccine/VBS/"

import ROOT
import os

ROOT.ROOT.EnableImplicitMT()
ROOT.gInterpreter.Declare('#include "postselection.h"')

pt_threshold = 0
# Provided by SWAN, but it could be created manually
#sc

# Use a Spark RDataFrame
#RDataFrame = ROOT.RDF.Experimental.Distributed.Spark.RDataFrame

'''
df = RDataFrame("Events",
                samples,
                npartitions=16,
                sparkcontext=sc)
'''
df = ROOT.RDataFrame("Events", "root://cms-xrd-global.cern.ch//store/user/apiccine/VBS/VBS_SSWW_LL_polarization_TuneCP5_13TeV-madgraph-pythia8/VBS_SSWW_LL_SM_2017/210420_123836/0000/tree_hadd_2.root")
df_goodvertex = df.Filter("Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && Flag_ecalBadCalibFilterV2", "Good vertex selection")
df_trigger = df_goodvertex.Filter("HLT_IsoMu27 || HLT_Mu50 || HLT_Ele35_WPTight_Gsf || HLT_Ele32_WPTight_Gsf_L1DoubleEG || HLT_Photon200", "Leptonic trigger")
df_passEleorMu
df_atleast2Jets = df_trigger.Filter("nJet>2", "At least two jets")
df_GoodJets = df_atleast2Jets.Define("GoodJets_idx", "GoodJets(Jet_jetId, Jet_eta, Jet_pt, Jet_puId)")
df_atleast2GoodJets = df_GoodJets.Filter("atleast2GoodJets(GoodJets_idx)", "At least two good jets")
df_VBSjets = df_atleast2GoodJets.Define("VBSJet_idx", "SelectVBSJets_invmass(Jet_pt, Jet_eta, Jet_phi, Jet_mass, GoodJets_idx)")
df_2VBSjets = df_VBSjets.Filter("VBSJet_idx[0] != VBSJet_idx[1]", "2 VBS jets")
df_VBSleadingjets = df_2VBSjets.Define("Leadingjet_pt", "GetLeading(Jet_pt, VBSJet_idx)")


h = df_VBSleadingjets.Histo1D(("Leadingjet_pt", "Leadingjet_pt", 100, 0, 1000),"Leadingjet_pt")


#df_? = df_goodvertex.Filter("Flag_eeBadScFilter && dataset == 'Data'", ?)
#df_atleast2Jets = df_goodvertex.Filter("nJet>2", "At least two jets")
#df_GoodJets = df_atleast2Jets.Define("GoodJets_idx", "GoodJets(Jet_jetId, Jet_eta, Jet_pt, Jet_puId)")
#df_atleast2GoodJets = df_GoodJets.Filter("atleast2GoodJets(GoodJets_idx)", "At least two good jets")
#df_VBSjets = df_atleast2GoodJets.Define("VBSJet_idx", "SelectVBSJets_invmass(Jet_pt, Jet_eta, Jet_phi, Jet_mass, GoodJets_idx)")

##### e/mu + tau process
df_selectElectron = df_atleast2GoodJets.Define("Electron_idx", "SelectElectron(Electron_pt, Electron_eta, Electron_phi, Electron_tightId, Electron_looseId, Electron_jetRelIso, Electron_pfRelIso04_all, Electron_mvaFall17V2Iso_WPL, mvaFall17V2Iso_WP90, Jet_pt, Jet_eta, Jet_phi, VBSJets_idx)")
df_selectMuon = df_selectElectron.Define("Muon_idx", "SelectMuon(Muon_pt, Muon_eta, Muon_phi, Muon_tightId, Muon_looseId, Muon_pfRelIso04_all, Jet_pt, Jet_eta, Jet_phi, VBSJets_idx)")
df_atLeast1Lepton = df_selectMuon.Filter("Electron_idx[1] != -1 || Muon_idx[1] != -1", "At least 1 at-least-loose lepton")
df_goodLepton = df_atLeast1Lepton.Define("GoodLeptonFamily", "DetermineGoodLepton(HLT_IsoMu27, HLT_Mu50, HLT_Ele35_WPTight_Gsf, HLT_Ele32_WPTight_Gsf_L1DoubleEG, HLT_Photon200, HLT_PFHT250, HLT_PFHT350, Electron_idx, Electron_pt, Electron_eta, Electron_mvaFall17V2Iso_WPL, Electron_jetRelIso, Muon_idx, Muon_pt, Muon_eta, Muon_pfRelIso04_all)")
df_compatibleLeptons = df_goodLepton.Filter("GoodLeptonFamily != -1 ", "Filter on leptons") 



df_selectTau = df_atleast1Lepton.Define("Lepton_idx", "SelectAndVetoTaus(Tau_eta, Tau_phi,\
                                         Tau_idDeepTau2017v2p1VSjet, rvec_f Tau_idDeepTau2017v2p1VSe, rvec_f Tau_idDeepTau2017v2p1VSmu, rvec_f Tau_idDecayModeNewDMs,\
                                         Lepton_eta, rvec_f Lepton_phi, rvec_i Lepton_pdgId, rvec_i Lepton_idx, rvec_f Jet_eta, rvec_f Jet_phi, rvec_i Jet_idx)")

df_1tau = df_selectTau.Filter("Lepton_idx[1] != -1", "Exactly 1 at least loose Tau")



