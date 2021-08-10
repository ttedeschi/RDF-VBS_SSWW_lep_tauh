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


df_goodvertex = df.Filter("Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_HBHENoiseFilter \
                           && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter \
                           && Flag_BadPFMuonFilter && Flag_ecalBadCalibFilterV2", 
                           "Good vertex selection")
#df_? = df_goodvertex.Filter("Flag_eeBadScFilter && dataset == 'Data'", ?)

df_atleast2Jets = df_goodvertex.Filter("nJet>2", "At least two jets")
df_GoodJets = df_atleast2Jets.Define("GoodJets_idx", "GoodJets(Jet_jetId, Jet_eta, Jet_pt, Jet_puId)")
df_atleast2GoodJets = df_GoodJets.Filter("atleast2GoodJets(GoodJets_idx)", "At least two good jets")
df_VBSjets = df_atleast2GoodJets.Define("VBSJet_idx", "SelectVBSJets_invmass(Jet_pt, Jet_eta, Jet_phi, Jet_mass, GoodJets_idx)")

##### e/mu + tau process
df_selectLepton = df_atleast2GoodJets.Define("Lepton_idx", 
                                               "SelectLepton(lepton_pt, lepton_eta, lepton_phi,lepton_pdgId, lepton_tightId, lepton_looseId, \
                                                lepton_jetRelIso, lepton_pfRelIso04_all,\
                                                mvaFall17V2Iso_WPL, mvaFall17V2Iso_WP90, \
                                                jet_pt, jet_eta, jet_phi, VBSJets_idx)")

df_atLeast1Lepton = df_selectLepton.Filter("Lepton_idx[1] != -1", "At least 1 at least loose lepton")

df_selectTau = df_atleast1Lepton.Define("Lepton_idx", "SelectAndVetoTaus(Tau_eta, Tau_phi,\
                                         Tau_idDeepTau2017v2p1VSjet, rvec_f Tau_idDeepTau2017v2p1VSe, rvec_f Tau_idDeepTau2017v2p1VSmu, rvec_f Tau_idDecayModeNewDMs,\
                                         Lepton_eta, rvec_f Lepton_phi, rvec_i Lepton_pdgId, rvec_i Lepton_idx, rvec_f Jet_eta, rvec_f Jet_phi, rvec_i Jet_idx)")

df_1tau = df_selectTau.Filter("Lepton_idx[1] != -1", "Exactly 1 at least loose Tau")



