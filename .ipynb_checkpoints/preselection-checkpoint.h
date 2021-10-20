/// Header file with functions needed to execute the Python version
/// preselection step of the analysis. The header is declared to the
/// ROOT C++ interpreter prior to the start of the analysis via the
/// `ROOT.gInterpreter.Declare()` function.
///

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TFile.h"
#include "TH2D.h"
#include "TLatex.h"
#include "Math/Vector4D.h"
#include "TStyle.h"
#include <map>

using namespace ROOT::VecOps;
using RNode = ROOT::RDF::RNode;
using rvec_f = const RVec<float> &;
using rvec_i = const RVec<int> &;
using rvec_b = const RVec<bool> &;

RVec<float> MHT_pt_phi(rvec_f Electron_pt, rvec_f Electron_eta, rvec_f Electron_phi, rvec_f Electron_mass, rvec_f Electron_miniPFRelIso_all, rvec_f Muon_pt, rvec_f Muon_eta, rvec_f Muon_phi, rvec_f Muon_mass, rvec_f Muon_miniPFRelIso_all, rvec_f Jet_pt, rvec_f Jet_eta, rvec_f Jet_phi, rvec_f Jet_mass, rvec_i Jet_muonIdx1, rvec_i Jet_muonIdx2, rvec_i Jet_electronIdx1, rvec_i Jet_electronIdx2, int nJet){
    
    RVec<float> MHT_pt_phi(2);
    
    ROOT::Math::PtEtaPhiMVector mht;
    
    for (size_t j = 0; j < Muon_pt.size(); j++){
        if (Muon_pt[j] > 20 && Muon_miniPFRelIso_all[j] < 0.2){
            ROOT::Math::PtEtaPhiMVector p(Muon_pt[j], Muon_eta[j], Muon_phi[j], Muon_mass[j]);
            mht = mht + p;
        }
    }
    for (size_t j = 0; j < Electron_pt.size(); j++){
        if (Electron_pt[j] > 20 && Electron_miniPFRelIso_all[j] < 0.2){
            ROOT::Math::PtEtaPhiMVector p(Electron_pt[j], Electron_eta[j], Electron_phi[j], Electron_mass[j]);
            mht = mht + p;
        }
    }
    //RVec<int> goodjets(nJets)
    for (size_t i = 0; i < Jet_pt.size(); i++) {
        if (Jet_pt[i] > 40) continue;
        if (Jet_muonIdx1[i] != -1 and Jet_muonIdx1[i] < nJet){
            if (Muon_pt[Jet_muonIdx1[i]] > 20 && Muon_miniPFRelIso_all[Jet_muonIdx1[i]] < 0.2) continue;  //prefer the muon
        }
        if (Jet_muonIdx2[i] != -1 and Jet_muonIdx2[i] < nJet){
            if (Muon_pt[Jet_muonIdx2[i]] > 20 && Muon_miniPFRelIso_all[Jet_muonIdx2[i]] < 0.2) continue;  //prefer the muon
        }
        if (Jet_electronIdx1[i] != -1 and Jet_electronIdx1[i] < nJet){
            if (Electron_pt[Jet_electronIdx1[i]] > 20 && Electron_miniPFRelIso_all[Jet_electronIdx1[i]] < 0.2) continue; //prefer the electron
        }
        if (Jet_electronIdx2[i] != -1 and Jet_electronIdx2[i] < nJet){
            if (Electron_pt[Jet_electronIdx2[i]] > 20 && Electron_miniPFRelIso_all[Jet_electronIdx2[i]] < 0.2) continue; //prefer the electron
        }
        //goodjets[i] = 1;
        ROOT::Math::PtEtaPhiMVector p(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
        mht = mht + p;    
    }
    
    MHT_pt_phi[0] = mht.Pt();
    MHT_pt_phi[1] = -mht.Phi();
    
    return MHT_pt_phi;
}

float getFirst(rvec_f a){
    return a[0];
}

float getSecond(rvec_f a){
    return a[1];
}


float GetYear(unsigned int slot, const ROOT::RDF::RSampleInfo &id){
    Int_t year;
    if (id.Contains("RunIISummer16NanoAODv7")){
        year = 2016;
    } else if (id.Contains("RunIIFall17NanoAODv7")){
        year = 2017;
    } else if (id.Contains("RunIIAutumn18NanoAODv7")){
        year = 2018;
    }
    return year;
}

float MET_HLT_Filter(Int_t Year, Bool_t Flag_goodVertices, Bool_t Flag_HBHENoiseFilter, Bool_t Flag_HBHENoiseIsoFilter, Bool_t Flag_EcalDeadCellTriggerPrimitiveFilter, Bool_t Flag_BadPFMuonFilter, Bool_t Flag_globalSuperTightHalo2016Filter, Bool_t HLT_Ele27_WPTight_Gsf, Bool_t HLT_Ele32_WPTight_Gsf, Bool_t HLT_IsoMu24, Bool_t HLT_IsoMu27, Bool_t HLT_Mu50, Bool_t HLT_Ele35_WPTight_Gsf, Bool_t HLT_Ele32_WPTight_Gsf_L1DoubleEG, Bool_t HLT_Photon200){
    bool good_MET = Flag_goodVertices && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter;
    bool good_HLT;
    bool HLT_IsoTkMu24 = true;
    if (Year == 2016){
        good_HLT = (HLT_Ele27_WPTight_Gsf || HLT_Ele32_WPTight_Gsf || HLT_IsoMu24 || HLT_IsoTkMu24) && Flag_globalSuperTightHalo2016Filter;
    } else if (Year == 2017){
        good_HLT = (HLT_IsoMu27 || HLT_Mu50 || HLT_Ele35_WPTight_Gsf || HLT_Ele32_WPTight_Gsf_L1DoubleEG || HLT_Photon200); // or HLT.PFHT250 or HLT.PFHT350)
    } else if (Year == 2018){
        good_HLT = (HLT_IsoMu27 || HLT_Mu50 || HLT_Ele35_WPTight_Gsf || HLT_Ele32_WPTight_Gsf_L1DoubleEG || HLT_Photon200); // or HLT.PFHT250 or HLT.PFHT350)
    }
    return good_MET && good_HLT;
}

//float MET_HLT_Filter(unsigned int slot, const ROOT::RDF::RSampleInfo &id, Bool_t Flag_goodVertices, Bool_t Flag_HBHENoiseFilter, Bool_t Flag_HBHENoiseIsoFilter, Bool_t Flag_EcalDeadCellTriggerPrimitiveFilter, Bool_t Flag_BadPFMuonFilter, Bool_t Flag_globalSuperTightHalo2016Filter, Bool_t HLT_Ele27_WPTight_Gsf, Bool_t HLT_Ele32_WPTight_Gsf, Bool_t HLT_IsoMu24, Bool_t HLT_IsoTkMu24, Bool_t HLT_IsoMu27, Bool_t HLT_Mu50, Bool_t HLT_Ele35_WPTight_Gsf, Bool_t HLT_Ele32_WPTight_Gsf_L1DoubleEG, Bool_t HLT_Photon200){
//    bool good_MET = Flag_goodVertices && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter;
//    bool good_HLT;
//    if (id.Contains("RunIISummer16NanoAODv7")){
//        good_HLT = (HLT_Ele27_WPTight_Gsf || HLT_Ele32_WPTight_Gsf || HLT_IsoMu24 || HLT_IsoTkMu24) && Flag_globalSuperTightHalo2016Filter;
//    } else if (id.Contains("RunIIFall17NanoAODv7")){
//        good_HLT = (HLT_IsoMu27 || HLT_Mu50 || HLT_Ele35_WPTight_Gsf || HLT_Ele32_WPTight_Gsf_L1DoubleEG || HLT_Photon200); // or HLT.PFHT250 or HLT.PFHT350)
//    } else if (id.Contains("RunIIAutumn18NanoAODv7")){
//        good_HLT = (HLT_IsoMu27 || HLT_Mu50 || HLT_Ele35_WPTight_Gsf || HLT_Ele32_WPTight_Gsf_L1DoubleEG || HLT_Photon200); // or HLT.PFHT250 or HLT.PFHT350)
//    }
//    return good_MET && good_HLT;
//}