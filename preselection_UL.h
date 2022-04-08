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
#include "TH2F.h"
#include "TLatex.h"
#include "Math/Vector4D.h"
#include "TStyle.h"
#include <map>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <TH1.h>
#include <TH2.h>
#include <TH2F.h>
#include <cmath>
#include <curl/curl.h>
#include <stdio.h>

using namespace ROOT::VecOps;
using RNode = ROOT::RDF::RNode;
using rvec_f = const RVec<float> &;
using rvec_i = const RVec<int> &;
using rvec_b = const RVec<bool> &;

const string remote_storage = "https://vbs-pg-support.web.cern.ch/nanoAOD-tools/";

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
        //if (Jet_pt[i] > 40) continue;
        if (Jet_pt[i] < 40) continue;
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


//float GetYear(unsigned int slot, const ROOT::RDF::RSampleInfo &id){
//    Int_t year;
//    if (id.Contains("RunIISummer16NanoAODv7")){
//        year = 2016;
//    } else if (id.Contains("RunIIFall17NanoAODv7")){
//        year = 2017;
//    } else if (id.Contains("RunIIAutumn18NanoAODv7")){
//        year = 2018;
//    }
//    return year;
//}

string GetYear(unsigned int slot, const ROOT::RDF::RSampleInfo &id){
    string year="2017";
    if (id.Contains("RunIISummer16NanoAODv7")){
        year = "2016";
    } else if (id.Contains("RunIIFall17NanoAODv7")){
        year = "2017";
    } else if (id.Contains("RunIIAutumn18NanoAODv7")){
        year = "2018";
    }
    return year;
}

map<string, int> dataset_map = {{"/WWTo2L2Nu_DoubleScattering_13TeV-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM", 193}, {"/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_EXT_102X_mc2017_realistic_v8-v1/NANOAODSIM", 194}, {"/WWZTo3L1Nu2Q_4f_TuneCP5_13TeV_amcatnlo_pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM", 195}, {"/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_EXT_102X_mc2017_realistic_v8-v1/NANOAODSIM", 196}, {"/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_EXT_102X_mc2017_realistic_v8-v1/NANOAODSIM", 197}, {"/WWG_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_EXT_102X_mc2017_realistic_v8-v1/NANOAODSIM", 198}, {"/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_EXT_102X_mc2017_realistic_v8_ext1-v1/NANOAODSIM", 165}, {"/TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8_ext1-v1/NANOAODSIM", 166}, {"/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_EXT_102X_mc2017_realistic_v8-v1/NANOAODSIM", 167}, {"/TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM", 168}, {"/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_EXT_102X_mc2017_realistic_v8-v1/NANOAODSIM", 169}, {"/tZq_ll_4f_ckm_NLO_TuneCP5_PSweights_13TeV-amcatnlo-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_new_pmx_102X_mc2017_realistic_v8-v1/NANOAODSIM", 170}, {"/WpWpJJ_QCD_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM", 1}, {"/LLAJJ_EWK_MLL-50_MJJ-120_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8_ext1-v1/NANOAODSIM", 162}, {"/LNuAJJ_EWK_MJJ-120_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM", 163}, {"/ZZTo2L2Nu_13TeV_powheg_pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM", 115}, {"/ZZJJTo4L_EWK_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM", 116}, {"/ZZJJTo4L_QCD_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM", 117}, {"/GluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_EXT_102X_mc2017_realistic_v8-v1/NANOAODSIM", 118}, {"/GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_EXT_102X_mc2017_realistic_v8-v1/NANOAODSIM", 119}, {"/GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_EXT_102X_mc2017_realistic_v8-v1/NANOAODSIM", 120}, {"/GluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_EXT_102X_mc2017_realistic_v8-v1/NANOAODSIM", 121}, {"/GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_EXT_102X_mc2017_realistic_v8-v1/NANOAODSIM", 122}, {"/GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_EXT_102X_mc2017_realistic_v8-v1/NANOAODSIM", 123}, {"/GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_EXT_102X_mc2017_realistic_v8-v1/NANOAODSIM", 124}, {"/GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM", 125}, {"/WZ_TuneCP5_13TeV-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_PU2017_EXT_102X_mc2017_realistic_v8-v1/NANOAODSIM", 199}, {"/VBS_SSWW_LL_polarization_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM", 3}, {"/VBS_SSWW_TL_polarization_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM", 4}, {"/VBS_SSWW_TT_polarization_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM", 5}, {"/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_new_pmx_102X_mc2017_realistic_v8-v1/NANOAODSIM", 151}, {"/DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_new_pmx_102X_mc2017_realistic_v8_ext1-v1/NANOAODSIM", 152}, {"/DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM", 153}, {"/DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_v2_102X_mc2017_realistic_v8-v1/NANOAODSIM", 154}, {"/DYJetsToLL_M-5to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM", 155}, {"/WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM", 172}, {"/GluGluToWWToENEN_13TeV_MCFM701_pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM", 173}, {"/GluGluToWWToENMN_13TeV_MCFM701_pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM", 174}, {"/GluGluToWWToENTN_13TeV_MCFM701_pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM", 175}, {"/GluGluToWWToMNEN_13TeV_MCFM701_pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM", 176}, {"/GluGluToWWToMNMN_13TeV_MCFM701_pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM", 177}, {"/GluGluToWWToMNTN_13TeV_MCFM701_pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM", 178}, {"/GluGluToWWToTNEN_13TeV_MCFM701_pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM", 179}, {"/GluGluToWWToTNMN_13TeV_MCFM701_pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM", 180}, {"/GluGluToWWToTNTN_13TeV_MCFM701_pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM", 181}, {"/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM", 182}, {"/GluGluHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM", 183}, {"/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8_ext3-v1/NANOAODSIM", 185}, {"/GluGluHToTauTau_M125_13TeV_powheg_pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_new_pmx_102X_mc2017_realistic_v8-v1/NANOAODSIM", 186}, {"/VBFHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM", 187}, {"/VBFHToTauTau_M125_13TeV_amcatnloFXFX_pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM", 188}, {"/ttHToNonbb_M125_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM", 189}, {"/VHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM", 190}, {"/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_new_pmx_102X_mc2017_realistic_v8-v1/NANOAODSIM", 130}, {"/SingleElectron/Run2017B-02Apr2020-v1/NANOAOD", 213}, {"/SingleElectron/Run2017C-02Apr2020-v1/NANOAOD", 214}, {"/SingleElectron/Run2017D-02Apr2020-v1/NANOAOD", 215}, {"/SingleElectron/Run2017E-02Apr2020-v1/NANOAOD", 216}, {"/SingleElectron/Run2017F-02Apr2020-v1/NANOAOD", 217}, {"/SingleMuon/Run2017B-02Apr2020-v1/NANOAOD", 202}, {"/SingleMuon/Run2017C-02Apr2020-v1/NANOAOD", 203}, {"/SingleMuon/Run2017D-02Apr2020-v1/NANOAOD", 204}, {"/SingleMuon/Run2017E-02Apr2020-v1/NANOAOD", 205}, {"/SingleMuon/Run2017F-02Apr2020-v1/NANOAOD", 206}};

string tokenize(string s, string del = " ")
{
    int start = 0;
    int end = s.find(del);
    int i = 0;
    int start_ = 0;
    int end_ = 0;
    string string1, string2, string3, string4;
    while (end != -1) {
        //cout << s.substr(start, end - start) << endl;
        start = end + del.size();
        end = s.find(del, start);
        if (i == 5) string1 = s.substr(start, end - start);
        else if (i==6) string2 = s.substr(start, end - start);
        else if (i==7) string3 = s.substr(start, end - start);
        else if (i==8) string4 = s.substr(start, end - start);
        //cout<<i<<" "<<s.substr(start, end - start)<<endl;
        i++;
    }
    string result;
    result += string("/") + string2 + string("/") + string1 + string("-") + string4 + string("/") + string3;
    return result;
}


int GetSample(unsigned int slot, const ROOT::RDF::RSampleInfo &id){
    return dataset_map[tokenize(id.AsString(), "/")];
}

float MET_HLT_Filter(string Year, Bool_t Flag_goodVertices, Bool_t Flag_HBHENoiseFilter, Bool_t Flag_HBHENoiseIsoFilter, Bool_t Flag_EcalDeadCellTriggerPrimitiveFilter, Bool_t Flag_BadPFMuonFilter, Bool_t Flag_globalSuperTightHalo2016Filter, Bool_t HLT_Ele27_WPTight_Gsf, Bool_t HLT_Ele32_WPTight_Gsf, Bool_t HLT_IsoMu24, Bool_t HLT_IsoMu27, Bool_t HLT_Mu50, Bool_t HLT_Ele35_WPTight_Gsf, Bool_t HLT_Ele32_WPTight_Gsf_L1DoubleEG, Bool_t HLT_Photon200){
    bool good_MET = Flag_goodVertices && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter;
    bool good_HLT;
    bool HLT_IsoTkMu24 = true;
    if (Year == "2016"){
        good_HLT = (HLT_Ele27_WPTight_Gsf || HLT_Ele32_WPTight_Gsf || HLT_IsoMu24 || HLT_IsoTkMu24) && Flag_globalSuperTightHalo2016Filter;
    } else if (Year == "2017"){
        good_HLT = (HLT_IsoMu27 || HLT_Mu50 || HLT_Ele35_WPTight_Gsf || HLT_Ele32_WPTight_Gsf_L1DoubleEG || HLT_Photon200); // or HLT.PFHT250 or HLT.PFHT350)
    } else if (Year == "2018"){
        good_HLT = (HLT_IsoMu27 || HLT_Mu50 || HLT_Ele35_WPTight_Gsf || HLT_Ele32_WPTight_Gsf_L1DoubleEG || HLT_Photon200); // or HLT.PFHT250 or HLT.PFHT350)
    }
    return good_MET && good_HLT;
}

float GetEventHT(rvec_f Jet_pt, rvec_f Jet_eta, rvec_f Jet_phi, rvec_f Jet_mass){
    ROOT::Math::PtEtaPhiMVector eventSum;
    for (size_t i = 0; i < Jet_pt.size(); i++) {
        ROOT::Math::PtEtaPhiMVector p(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
        eventSum = eventSum + p;    
    }
    return eventSum.Pt();
}

class WeightCalculatorFromHistogram {
 public:
  WeightCalculatorFromHistogram() {}
  // get the weight from the bin content of the passed histogram
  WeightCalculatorFromHistogram(TH1 *histogram, bool verbose=false) : histogram_(histogram), verbose_(verbose) {}
  // get the weight from the bin content of the ratio hist/targethist
  WeightCalculatorFromHistogram(TH1 *hist, TH1* targethist, bool norm=true, bool fixLargeWeights=true, bool verbose=false);
  ~WeightCalculatorFromHistogram() {}
  
  float getWeight(float x, float y=0) const;
  float getWeightErr(float x, float y=0) const;
  
 private:
  std::vector<float> loadVals(TH1 *hist, bool norm=true);
  TH1* ratio(TH1 *hist, TH1* targethist, bool fixLargeWgts);
  void fixLargeWeights(std::vector<float> &weights, float maxshift=0.0025,float hardmax=3);
  float checkIntegral(std::vector<float> wgt1, std::vector<float> wgt2);

  TH1* histogram_;
  std::vector<float> refvals_,targetvals_;
  bool verbose_;
  bool norm_;
};


WeightCalculatorFromHistogram::WeightCalculatorFromHistogram(TH1 *hist, TH1* targethist, bool norm, bool fixLargeWeights, bool verbose) {
  norm_ = norm;
  verbose_ = verbose;
  if(hist->GetNcells()!=targethist->GetNcells()) {
    std::cout << "ERROR! Numerator and denominator histograms have different number of bins!" << std::endl;
    histogram_=0;
  } else {
    for(int i=0; i<(int)hist->GetNcells(); ++i) {
      refvals_.push_back(hist->GetBinContent(i));
      targetvals_.push_back(targethist->GetBinContent(i));
    }
    histogram_ = ratio(hist,targethist,fixLargeWeights);
  }
}

float WeightCalculatorFromHistogram::getWeight(float x, float y) const {
  if(histogram_==NULL) {
    std::cout << "ERROR! The weights input histogram is not loaded. Returning weight 0!" << std::endl;
    return 0.;
  }
  if(!histogram_->InheritsFrom("TH2")) {
    int bin = std::max(1, std::min(histogram_->GetNbinsX(), histogram_->GetXaxis()->FindBin(x)));
    return histogram_->GetBinContent(bin);
  } else {
    int binx = std::max(1, std::min(histogram_->GetNbinsX(), histogram_->GetXaxis()->FindBin(x)));
    int biny = std::max(1, std::min(histogram_->GetNbinsY(), histogram_->GetYaxis()->FindBin(y)));
    return histogram_->GetBinContent(binx,biny);
  }
}

float WeightCalculatorFromHistogram::getWeightErr(float x, float y) const {
  if(histogram_==NULL) {
    std::cout << "ERROR! The weights input histogram is not loaded. Returning weight error 1!" << std::endl;
    return 1.;
  }
  if(!histogram_->InheritsFrom("TH2")) {
    int bin = std::max(1, std::min(histogram_->GetNbinsX(), histogram_->GetXaxis()->FindBin(x)));
    return histogram_->GetBinError(bin);
  } else {
    int binx = std::max(1, std::min(histogram_->GetNbinsX(), histogram_->GetXaxis()->FindBin(x)));
    int biny = std::max(1, std::min(histogram_->GetNbinsY(), histogram_->GetYaxis()->FindBin(y)));
    return histogram_->GetBinError(binx,biny);
  }
}

std::vector<float> WeightCalculatorFromHistogram::loadVals(TH1 *hist, bool norm) {
  int nbins=hist->GetNcells();
  std::vector<float> vals;
  for(int i=0; i<nbins; ++i) {
    double bc=hist->GetBinContent(i);
    //double val = (i>0 && bc==0 && hist->GetBinContent(i-1)>0 && hist->GetBinContent(i+1)>0) ? 0.5*(hist->GetBinContent(i-1)+hist->GetBinContent(i+1)) : bc;
    vals.push_back(std::max(bc,0.));
  }
  if(verbose_) std::cout << "Normalization of " << hist->GetName() << ": " << hist->Integral() << std::endl;
  if(norm) {
    float scale = 1.0/hist->Integral();
    for(int i=0; i<nbins; ++i) vals[i] *= scale;
  }
  return vals;
}

TH1* WeightCalculatorFromHistogram::ratio(TH1 *hist, TH1* targethist, bool fixLargeWgts) {
  TH1 *ret = (TH1*)hist->Clone("hweights");
  ret->SetDirectory(0);

  std::vector<float> vals = loadVals(hist,norm_);
  std::vector<float> targetvals = loadVals(targethist,norm_);
  std::vector<float> weights;
  int nbins = vals.size();
  if(verbose_) std::cout << "Weights for variable " << hist->GetName() << " with a number of bins equal to " << nbins << ":" << std::endl;
  for(int i=0; i<nbins; ++i) {
    float weight = vals[i] !=0 ? targetvals[i]/vals[i] : 1.;
    if(verbose_) std::cout <<  std::setprecision(3) << weight << " ";
    weights.push_back(weight);
  }
  if(verbose_) std::cout << "." << std::endl;
  if(fixLargeWgts) fixLargeWeights(weights);
  if(verbose_) std::cout << "Final weights: " << std::endl;
  for(int i=0; i<(int)weights.size(); ++i) {
    ret->SetBinContent(i,weights[i]);
    if(verbose_) std::cout << std::setprecision(3) << weights[i] << " ";
  }
  if(verbose_) std::cout << "." << std::endl;
  return ret;
}

float WeightCalculatorFromHistogram::checkIntegral(std::vector<float> wgt1, std::vector<float> wgt2) {
  float myint=0;
  float refint=0;
  for(int i=0; i<(int)wgt1.size(); ++i) {
    myint += wgt1[i]*refvals_[i];
    refint += wgt2[i]*refvals_[i];
  }
  return (myint-refint)/refint;
}

void WeightCalculatorFromHistogram::fixLargeWeights(std::vector<float> &weights, float maxshift,float hardmax) {
  float maxw = std::min(*(std::max_element(weights.begin(),weights.end())),float(5.));
  std::vector<float> cropped;
  while (maxw > hardmax) {
    cropped.clear();  
    for(int i=0; i<(int)weights.size(); ++i) cropped.push_back(std::min(maxw,weights[i]));
    float shift = checkIntegral(cropped,weights);
    if(verbose_) std::cout << "For maximum weight " << maxw << ": integral relative change: " << shift << std::endl;
    if(fabs(shift) > maxshift) break;
    maxw *= 0.95;
  }
  maxw /= 0.95;
  if (cropped.size()>0) {
      for(int i=0; i<(int)weights.size(); ++i) cropped[i] = std::min(maxw,weights[i]);
      float normshift = checkIntegral(cropped,weights);
      for(int i=0; i<(int)weights.size(); ++i) weights[i] = cropped[i]*(1-normshift);
  }
}

class LeptonEfficiencyCorrector {
 public:

  LeptonEfficiencyCorrector() {effmaps_.clear();}
  LeptonEfficiencyCorrector(std::vector<std::string> files, std::vector<std::string> histos);
  ~LeptonEfficiencyCorrector() {}

  void setLeptons(int nLep, int *lepPdgId, float *lepPt, float *lepEta);

  float getSF(int pdgid, float pt, float eta);
  float getSFErr(int pdgid, float pt, float eta);
  const std::vector<float> & run();

private:
  std::vector<TH2F *> effmaps_;
  std::vector<float> ret_;
  int nLep_;
  float *Lep_eta_, *Lep_pt_;
  int *Lep_pdgId_;
};

LeptonEfficiencyCorrector:: LeptonEfficiencyCorrector(std::vector<std::string> files, std::vector<std::string> histos) {
  effmaps_.clear();
  if(files.size()!=histos.size()) {
    std::cout << "ERROR! There should be one histogram per input file! Returning 0 as SF." << std::endl;
    return;
  }

  for(int i=0; i<(int)files.size();++i) {
    TFile *f = TFile::Open(files[i].c_str(),"read");
    if(!f) {
      std::cout << "WARNING! File " << files[i] << " cannot be opened. Skipping this scale factor " << std::endl;
      continue;
    }
    TH2F *hist = (TH2F*)(f->Get(histos[i].c_str()))->Clone(("eff_"+histos[i]).c_str());
    hist->SetDirectory(0);
    if(!hist) {
      std::cout << "ERROR! Histogram " << histos[i] << " not in file " << files[i] << ". Not considering this SF. " << std::endl;
      continue;
    } else {
      std::cout << "Loading histogram " << histos[i] << " from file " << files[i] << "... " << std::endl;
    }
    effmaps_.push_back(hist);
    f->Close();
  }
}

void LeptonEfficiencyCorrector::setLeptons(int nLep, int *lepPdgId, float *lepPt, float *lepEta) {
  nLep_ = nLep; Lep_pdgId_ = lepPdgId; Lep_pt_ = lepPt; Lep_eta_ = lepEta;
}

float LeptonEfficiencyCorrector::getSF(int pdgid, float pt, float eta) {
  float out=1.;
  //float x = abs(pdgid)==13 ? pt : eta;
  //float y = abs(pdgid)==13 ? fabs(eta) : pt;
  //float x = abs(pdgid)==13 ? fabs(eta) : eta;
  //float y = abs(pdgid)==13 ? pt : pt;
  float x = abs(pdgid)==13 ? pt : eta;
  float y = abs(pdgid)==13 ? fabs(eta) : pt;
  for(std::vector<TH2F*>::iterator hist=effmaps_.begin(); hist<effmaps_.end(); ++hist) {
    WeightCalculatorFromHistogram wc(*hist);
    out *= wc.getWeight(x,y);
  }
  return out;
}

float LeptonEfficiencyCorrector::getSFErr(int pdgid, float pt, float eta) {
  float out=1.;
  //float x = abs(pdgid)==13 ? pt : eta;
  //float y = abs(pdgid)==13 ? fabs(eta) : pt;
  //float x = abs(pdgid)==13 ? fabs(eta) : eta;
  //float y = abs(pdgid)==13 ? pt : pt;
  float x = pt;
  float y = abs(pdgid)==13 ? fabs(eta) : eta;
  for(std::vector<TH2F*>::iterator hist=effmaps_.begin(); hist<effmaps_.end(); ++hist) {
    WeightCalculatorFromHistogram wc(*hist);
    out *= wc.getWeightErr(x,y);
  }
  return out;
}

const std::vector<float> & LeptonEfficiencyCorrector::run() {
  ret_.clear();
  for (int iL = 0, nL = nLep_; iL < nL; ++iL) {
    ret_.push_back(getSF((Lep_pdgId_)[iL], (Lep_pt_)[iL], (Lep_eta_)[iL]));
  }
  return ret_;
}

//string path = "https://ttedesch.web.cern.ch/ttedesch/nanoAOD-tools/python/postprocessing/data/leptonSF/";
string path = "python/postprocessing/data/leptonSF/"

std::vector<std::string> mu_f_2016{remote_storage + path + "Mu_RunBCDEFGH_SF_ID_2016_syst.root"};
std::vector<std::string> el_f_2016{remote_storage + path + "EGM2D_RECO_SF_2016.root", remote_storage + path + "2016LegacyReReco_ElectronMVA90noiso_Fall17V2.root"};
std::vector<std::string> mu_h_2016{"NUM_TightID_DEN_genTracks_eta_pt"};
std::vector<std::string> el_h_2016{"EGamma_SF2D", "EGamma_SF2D"};
std::vector<std::string> mu_f_2017{remote_storage + path + "Muon_RunBCDEF_SF_ID_2017.root", remote_storage + path + "Muon_RunBCDEF_SF_ISO_2017.root"};
std::vector<std::string> el_f_2017{remote_storage + path + "EGM2D_2017_passingRECO_highEt.root", path + "Electron_MVA90_2017.root"};
std::vector<std::string> mu_h_2017{"NUM_TightID_DEN_genTracks_pt_abseta", "NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta"};
std::vector<std::string> el_h_2017{"EGamma_SF2D", "EGamma_SF2D"};
std::vector<std::string> mu_f_2018{remote_storage + path + "Muon_RunBCDEF_SF_ID_2018.root", remote_storage + path + "Muon_RunBCDEF_SF_ISO_2018.root"};
std::vector<std::string> el_f_2018{remote_storage + path + "EGM2D_passingRECO_2018All.root", remote_storage + path + "2018_ElectronMVA90Iso.root"};
std::vector<std::string> mu_h_2018{"NUM_TightID_DEN_TrackerMuons_pt_abseta", "NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta"};
std::vector<std::string> el_h_2018{"EGamma_SF2D", "EGamma_SF2D"};

LeptonEfficiencyCorrector worker_mu_2016(mu_f_2016, mu_h_2016);
LeptonEfficiencyCorrector worker_el_2016(el_f_2016, el_h_2016);
LeptonEfficiencyCorrector worker_mu_2017(mu_f_2017, mu_h_2017);
LeptonEfficiencyCorrector worker_el_2017(el_f_2017, el_h_2017);
LeptonEfficiencyCorrector worker_mu_2018(mu_f_2018, mu_h_2018);
LeptonEfficiencyCorrector worker_el_2018(el_f_2018, el_h_2018);



RVec<float> ElectronSFs(rvec_f Electron_pt, rvec_f Electron_eta, rvec_i Electron_pdgId, string Year){
    
    LeptonEfficiencyCorrector worker_el;

    if (Year == "2016") {
        worker_el = worker_el_2016;
    }
    if (Year == "2017"){
        worker_el = worker_el_2017;
    }
    
    if (Year == "2018"){
        worker_el = worker_el_2018;
    }
    
    RVec<float> sf_el(Electron_pt.size());
    RVec<float> sferr_el(Electron_pt.size());
    
    for (size_t j = 0; j < Electron_pt.size(); j++) sf_el[j] =  worker_el.getSF(Electron_pdgId[j], Electron_pt[j], Electron_eta[j]);
    for (size_t j = 0; j < Electron_pt.size(); j++) sferr_el[j] =  worker_el.getSFErr(Electron_pdgId[j], Electron_pt[j], Electron_eta[j]);
    
    RVec<float> Electron_effSF_errUp(Electron_pt.size());
    for (size_t j = 0; j < Electron_pt.size(); j++) Electron_effSF_errUp[j] = sferr_el[j] + sf_el[j];
    
    RVec<float> Electron_effSF_errDown(Electron_pt.size());
    for (size_t j = 0; j < Electron_pt.size(); j++) Electron_effSF_errDown[j] = sferr_el[j] - sf_el[j];
    
    RVec<float> result;
    for (int j = 0; j < Electron_pt.size(); j++){
        result.emplace_back(sf_el[j]);
        result.emplace_back(Electron_effSF_errUp[j]);
        result.emplace_back(Electron_effSF_errDown[j]);
    }
    
    return result;
}

RVec<float> MuonSFs(rvec_f Muon_pt, rvec_f Muon_eta, rvec_i Muon_pdgId, string Year){
    
    LeptonEfficiencyCorrector worker_mu;

    if (Year == "2016") {
        worker_mu = worker_mu_2016;
    }
    
    if (Year == "2017"){
        worker_mu = worker_mu_2017;
    }
    
    if (Year == "2018"){
        worker_mu = worker_mu_2018;
    }
    
    RVec<float> sf_mu(Muon_pt.size());
    RVec<float> sferr_mu(Muon_pt.size()); 
    
    for (size_t j = 0; j < Muon_pt.size(); j++) sf_mu[j] =  worker_mu.getSF(Muon_pdgId[j], Muon_pt[j], Muon_eta[j]);
    for (size_t j = 0; j < Muon_pt.size(); j++) sferr_mu[j] =  worker_mu.getSFErr(Muon_pdgId[j], Muon_pt[j], Muon_eta[j]);
    
    RVec<float> Muon_effSF_errUp(Muon_pt.size());
    for (size_t j = 0; j < Muon_pt.size(); j++) Muon_effSF_errUp[j] = sferr_mu[j] + sf_mu[j];
    
    RVec<float> Muon_effSF_errDown(Muon_pt.size());
    for (size_t j = 0; j < Muon_pt.size(); j++) Muon_effSF_errDown[j] = sferr_mu[j] - sf_mu[j];
    
    RVec<float> result;
    for (int j = 0; j < Muon_pt.size(); j++){
        result.emplace_back(sf_mu[j]);
        result.emplace_back(Muon_effSF_errUp[j]);
        result.emplace_back(Muon_effSF_errDown[j]);
    }
    
    return result;
}

//string remote_storage = "https://ttedesch.web.cern.ch/ttedesch/nanoAOD-tools/python/postprocessing/data/";
string path_pu = "python/postprocessing/data/pileup/";

void set_null_directory(TH2F *histo){
    histo->SetDirectory(NULL);
}


TFile *pufile_data2016 = TFile::Open(TString(remote_storage) + TString(path_pu) + TString("PileupData_GoldenJSON_Full2016.root"));
TH1 *histo_target_2016 = (TH1*)pufile_data2016->Get("pileup");
TH1 *histo_target_2016_plus = (TH1*)pufile_data2016->Get("pileup_plus");
TH1 *histo_target_2016_minus = (TH1*)pufile_data2016->Get("pileup_minus");
TFile *pufile_mc2016 = TFile::Open(TString(remote_storage) + TString(path_pu) + TString("pileup_profile_Summer16.root"));
TH1 *histo_2016 = (TH1*)pufile_mc2016->Get("pu_mc");

TFile *pufile_data2017 = TFile::Open(TString(remote_storage) + TString(path_pu) + TString("PileupHistogram-goldenJSON-13tev-2017-99bins_withVar.root"));
TH1 *histo_target_2017 = (TH1*)pufile_data2017->Get("pileup");
TH1 *histo_target_2017_plus = (TH1*)pufile_data2017->Get("pileup_plus");
TH1 *histo_target_2017_minus = (TH1*)pufile_data2017->Get("pileup_minus");
TFile *pufile_mc2017 = TFile::Open(TString(remote_storage) + TString(path_pu) + TString("pileup/mcPileup2017.root"));
TH1 *histo_2017 = (TH1*)pufile_mc2017->Get("pu_mc");

TFile *pufile_data2018 = TFile::Open(TString(remote_storage) + TString(path_pu) + TString("PileupHistogram-goldenJSON-13tev-2018-100bins_withVar.root"));
TH1 *histo_target_2018 = (TH1*)pufile_data2018->Get("pileup");
TH1 *histo_target_2018_plus = (TH1*)pufile_data2018->Get("pileup_plus");
TH1 *histo_target_2018_minus = (TH1*)pufile_data2018->Get("pileup_minus");
TFile *pufile_mc2018 = TFile::Open(TString(remote_storage) + TString(path_pu) + TString("mcPileup2018.root"));
TH1 *histo_2018 = (TH1*)pufile_mc2018->Get("pu_mc");
//close files??

RVec<float> puWeight(string Year, int Pileup_nTrueInt){
    RVec<float> result(3);
    if (Year == "2016"){ 
        WeightCalculatorFromHistogram worker_2016(histo_2016, histo_target_2016, true, true, false);
        WeightCalculatorFromHistogram worker_2016_plus(histo_2016, histo_target_2016_plus, true, true, false);
        WeightCalculatorFromHistogram worker_2016_minus(histo_2016, histo_target_2016_minus, true, true, false);
        if(Pileup_nTrueInt < histo_2016->GetNbinsX()) result[0] = worker_2016.getWeight(Pileup_nTrueInt);
        else result[0] = 1;
        if(Pileup_nTrueInt < histo_2016->GetNbinsX()) result[1] = worker_2016_plus.getWeight(Pileup_nTrueInt);
        else result[1] = 1;
        if(Pileup_nTrueInt < histo_2016->GetNbinsX()) result[2] = worker_2016_minus.getWeight(Pileup_nTrueInt);
        else result[2] = 1;
    }
    else if (Year == "2017"){
        WeightCalculatorFromHistogram worker_2017(histo_2017, histo_target_2017, true, true, false);
        WeightCalculatorFromHistogram worker_2017_plus(histo_2017, histo_target_2017_plus, true, true, false);
        WeightCalculatorFromHistogram worker_2017_minus(histo_2017, histo_target_2017_minus, true, true, false);
        if(Pileup_nTrueInt < histo_2017->GetNbinsX()) result[0] = worker_2017.getWeight(Pileup_nTrueInt);
        else result[0] = 1;
        if(Pileup_nTrueInt < histo_2017->GetNbinsX()) result[1] = worker_2017_plus.getWeight(Pileup_nTrueInt);
        else result[1] = 1;
        if(Pileup_nTrueInt < histo_2017->GetNbinsX()) result[2] = worker_2017_minus.getWeight(Pileup_nTrueInt);
        else result[2] = 1;
    }
    else if (Year == "2018"){
        WeightCalculatorFromHistogram worker_2018(histo_2018, histo_target_2018, true, true, false);
        WeightCalculatorFromHistogram worker_2018_plus(histo_2018, histo_target_2018_plus, true, true, false);
        WeightCalculatorFromHistogram worker_2018_minus(histo_2018, histo_target_2018_minus, true, true, false);
        if(Pileup_nTrueInt < histo_2018->GetNbinsX()) result[0] = worker_2018.getWeight(Pileup_nTrueInt);
        else result[0] = 1;
        if(Pileup_nTrueInt < histo_2018->GetNbinsX()) result[1] = worker_2018_plus.getWeight(Pileup_nTrueInt);
        else result[1] = 1;
        if(Pileup_nTrueInt < histo_2018->GetNbinsX()) result[2] = worker_2018_minus.getWeight(Pileup_nTrueInt);
        else result[2] = 1;
    }
    return result;    
}






//string remote_storage = "https://ttedesch.web.cern.ch/ttedesch/nanoAOD-tools/";
//string remote_storage_ = "https://ttedesch.web.cern.ch/ttedesch/nanoAOD-tools/";
string path_pf =  "data/prefire_maps/";

TFile *L1PrefiringMaps = TFile::Open(TString(remote_storage) + TString(path_pf) + TString("L1PrefiringMaps.root"));
TH2F * L1prefiring_jetptvseta_2017BtoF = (TH2F *) L1PrefiringMaps->Get("L1prefiring_jetptvseta_2017BtoF");
TH2F * L1prefiring_photonptvseta_2017BtoF = (TH2F *) L1PrefiringMaps->Get("L1prefiring_photonptvseta_2017BtoF");


//TH2F * L1prefiring_jetptvseta_UL2017BtoF = (TH2F *) L1PrefiringMaps->Get("L1prefiring_jetptvseta_UL2017BtoF");
//TH2F * L1prefiring_photonptvseta_UL2017BtoF = (TH2F *) L1PrefiringMaps->Get("L1prefiring_photonptvseta_UL2017BtoF");

float GetPrefireProbability(TH2F *Map, float eta, float pt, float maxpt, int variation){
        float x = maxpt - 0.01;
        int bin = Map->FindBin(eta, min(pt, x));
        float pref_prob = Map->GetBinContent(bin);

        float stat = Map->GetBinError(bin);  //bin statistical uncertainty
        float syst = 0.2 * pref_prob;  //20% of prefire rate

        float r = sqrt(stat * stat + syst * syst);
        float y = pref_prob + r;
        float z = pref_prob - r;
        if (variation == 1) pref_prob = min(y , float(1.0));
        else if (variation == -1) pref_prob = max(z , float(0.0));

        return pref_prob;
}

float EGvalue(rvec_f Photon_pt, rvec_f Photon_eta, rvec_i Photon_jetIdx, rvec_i Photon_electronIdx, rvec_f Electron_pt, rvec_f Electron_eta, rvec_i Electron_jetIdx, rvec_i Electron_photonIdx, int jid, int s){
        float phopf = 1.0;
        vector<int> PhotonInJet;
        float JetMinPt = 20;  //Min/Max Values may need to be fixed for new maps
        float JetMaxPt = 500;
        float JetMinEta = 2.0;
        float JetMaxEta = 3.0;
        float PhotonMinPt = 20;
        float PhotonMaxPt = 500;
        float PhotonMinEta = 2.0;
        float PhotonMaxEta = 3.0;

        TH2F *photon_map = L1prefiring_photonptvseta_2017BtoF;

        for(int i = 0; i < Photon_pt.size(); i++){
            if (Photon_jetIdx[i] == jid){
                if (Photon_pt[i] >= PhotonMinPt && abs(Photon_eta[i]) <= PhotonMaxEta && abs(Photon_eta[i]) >= PhotonMinEta){
                    float phopf_temp = 1 - GetPrefireProbability(photon_map, Photon_eta[i], Photon_pt[i], PhotonMaxPt, s);
                    float elepf_temp = 1.0;
                    if (Photon_electronIdx[i] > -1){
                        if (Electron_pt[Photon_electronIdx[i]] >= PhotonMinPt && abs(Electron_eta[Photon_electronIdx[i]]) <= PhotonMaxEta and abs(Electron_eta[Photon_electronIdx[i]]) >= PhotonMinEta) elepf_temp = 1 - GetPrefireProbability(photon_map, Electron_eta[Photon_electronIdx[i]], Electron_pt[Photon_electronIdx[i]], PhotonMaxPt, s);
                    }
                    phopf = phopf * min(phopf_temp, elepf_temp);
                    PhotonInJet.push_back(i);
                }
            }
        }
        for(int i = 0; i < Electron_pt.size(); i++){
            if (Electron_jetIdx[i] == jid && find(PhotonInJet.begin(), PhotonInJet.end(), Electron_photonIdx[i]) == PhotonInJet.end()){
                if (Electron_pt[i] >= PhotonMinPt && abs(Electron_eta[i]) <= PhotonMaxEta && abs(Electron_eta[i]) >= PhotonMinEta) phopf = phopf * ( 1 - GetPrefireProbability(photon_map, Electron_eta[i], Electron_pt[i], PhotonMaxPt, s));
            }
        }
        return phopf;
}

RVec<float> PrefCorr(rvec_f Photon_pt, rvec_f Photon_eta, rvec_i Photon_jetIdx, rvec_i Photon_electronIdx, rvec_f Electron_pt, rvec_f Electron_eta, rvec_i Electron_jetIdx, rvec_i Electron_photonIdx, rvec_f Jet_pt, rvec_f Jet_eta){
    RVec<float> result(3);

    float JetMinPt = 20;  //Min/Max Values may need to be fixed for new maps
    float JetMaxPt = 500;
    float JetMinEta = 2.0;
    float JetMaxEta = 3.0;
    float PhotonMinPt = 20;
    float PhotonMaxPt = 500;
    float PhotonMinEta = 2.0;
    float PhotonMaxEta = 3.0;
    TH2F *jet_map =  L1prefiring_jetptvseta_2017BtoF;

    float prefw = 1.0;

    vector<int> v{0,1,-1};
    for (int j = 0; j < v.size(); j++){
        prefw = 1.0; // new
        int s = v[j];
        for(int i = 0; i < Jet_pt.size(); i++){
            float jetpf = 1.0;
            //PhotonInJet = []

            if(Jet_pt[i] >= JetMinPt && abs(Jet_eta[i]) <= JetMaxEta && abs(Jet_eta[i]) >= JetMinEta){
                jetpf = jetpf * ( 1 - GetPrefireProbability(jet_map, Jet_eta[i], Jet_pt[i], JetMaxPt, s));
                float phopf = EGvalue(Photon_pt, Photon_eta, Photon_jetIdx, Photon_electronIdx, Electron_pt, Electron_eta, Electron_jetIdx, Electron_photonIdx, i, s);
                prefw = prefw * min(jetpf, phopf);
            }
        }
        //Then loop over all photons/electrons not associated to jets
        prefw = prefw * EGvalue(Photon_pt, Photon_eta, Photon_jetIdx, Photon_electronIdx, Electron_pt, Electron_eta, Electron_jetIdx, Electron_photonIdx, -1, s);
        result[j] = prefw;
    }
    return result;

}

/*
def mk_safe(fct, *args):
    try:
        return fct(*args)
    except Exception as e:
        if any('Error in function boost::math::erf_inv' in arg for arg in e.args):
            print('WARNING: catching exception and returning -1. Exception arguments: %s' % e.args)
            return -1.
        else:
            raise e
*/


//RoccoR roccor_2016("RoccoR2016.txt");
//RoccoR roccor_2017("RoccoR2017.txt");
//RoccoR roccor_2018("RoccoR2018.txt");
RoccoR roccor_2016aUL("RoccoR2016aUL.txt");
RoccoR roccor_2016bUL("RoccoR2016bUL.txt");
RoccoR roccor_2017UL("RoccoR2017UL.txt");
RoccoR roccor_2018UL("RoccoR2018UL.txt");

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0, 1);//uniform distribution between 0 and 1


RVec<float> muonScaleRes(rvec_f Muon_pt, rvec_f Muon_eta, rvec_f Muon_phi, rvec_i Muon_charge, rvec_i Muon_nTrackerLayers, rvec_i Muon_genPartIdx, rvec_f GenPart_pt, string era){
    RVec<float> result;
    RVec<float> pt_corr, pt_err;
    RoccoR roccor;
    
    if (era == "2016") roccor = roccor_2016;
    else if (era == "2017") roccor = roccor_2017;
    else if (era == "2018") roccor = roccor_2018;
    //else if (era == "2016aUL") roccor = roccor_2016aUL;
    //else if (era == "2016bUL") roccor = roccor_2016bUL;
    //else if (era == "2017UL") roccor = roccor_2017UL;
    //else roccor = roccor_2018UL;
    bool isMC = true;
    
    if (isMC == true){
        for (int j = 0; j < Muon_pt.size(); j++){
            int genIdx = Muon_genPartIdx[j];
            if(genIdx >= 0 && genIdx < GenPart_pt.size()){
                //pt_corr.emplace_back(Muon_pt[j] * mk_safe(roccor.kSpreadMC, Muon_charge[j], Muon_pt[j], Muon_eta[j], Muon_phi[j], GenPart_pt[genIdx]));
                pt_corr.emplace_back(Muon_pt[j] * roccor.kSpreadMC(Muon_charge[j], Muon_pt[j], Muon_eta[j], Muon_phi[j], GenPart_pt[genIdx]));
                //pt_err.emplace_back(Muon_pt[j] * mk_safe(roccor.kSpreadMCerror, Muon_charge[j], Muon_pt[j], Muon_eta[j], Muon_phi[j], GenPart_pt[genIdx]));
                pt_err.emplace_back(Muon_pt[j] * roccor.kSpreadMCerror(Muon_charge[j], Muon_pt[j], Muon_eta[j], Muon_phi[j], GenPart_pt[genIdx]));
            }
            else{
                float u1 = dis(gen);
                //pt_corr.emplace_back(Muon_pt[j] * mk_safe(roccor.kSmearMC, Muon_charge[j], Muon_pt[j], Muon_eta[j], Muon_phi[j], Muon_nTrackerLayers[j], u1));
                pt_corr.emplace_back(Muon_pt[j] * roccor.kSmearMC(Muon_charge[j], Muon_pt[j], Muon_eta[j], Muon_phi[j], Muon_nTrackerLayers[j], u1));
                //pt_err.emplace_back(Muon_pt[j] * mk_safe(roccor.kSmearMCerror, Muon_charge[j], Muon_pt[j], Muon_eta[j], Muon_phi[j], Muon_nTrackerLayers[j], u1));
                pt_err.emplace_back(Muon_pt[j] * roccor.kSmearMCerror(Muon_charge[j], Muon_pt[j], Muon_eta[j], Muon_phi[j], Muon_nTrackerLayers[j], u1));
            }
        }
    }
    else{
        for (int j = 0; j < Muon_pt.size(); j++){
            //pt_corr.emplace_back(Muon_pt[j] * mk_safe(roccor.kScaleDT, Muon_charge[j], Muon_pt[j], Muon_eta[j], Muon_phi[j]));
            pt_corr.emplace_back(Muon_pt[j] * roccor.kScaleDT(Muon_charge[j], Muon_pt[j], Muon_eta[j], Muon_phi[j]));
            //pt_err.emplace_back(Muon_pt[j] * mk_safe(roccor.kScaleDTerror, Muon_charge[j], Muon_pt[j], Muonu_eta[j], Muon_phi[j]));
            pt_err.emplace_back(Muon_pt[j] * roccor.kScaleDTerror(Muon_charge[j], Muon_pt[j], Muon_eta[j], Muon_phi[j]));
        }
    }
    
    RVec<float> pt_corr_up, pt_corr_down;
    
    for (int j = 0; j < Muon_pt.size(); j++){
        pt_corr_up.emplace_back(max(pt_corr[j] + pt_err[j], float(0.0)));
        pt_corr_down.emplace_back(max(pt_corr[j] - pt_err[j], float(0.0)));
    }
    
    for (int j = 0; j < Muon_pt.size(); j++){
        result.emplace_back(pt_corr[j]);
        result.emplace_back(pt_corr_up[j]);
        result.emplace_back(pt_corr_down[j]);
    }
    
    return result;
}

RVec<float> muonScaleRes_data(rvec_f Muon_pt, rvec_f Muon_eta, rvec_f Muon_phi, rvec_i Muon_charge, string era){
    RVec<float> result;
    RVec<float> pt_corr, pt_err;
    RoccoR roccor;
    
    if (era == "2016") roccor = roccor_2016;
    else if (era == "2017") roccor = roccor_2017;
    else if (era == "2018") roccor = roccor_2018;
    //else if (era == "2016aUL") roccor = roccor_2016aUL;
    //else if (era == "2016bUL") roccor = roccor_2016bUL;
    //else if (era == "2017UL") roccor = roccor_2017UL;
    //else roccor = roccor_2018UL;
    
    for (int j = 0; j < Muon_pt.size(); j++){
            //pt_corr.emplace_back(Muon_pt[j] * mk_safe(roccor.kScaleDT, Muon_charge[j], Muon_pt[j], Muon_eta[j], Muon_phi[j]));
            pt_corr.emplace_back(Muon_pt[j] * roccor.kScaleDT(Muon_charge[j], Muon_pt[j], Muon_eta[j], Muon_phi[j]));
            //pt_err.emplace_back(Muon_pt[j] * mk_safe(roccor.kScaleDTerror, Muon_charge[j], Muon_pt[j], Muonu_eta[j], Muon_phi[j]));
            pt_err.emplace_back(Muon_pt[j] * roccor.kScaleDTerror(Muon_charge[j], Muon_pt[j], Muon_eta[j], Muon_phi[j]));
    }
    
    RVec<float> pt_corr_up, pt_corr_down;
    
    for (int j = 0; j < Muon_pt.size(); j++){
        pt_corr_up.emplace_back(max(pt_corr[j] + pt_err[j], float(0.0)));
        pt_corr_down.emplace_back(max(pt_corr[j] - pt_err[j], float(0.0)));
    }
    
    for (int j = 0; j < Muon_pt.size(); j++){
        result.emplace_back(pt_corr[j]);
        result.emplace_back(pt_corr_up[j]);
        result.emplace_back(pt_corr_down[j]);
    }
    
    return result;
}

RVec<float> getMatrixColumn(const RVec<RVec<float>> & matrix, int column_index){
    RVec<float> result;
    for (int i = 0; i < matrix.size(); i++) result.emplace_back(matrix[i][column_index]);
    return result;
}

RVec<float> getFlattenedMatrixColumn(rvec_f flattened_matrix, int nColumns, int column_index){
    RVec<float> result;
    for (int i = 0; i < flattened_matrix.size()/nColumns; i++) result.emplace_back(flattened_matrix[column_index + i*nColumns]);
    return result;
}

float htProducer(int nJet, rvec_f Jet_pt){
    float ht(0.0);
    for (unsigned i=0; i<nJet; i++){
      ht += Jet_pt[i];
    }
    return ht;
}



///////////////////////////////////////////////////////////////////////////////////////////////BTAG utils//////////////////////////////////////////////777
#ifndef BTagEntry_H
#define BTagEntry_H

/**
 *
 * BTagEntry
 *
 * Represents one pt- or discriminator-dependent calibration function.
 *
 * measurement_type:    e.g. comb, ttbar, di-mu, boosted, ...
 * sys_type:            e.g. central, plus, minus, plus_JEC, plus_JER, ...
 *
 * Everything is converted into a function, as it is easiest to store it in a
 * txt or json file.
 *
 ************************************************************/

#include <string>
#include <TF1.h>
#include <TH1.h>


class BTagEntry
{
public:
  enum OperatingPoint {
    OP_LOOSE=0,
    OP_MEDIUM=1,
    OP_TIGHT=2,
    OP_RESHAPING=3,
  };
  enum JetFlavor {
    FLAV_B=0,
    FLAV_C=1,
    FLAV_UDSG=2,
  };
  struct Parameters {
    OperatingPoint operatingPoint;
    std::string measurementType;
    std::string sysType;
    JetFlavor jetFlavor;
    float etaMin;
    float etaMax;
    float ptMin;
    float ptMax;
    float discrMin;
    float discrMax;

    // default constructor
    Parameters(
      OperatingPoint op=OP_TIGHT,
      std::string measurement_type="comb",
      std::string sys_type="central",
      JetFlavor jf=FLAV_B,
      float eta_min=-99999.,
      float eta_max=99999.,
      float pt_min=0.,
      float pt_max=99999.,
      float discr_min=0.,
      float discr_max=99999.
    );

  };

  BTagEntry() {}
  BTagEntry(const std::string &csvLine);
  BTagEntry(const std::string &func, Parameters p);
  BTagEntry(const TF1* func, Parameters p);
  BTagEntry(const TH1* histo, Parameters p);
  ~BTagEntry() {}
  static std::string makeCSVHeader();
  std::string makeCSVLine() const;
  static std::string trimStr(std::string str);

  // public, no getters needed
  std::string formula;
  Parameters params;

};

#endif  // BTagEntry_H


#ifndef BTagCalibration_H
#define BTagCalibration_H

/**
 * BTagCalibration
 *
 * The 'hierarchy' of stored information is this:
 * - by tagger (BTagCalibration)
 *   - by operating point or reshape bin
 *     - by jet parton flavor
 *       - by type of measurement
 *         - by systematic
 *           - by eta bin
 *             - as 1D-function dependent of pt or discriminant
 *
 ************************************************************/

#include <map>
#include <vector>
#include <string>
#include <istream>
#include <ostream>


class BTagCalibration
{
public:
  BTagCalibration() {}
  BTagCalibration(const std::string &tagger);
  BTagCalibration(const std::string &tagger, const std::string &filename);
  ~BTagCalibration() {}

  std::string tagger() const {return tagger_;}

  void addEntry(const BTagEntry &entry);
  const std::vector<BTagEntry>& getEntries(const BTagEntry::Parameters &par) const;

  void readCSV(std::istream &s);
  void readCSV(const std::string &s);
  void makeCSV(std::ostream &s) const;
  std::string makeCSV() const;

protected:
  static std::string token(const BTagEntry::Parameters &par);

  std::string tagger_;
  std::map<std::string, std::vector<BTagEntry> > data_;

};

#endif  // BTagCalibration_H


#ifndef BTagCalibrationReader_H
#define BTagCalibrationReader_H

/**
 * BTagCalibrationReader
 *
 * Helper class to pull out a specific set of BTagEntry's out of a
 * BTagCalibration. TF1 functions are set up at initialization time.
 *
 ************************************************************/

#include <memory>
#include <string>



class BTagCalibrationReader
{
public:
  class BTagCalibrationReaderImpl;

  BTagCalibrationReader() {}
  BTagCalibrationReader(BTagEntry::OperatingPoint op,
                        const std::string & sysType="central",
                        const std::vector<std::string> & otherSysTypes={});

  void load(const BTagCalibration & c,
            BTagEntry::JetFlavor jf,
            const std::string & measurementType="comb");

  double eval(BTagEntry::JetFlavor jf,
              float eta,
              float pt,
              float discr=0.) const;

  double eval_auto_bounds(const std::string & sys,
                          BTagEntry::JetFlavor jf,
                          float eta,
                          float pt,
                          float discr=0.) const;

  std::pair<float, float> min_max_pt(BTagEntry::JetFlavor jf,
                                     float eta,
                                     float discr=0.) const;
protected:
  std::shared_ptr<BTagCalibrationReaderImpl> pimpl;
};


#endif  // BTagCalibrationReader_H

#include <iostream>
#include <exception>
#include <algorithm>
#include <sstream>


BTagEntry::Parameters::Parameters(
  OperatingPoint op,
  std::string measurement_type,
  std::string sys_type,
  JetFlavor jf,
  float eta_min,
  float eta_max,
  float pt_min,
  float pt_max,
  float discr_min,
  float discr_max
):
  operatingPoint(op),
  measurementType(measurement_type),
  sysType(sys_type),
  jetFlavor(jf),
  etaMin(eta_min),
  etaMax(eta_max),
  ptMin(pt_min),
  ptMax(pt_max),
  discrMin(discr_min),
  discrMax(discr_max)
{
  std::transform(measurementType.begin(), measurementType.end(),
                 measurementType.begin(), ::tolower);
  std::transform(sysType.begin(), sysType.end(),
                 sysType.begin(), ::tolower);
}

BTagEntry::BTagEntry(const std::string &csvLine)
{
  // make tokens
  std::stringstream buff(csvLine);
  std::vector<std::string> vec;
  std::string token;
  while (std::getline(buff, token, ","[0])) {
    token = BTagEntry::trimStr(token);
    if (token.empty()) {
      continue;
    }
    vec.push_back(token);
  }
  if (vec.size() != 11) {
std::cerr << "ERROR in BTagCalibration: "
          << "Invalid csv line; num tokens != 11: "
          << csvLine;
throw std::exception();
  }

  // clean string values
  char chars[] = " \"\n";
  for (unsigned int i = 0; i < strlen(chars); ++i) {
    vec[1].erase(remove(vec[1].begin(),vec[1].end(),chars[i]),vec[1].end());
    vec[2].erase(remove(vec[2].begin(),vec[2].end(),chars[i]),vec[2].end());
    vec[10].erase(remove(vec[10].begin(),vec[10].end(),chars[i]),vec[10].end());
  }

  // make formula
  formula = vec[10];
  TF1 f1("", formula.c_str());  // compile formula to check validity
  if (f1.IsZombie()) {
std::cerr << "ERROR in BTagCalibration: "
          << "Invalid csv line; formula does not compile: "
          << csvLine;
throw std::exception();
  }

  // make parameters
  unsigned op = stoi(vec[0]);
  if (op > 3) {
std::cerr << "ERROR in BTagCalibration: "
          << "Invalid csv line; OperatingPoint > 3: "
          << csvLine;
throw std::exception();
  }
  unsigned jf = stoi(vec[3]);
  if (jf > 2) {
std::cerr << "ERROR in BTagCalibration: "
          << "Invalid csv line; JetFlavor > 2: "
          << csvLine;
throw std::exception();
  }
  params = BTagEntry::Parameters(
    BTagEntry::OperatingPoint(op),
    vec[1],
    vec[2],
    BTagEntry::JetFlavor(jf),
    stof(vec[4]),
    stof(vec[5]),
    stof(vec[6]),
    stof(vec[7]),
    stof(vec[8]),
    stof(vec[9])
  );
}

BTagEntry::BTagEntry(const std::string &func, BTagEntry::Parameters p):
  formula(func),
  params(p)
{
  TF1 f1("", formula.c_str());  // compile formula to check validity
  if (f1.IsZombie()) {
std::cerr << "ERROR in BTagCalibration: "
          << "Invalid func string; formula does not compile: "
          << func;
throw std::exception();
  }
}

BTagEntry::BTagEntry(const TF1* func, BTagEntry::Parameters p):
  formula(std::string(func->GetExpFormula("p").Data())),
  params(p)
{
  if (func->IsZombie()) {
std::cerr << "ERROR in BTagCalibration: "
          << "Invalid TF1 function; function is zombie: "
          << func->GetName();
throw std::exception();
  }
}

// Creates chained step functions like this:
// "<prevous_bin> : x<bin_high_bound ? bin_value : <next_bin>"
// e.g. "x<0 ? 1 : x<1 ? 2 : x<2 ? 3 : 4"
std::string th1ToFormulaLin(const TH1* hist) {
  int nbins = hist->GetNbinsX();
  TAxis const* axis = hist->GetXaxis();
  std::stringstream buff;
  buff << "x<" << axis->GetBinLowEdge(1) << " ? 0. : ";  // default value
  for (int i=1; i<nbins+1; ++i) {
    char tmp_buff[50];
    sprintf(tmp_buff,
            "x<%g ? %g : ",  // %g is the smaller one of %e or %f
            axis->GetBinUpEdge(i),
            hist->GetBinContent(i));
    buff << tmp_buff;
  }
  buff << 0.;  // default value
  return buff.str();
}

// Creates step functions making a binary search tree:
// "x<mid_bin_bound ? (<left side tree>) : (<right side tree>)"
// e.g. "x<2 ? (x<1 ? (x<0 ? 0:0.1) : (1)) : (x<4 ? (x<3 ? 2:3) : (0))"
std::string th1ToFormulaBinTree(const TH1* hist, int start=0, int end=-1) {
  if (end == -1) {                      // initialize
    start = 0.;
    end = hist->GetNbinsX()+1;
    TH1* h2 = (TH1*) hist->Clone();
    h2->SetBinContent(start, 0);  // kill underflow
    h2->SetBinContent(end, 0);    // kill overflow
    std::string res = th1ToFormulaBinTree(h2, start, end);
    delete h2;
    return res;
  }
  if (start == end) {                   // leave is reached
    char tmp_buff[20];
    sprintf(tmp_buff, "%g", hist->GetBinContent(start));
    return std::string(tmp_buff);
  }
  if (start == end - 1) {               // no parenthesis for neighbors
    char tmp_buff[70];
    sprintf(tmp_buff,
            "x<%g ? %g:%g",
            hist->GetXaxis()->GetBinUpEdge(start),
            hist->GetBinContent(start),
            hist->GetBinContent(end));
    return std::string(tmp_buff);
  }

  // top-down recursion
  std::stringstream buff;
  int mid = (end-start)/2 + start;
  char tmp_buff[25];
  sprintf(tmp_buff,
          "x<%g ? (",
          hist->GetXaxis()->GetBinUpEdge(mid));
  buff << tmp_buff
       << th1ToFormulaBinTree(hist, start, mid)
       << ") : ("
       << th1ToFormulaBinTree(hist, mid+1, end)
       << ")";
  return buff.str();
}

BTagEntry::BTagEntry(const TH1* hist, BTagEntry::Parameters p):
  params(p)
{
  int nbins = hist->GetNbinsX();
  TAxis const* axis = hist->GetXaxis();

  // overwrite bounds with histo values
  if (params.operatingPoint == BTagEntry::OP_RESHAPING) {
    params.discrMin = axis->GetBinLowEdge(1);
    params.discrMax = axis->GetBinUpEdge(nbins);
  } else {
    params.ptMin = axis->GetBinLowEdge(1);
    params.ptMax = axis->GetBinUpEdge(nbins);
  }

  // balanced full binary tree height = ceil(log(2*n_leaves)/log(2))
  // breakes even around 10, but lower values are more propable in pt-spectrum
  if (nbins < 15) {
    formula = th1ToFormulaLin(hist);
  } else {
    formula = th1ToFormulaBinTree(hist);
  }

  // compile formula to check validity
  TF1 f1("", formula.c_str());
  if (f1.IsZombie()) {
std::cerr << "ERROR in BTagCalibration: "
          << "Invalid histogram; formula does not compile (>150 bins?): "
          << hist->GetName();
throw std::exception();
  }
}

std::string BTagEntry::makeCSVHeader()
{
  return "OperatingPoint, "
         "measurementType, "
         "sysType, "
         "jetFlavor, "
         "etaMin, "
         "etaMax, "
         "ptMin, "
         "ptMax, "
         "discrMin, "
         "discrMax, "
         "formula \n";
}

std::string BTagEntry::makeCSVLine() const
{
  std::stringstream buff;
  buff << params.operatingPoint
       << ", " << params.measurementType
       << ", " << params.sysType
       << ", " << params.jetFlavor
       << ", " << params.etaMin
       << ", " << params.etaMax
       << ", " << params.ptMin
       << ", " << params.ptMax
       << ", " << params.discrMin
       << ", " << params.discrMax
       << ", \"" << formula
       << "\" \n";
  return buff.str();
}

std::string BTagEntry::trimStr(std::string str) {
  size_t s = str.find_first_not_of(" \n\r\t");
  size_t e = str.find_last_not_of (" \n\r\t");

  if((std::string::npos == s) || (std::string::npos == e))
    return "";
  else
    return str.substr(s, e-s+1);
}


#include <fstream>
#include <sstream>



BTagCalibration::BTagCalibration(const std::string &taggr):
  tagger_(taggr)
{}

BTagCalibration::BTagCalibration(const std::string &taggr,
                                 const std::string &filename):
  tagger_(taggr)
{
  std::ifstream ifs(filename);
  if (!ifs.good()) {
std::cerr << "ERROR in BTagCalibration: "
          << "input file not available: "
          << filename;
throw std::exception();
  }
  readCSV(ifs);
  ifs.close();
}

void BTagCalibration::addEntry(const BTagEntry &entry)
{
  data_[token(entry.params)].push_back(entry);
}

const std::vector<BTagEntry>& BTagCalibration::getEntries(
  const BTagEntry::Parameters &par) const
{
  std::string tok = token(par);
  if (!data_.count(tok)) {
std::cerr << "ERROR in BTagCalibration: "
          << "(OperatingPoint, measurementType, sysType) not available: "
          << tok;
throw std::exception();
  }
  return data_.at(tok);
}

void BTagCalibration::readCSV(const std::string &s)
{
  std::stringstream buff(s);
  readCSV(buff);
}

void BTagCalibration::readCSV(std::istream &s)
{
  std::string line;

  // firstline might be the header
  getline(s,line);
  if (line.find("OperatingPoint") == std::string::npos) {
    addEntry(BTagEntry(line));
  }

  //while (getline(s,line)) {
  int line_counter = 0;
  int max_lines = 1000;
  while (getline(s,line) && line_counter < max_lines) {
    line_counter++;
    line = BTagEntry::trimStr(line);
    if (line.empty()) {  // skip empty lines
      continue;
    }
    addEntry(BTagEntry(line));
  }
}

void BTagCalibration::makeCSV(std::ostream &s) const
{
  s << tagger_ << ";" << BTagEntry::makeCSVHeader();
  for (std::map<std::string, std::vector<BTagEntry> >::const_iterator i
           = data_.cbegin(); i != data_.cend(); ++i) {
    const std::vector<BTagEntry> &vec = i->second;
    for (std::vector<BTagEntry>::const_iterator j
             = vec.cbegin(); j != vec.cend(); ++j) {
      s << j->makeCSVLine();
    }
  }
}

std::string BTagCalibration::makeCSV() const
{
  std::stringstream buff;
  makeCSV(buff);
  return buff.str();
}

std::string BTagCalibration::token(const BTagEntry::Parameters &par)
{
  std::stringstream buff;
  buff << par.operatingPoint << ", "
       << par.measurementType << ", "
       << par.sysType;
  return buff.str();
}




class BTagCalibrationReader::BTagCalibrationReaderImpl
{
  friend class BTagCalibrationReader;

public:
  struct TmpEntry {
    float etaMin;
    float etaMax;
    float ptMin;
    float ptMax;
    float discrMin;
    float discrMax;
    TF1 func;
  };

private:
  BTagCalibrationReaderImpl(BTagEntry::OperatingPoint op,
                            const std::string & sysType,
                            const std::vector<std::string> & otherSysTypes={});

  void load(const BTagCalibration & c,
            BTagEntry::JetFlavor jf,
            std::string measurementType);

  double eval(BTagEntry::JetFlavor jf,
              float eta,
              float pt,
              float discr) const;

  double eval_auto_bounds(const std::string & sys,
                          BTagEntry::JetFlavor jf,
                          float eta,
                          float pt,
                          float discr) const;

  std::pair<float, float> min_max_pt(BTagEntry::JetFlavor jf,
                                     float eta,
                                     float discr) const;
 
  std::pair<float, float> min_max_eta(BTagEntry::JetFlavor jf,
                                     float discr) const;

  BTagEntry::OperatingPoint op_;
  std::string sysType_;
  std::vector<std::vector<TmpEntry> > tmpData_;  // first index: jetFlavor
  std::vector<bool> useAbsEta_;                  // first index: jetFlavor
  std::map<std::string, std::shared_ptr<BTagCalibrationReaderImpl>> otherSysTypeReaders_;
};


BTagCalibrationReader::BTagCalibrationReaderImpl::BTagCalibrationReaderImpl(
                                             BTagEntry::OperatingPoint op,
                                             const std::string & sysType,
                                             const std::vector<std::string> & otherSysTypes):
  op_(op),
  sysType_(sysType),
  tmpData_(3),
  useAbsEta_(3, true)
{
  for (const std::string & ost : otherSysTypes) {
    if (otherSysTypeReaders_.count(ost)) {
std::cerr << "ERROR in BTagCalibration: "
            << "Every otherSysType should only be given once. Duplicate: "
            << ost;
throw std::exception();
    }
    //otherSysTypeReaders_[ost] = std::auto_ptr<BTagCalibrationReaderImpl>(
    otherSysTypeReaders_[ost] = std::unique_ptr<BTagCalibrationReaderImpl>(
        new BTagCalibrationReaderImpl(op, ost)
    );
  }
}

void BTagCalibrationReader::BTagCalibrationReaderImpl::load(
                                             const BTagCalibration & c,
                                             BTagEntry::JetFlavor jf,
                                             std::string measurementType)
{
  //if (tmpData_[jf].size()) {
  if (!tmpData_[jf].empty()) { 
std::cerr << "ERROR in BTagCalibration: "
          << "Data for this jet-flavor is already loaded: "
          << jf;
throw std::exception();
  }

  BTagEntry::Parameters params(op_, measurementType, sysType_);
  const std::vector<BTagEntry> &entries = c.getEntries(params);

  for (const auto &be : entries) {
    if (be.params.jetFlavor != jf) {
      continue;
    }

    TmpEntry te;
    te.etaMin = be.params.etaMin;
    te.etaMax = be.params.etaMax;
    te.ptMin = be.params.ptMin;
    te.ptMax = be.params.ptMax;
    te.discrMin = be.params.discrMin;
    te.discrMax = be.params.discrMax;

    if (op_ == BTagEntry::OP_RESHAPING) {
      te.func = TF1("", be.formula.c_str(),
                    be.params.discrMin, be.params.discrMax);
    } else {
      te.func = TF1("", be.formula.c_str(),
                    be.params.ptMin, be.params.ptMax);
    }

    tmpData_[be.params.jetFlavor].push_back(te);
    if (te.etaMin < 0) {
      useAbsEta_[be.params.jetFlavor] = false;
    }
  }

  for (auto & p : otherSysTypeReaders_) {
    p.second->load(c, jf, measurementType);
  }
}

double BTagCalibrationReader::BTagCalibrationReaderImpl::eval(
                                             BTagEntry::JetFlavor jf,
                                             float eta,
                                             float pt,
                                             float discr) const
{
  bool use_discr = (op_ == BTagEntry::OP_RESHAPING);
  if (useAbsEta_[jf] && eta < 0) {
    eta = -eta;
  }

  // search linearly through eta, pt and discr ranges and eval
  // future: find some clever data structure based on intervals
  const auto &entries = tmpData_.at(jf);
  for (unsigned i=0; i<entries.size(); ++i) {
    const auto &e = entries.at(i);
    if (
      e.etaMin <= eta && eta <= e.etaMax                   // find eta
      && e.ptMin < pt && pt <= e.ptMax                    // check pt
    ){
      if (use_discr) {                                    // discr. reshaping?
        if (e.discrMin <= discr && discr < e.discrMax) {  // check discr
          return e.func.Eval(discr);
        }
      } else {
        return e.func.Eval(pt);
      }
    }
  }

  return 0.;  // default value
}

double BTagCalibrationReader::BTagCalibrationReaderImpl::eval_auto_bounds(
                                             const std::string & sys,
                                             BTagEntry::JetFlavor jf,
                                             float eta,
                                             float pt,
                                             float discr) const
{
  auto sf_bounds_eta = min_max_eta(jf, discr);
  bool eta_is_out_of_bounds = false;

  //if (sf_bounds_eta.first < 0) sf_bounds_eta.first = -sf_bounds_eta.second;   
  
  //added from CMSSW master, not present in CMSSW_10_0
  if (sf_bounds_eta.first < 0) sf_bounds_eta.first = -sf_bounds_eta.second;  
  if (useAbsEta_[jf] && eta < 0) {
      eta = -eta;
  } 
    
  if (eta <= sf_bounds_eta.first || eta > sf_bounds_eta.second ) {
    eta_is_out_of_bounds = true;
  }
   
  if (eta_is_out_of_bounds) {
    return 1.;
  }


   auto sf_bounds = min_max_pt(jf, eta, discr);
   float pt_for_eval = pt;
   bool is_out_of_bounds = false;

   if (pt <= sf_bounds.first) {
    pt_for_eval = sf_bounds.first + .0001;
    is_out_of_bounds = true;
  } else if (pt > sf_bounds.second) {
    pt_for_eval = sf_bounds.second - .0001;
    is_out_of_bounds = true;
  }

  // get central SF (and maybe return)
  double sf = eval(jf, eta, pt_for_eval, discr);
  if (sys == sysType_) {
    return sf;
  }

  // get sys SF (and maybe return)
  if (!otherSysTypeReaders_.count(sys)) {
std::cerr << "ERROR in BTagCalibration: "
        << "sysType not available (maybe not loaded?): "
        << sys;
throw std::exception();
  }
  double sf_err = otherSysTypeReaders_.at(sys)->eval(jf, eta, pt_for_eval, discr);
  if (!is_out_of_bounds) {
    return sf_err;
  }

  // double uncertainty on out-of-bounds and return
  sf_err = sf + 2*(sf_err - sf);
  return sf_err;
}

std::pair<float, float> BTagCalibrationReader::BTagCalibrationReaderImpl::min_max_pt(
                                               BTagEntry::JetFlavor jf,
                                               float eta,
                                               float discr) const
{
  bool use_discr = (op_ == BTagEntry::OP_RESHAPING);
  if (useAbsEta_[jf] && eta < 0) {
    eta = -eta;
  }

  const auto &entries = tmpData_.at(jf);
  float min_pt = -1., max_pt = -1.;
  for (const auto & e: entries) {
    if (
      e.etaMin <= eta && eta <=e.etaMax                   // find eta
    ){
      if (min_pt < 0.) {                                  // init
        min_pt = e.ptMin;
        max_pt = e.ptMax;
        continue;
      }

      if (use_discr) {                                    // discr. reshaping?
        if (e.discrMin <= discr && discr < e.discrMax) {  // check discr
          min_pt = min_pt < e.ptMin ? min_pt : e.ptMin;
          max_pt = max_pt > e.ptMax ? max_pt : e.ptMax;
        }
      } else {
        min_pt = min_pt < e.ptMin ? min_pt : e.ptMin;
        max_pt = max_pt > e.ptMax ? max_pt : e.ptMax;
      }
    }
  }

  return std::make_pair(min_pt, max_pt);
}

std::pair<float, float> BTagCalibrationReader::BTagCalibrationReaderImpl::min_max_eta(
                                               BTagEntry::JetFlavor jf,
                                               float discr) const
{
  bool use_discr = (op_ == BTagEntry::OP_RESHAPING);

  const auto &entries = tmpData_.at(jf);
  float min_eta = 0., max_eta = 0.;
  for (const auto & e: entries) {

      if (use_discr) {                                    // discr. reshaping?
        if (e.discrMin <= discr && discr < e.discrMax) {  // check discr
          min_eta = min_eta < e.etaMin ? min_eta : e.etaMin;
          max_eta = max_eta > e.etaMax ? max_eta : e.etaMax;
        }
      } else {
        min_eta = min_eta < e.etaMin ? min_eta : e.etaMin;
        max_eta = max_eta > e.etaMax ? max_eta : e.etaMax;
      }
    }


  return std::make_pair(min_eta, max_eta);
}


BTagCalibrationReader::BTagCalibrationReader(BTagEntry::OperatingPoint op,
                                             const std::string & sysType,
                                             const std::vector<std::string> & otherSysTypes):
  pimpl(new BTagCalibrationReaderImpl(op, sysType, otherSysTypes)) {}

void BTagCalibrationReader::load(const BTagCalibration & c,
                                 BTagEntry::JetFlavor jf,
                                 const std::string & measurementType)
{
  pimpl->load(c, jf, measurementType);
}

double BTagCalibrationReader::eval(BTagEntry::JetFlavor jf,
                                   float eta,
                                   float pt,
                                   float discr) const
{
  return pimpl->eval(jf, eta, pt, discr);
}

double BTagCalibrationReader::eval_auto_bounds(const std::string & sys,
                                               BTagEntry::JetFlavor jf,
                                               float eta,
                                               float pt,
                                               float discr) const
{
  return pimpl->eval_auto_bounds(sys, jf, eta, pt, discr);
}

std::pair<float, float> BTagCalibrationReader::min_max_pt(BTagEntry::JetFlavor jf,
                                                          float eta,
                                                          float discr) const
{
  return pimpl->min_max_pt(jf, eta, discr);
}


//# load libraries for accessing b-tag scale factors (SFs) from conditions database
//for library in ["libCondFormatsBTauObjects", "libCondToolsBTau"]:
//    if library not in ROOT.gSystem.GetLibraries():
//        print("Load Library '%s'" % library.replace("lib", ""))
//        ROOT.gSystem.Load(library)

// initialize BTagCalibrationReader
// (cf. https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagCalibration )

//BTagCalibration calibration_Legacy2016("deepjet", "DeepJet_2016LegacySF_V1.csv");
//BTagCalibration calibration_2017("deepjet", "DeepFlavour_94XSF_V3_B_F.csv");
//BTagCalibration calibration_2018("deepjet", "DeepJet_102XSF_V1.csv");
BTagCalibration calibration_UL2016APV("deepjet", "DeepJet_106XUL16preVFPSF_v1_new.csv");
BTagCalibration calibration_UL2016("deepjet", "DeepJet_106XUL16postVFPSF_v2_new.csv");
BTagCalibration calibration_UL2017("deepjet", "DeepJet_106XUL17_v3_new.csv");
BTagCalibration calibration_UL2018("deepjet", "DeepJet_106XUL18_v2_new.csv");

vector<string> v_systs{"up", "down"};
//BTagCalibrationReader reader_0_Legacy2016(BTagEntry::OP_LOOSE,"central",v_systs); //0 is wp_btv
//BTagCalibrationReader reader_1_Legacy2016(BTagEntry::OP_MEDIUM,"central",v_systs);
//BTagCalibrationReader reader_2_Legacy2016(BTagEntry::OP_TIGHT,"central",v_systs); 
//BTagCalibrationReader reader_0_2017(BTagEntry::OP_LOOSE,"central",v_systs); //0 is wp_btv
//BTagCalibrationReader reader_1_2017(BTagEntry::OP_MEDIUM,"central",v_systs);
//BTagCalibrationReader reader_2_2017(BTagEntry::OP_TIGHT,"central",v_systs); 
//BTagCalibrationReader reader_0_2018(BTagEntry::OP_LOOSE,"central",v_systs); //0 is wp_btv
//BTagCalibrationReader reader_1_2018(BTagEntry::OP_MEDIUM,"central",v_systs);
BTagCalibrationReader reader_2_2018(BTagEntry::OP_TIGHT,"central",v_systs); 
BTagCalibrationReader reader_0_UL2016APV(BTagEntry::OP_LOOSE,"central",v_systs); //0 is wp_btv
BTagCalibrationReader reader_1_UL2016APV(BTagEntry::OP_MEDIUM,"central",v_systs);
BTagCalibrationReader reader_2_UL2016APV(BTagEntry::OP_TIGHT,"central",v_systs); 
BTagCalibrationReader reader_0_UL2016(BTagEntry::OP_LOOSE,"central",v_systs); //0 is wp_btv
BTagCalibrationReader reader_1_UL2016(BTagEntry::OP_MEDIUM,"central",v_systs);
BTagCalibrationReader reader_2_UL2016(BTagEntry::OP_TIGHT,"central",v_systs); 
BTagCalibrationReader reader_0_UL2017(BTagEntry::OP_LOOSE,"central",v_systs); //0 is wp_btv
BTagCalibrationReader reader_1_UL2017(BTagEntry::OP_MEDIUM,"central",v_systs);
BTagCalibrationReader reader_2_UL2017(BTagEntry::OP_TIGHT,"central",v_systs); 
BTagCalibrationReader reader_0_UL2018(BTagEntry::OP_LOOSE,"central",v_systs); //0 is wp_btv
BTagCalibrationReader reader_1_UL2018(BTagEntry::OP_MEDIUM,"central",v_systs);
BTagCalibrationReader reader_2_UL2018(BTagEntry::OP_TIGHT,"central",v_systs); 


/*
    
reader_0_UL2016APV.load(calibration_UL2016APV, 0, "comb"); //0 is flavor_btv
reader_0_UL2016APV.load(calibration_UL2016APV, 1, "comb");
reader_0_UL2016APV.load(calibration_UL2016APV, 2, "incl");
reader_1_UL2016APV.load(calibration_UL2016APV, 0, "comb"); //0 is flavor_btv
reader_1_UL2016APV.load(calibration_UL2016APV, 1, "comb");
reader_1_UL2016APV.load(calibration_UL2016APV, 2, "incl");
reader_2_UL2016APV.load(calibration_UL2016APV, 0, "comb"); //0 is flavor_btv
reader_2_UL2016APV.load(calibration_UL2016APV, 1, "comb");
reader_2_UL2016APV.load(calibration_UL2016APV, 2, "incl");
    
reader_0_UL2016.load(calibration_UL2016, 0, "comb"); //0 is flavor_btv
reader_0_UL2016.load(calibration_UL2016, 1, "comb");
reader_0_UL2016.load(calibration_UL2016, 2, "incl");
reader_1_UL2016.load(calibration_UL2016, 0, "comb"); //0 is flavor_btv
reader_1_UL2016.load(calibration_UL2016, 1, "comb");
reader_1_UL2016.load(calibration_UL2016, 2, "incl");
reader_2_UL2016.load(calibration_UL2016, 0, "comb"); //0 is flavor_btv
reader_2_UL2016.load(calibration_UL2016, 1, "comb");
reader_2_UL2016.load(calibration_UL2016, 2, "incl");
    
reader_0_UL2017.load(calibration_UL2017, 0, "comb"); //0 is flavor_btv
reader_0_UL2017.load(calibration_UL2017, 1, "comb");
reader_0_UL2017.load(calibration_UL2017, 2, "incl");
reader_1_UL2017.load(calibration_UL2017, 0, "comb"); //0 is flavor_btv
reader_1_UL2017.load(calibration_UL2017, 1, "comb");
reader_1_UL2017.load(calibration_UL2017, 2, "incl");
reader_2_UL2017.load(calibration_UL2017, 0, "comb"); //0 is flavor_btv
reader_2_UL2017.load(calibration_UL2017, 1, "comb");
reader_2_UL2017.load(calibration_UL2017, 2, "incl");
    
reader_0_UL2018.load(calibration_UL2018, 0, "comb"); //0 is flavor_btv
reader_0_UL2018.load(calibration_UL2018, 1, "comb");
reader_0_UL2018.load(calibration_UL2018, 2, "incl");
reader_1_UL2018.load(calibration_UL2018, 0, "comb"); //0 is flavor_btv
reader_1_UL2018.load(calibration_UL2018, 1, "comb");
reader_1_UL2018.load(calibration_UL2018, 2, "incl");
reader_2_UL2018.load(calibration_UL2018, 0, "comb"); //0 is flavor_btv
reader_2_UL2018.load(calibration_UL2018, 1, "comb");
reader_2_UL2018.load(calibration_UL2018, 2, "incl");


int getFlavorBTV(int flavor){
    //Maps hadronFlavor to BTV flavor:
    //Note the flavor convention: hadronFlavor is b = 5, c = 4, f = 0
    //Convert them to the btagging group convention of 0, 1, 2
    int flavor_btv;
    if (abs(flavor) == 5) flavor_btv = 0;
    else if (abs(flavor) == 4) flavor_btv = 1;
    else if (abs(flavor) == 0 || abs(flavor) == 1 || abs(flavor) == 2 || abs(flavor) == 3 || abs(flavor) == 21) flavor_btv = 2;
    else{
        cout<<"WARNING: Unknown flavor "<<flavor<<"setting b-tagging SF to -1!"<<endl;
        return -1
        }
    return flavor_btv;
}

RVec<RVec<float>> btagSF(string sample, string WP, rvec_f Jet_pt, rvec_f Jet_eta, rvec_i Jet_hadronFlavour, rvec_f Jet_btagDeepFlavB, string era, string wp){
    RVec<RVec<float>> result;
    float max_abs_eta = 2.4;
    string discr;
    //if (algo == "csvv2") discr = "btagCSVV2";
    //else if (algo == "deepcsv") discr = "btagDeepB";
    //else if (algo == "cmva") discr = "btagCMVA";
    //else if (algo == "deepjet") discr = "btagDeepFlavB";
    //else cout<<"ERROR: Invalid algorithm "<< algo<<endl;
    
    BTagCalibrationReader reader;
    if (era == "Legacy2016"){
        if (wp == "L") reader = reader_0_Legacy2016;
        else if (wp == "M") reader = reader_1_Legacy2016;
        else reader = reader_2_Legacy2016;
    }
    else if (era == "2017"){
        if (wp == "L") reader = reader_0_2017;
        else if (wp == "M") reader = reader_1_2017;
        else reader = reader_2_2017;
    }
    else if (era == "2018"){
        if (wp == "L") reader = reader_0_2018;
        else if (wp == "M") reader = reader_1_2018;
        else reader = reader_2_2018;
    }
    else if (era == "UL2016APV"){
        if (wp == "L") reader = reader_0_UL2016APV;
        else if (wp == "M") reader = reader_1_UL2016APV;
        else reader = reader_2_UL2016APV;
    }
    else if (era == "UL2016"){
        if (wp == "L") reader = reader_0_UL2016;
        else if (wp == "M") reader = reader_1_UL2016;
        else reader = reader_2_UL2016;
    }
    else if (era == "UL2017"){
        if (wp == "L") reader = reader_0_UL2017;
        else if (wp == "M") reader = reader_1_UL2017;
        else reader = reader_2_UL2017;
    }
    else {
        if (wp == "L") reader = reader_0_UL2018;
        else if (wp == "M") reader = reader_1_UL2018;
        else reader = reader_2_UL2018;
    }
    
    for(int i = 0; i<Jet_pt.size(); i++){
        RVec<float> sfs;
        float pt = Jet_pt[i];
        float eta = Jet_eta[i];
        int flavor_btv_int = getFlavorBTV(Jet_hadronFlavour[i]);
        
        float discr = Jet_btagDeepFlavB[i];
        
        float epsilon = 1.e-3;
        
        if (eta <= -max_abs_eta) eta = -max_abs_eta + epsilon;
        if (eta >= +max_abs_eta) eta = +max_abs_eta - epsilon;

        if(flavor_btv == 0){
            sf = reader.eval_auto_bounds("central", BTagEntry::FLAV_B, eta, pt);
            sf_up = reader.eval_auto_bounds("up", BTagEntry::FLAV_B, eta, pt);
            sf_down = reader.eval_auto_bounds("down", BTagEntry::FLAV_B, eta, pt);
        }
        else if(flavor_btv == 1){
            sf = reader.eval_auto_bounds("central", BTagEntry::FLAV_C, eta, pt);
            sf_up = reader.eval_auto_bounds("up", BTagEntry::FLAV_C, eta, pt);
            sf_down = reader.eval_auto_bounds("down", BTagEntry::FLAV_C, eta, pt);
        }
        else {
            sf = reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, eta, pt);
            sf_up = reader.eval_auto_bounds("up", BTagEntry::FLAV_UDSG, eta, pt);
            sf_down = reader.eval_auto_bounds("down", BTagEntry::FLAV_UDSG, eta, pt);
        }
        
        float sf, sf_up, sf_down;
        
        // check if SF is OK
        if (sf < 0.01) sf = 1.;
        if (sf_up < 0.01) sf_up = 1.;
        if (sf_down < 0.01) sf_down = 1.;
        
        sfs[0] = sf;
        sfs[1] = sf_up;
        sfs[2] = sf_down;
        
        result.emplace_back(sfs);
    }
    
    return result;
}
*/

/*
if isMC:
        f.write("metCorrector = createJMECorrector(isMC="+str(isMC)+", dataYear=\""+str(sample.year)+"\", jesUncert='All', applyHEMfix=True)\n")
        f.write("fatJetCorrector = createJMECorrector(isMC="+str(isMC)+", dataYear=\""+str(sample.year)+"\", jesUncert='All', applyHEMfix=True, jetType = 'AK8PFPuppi')\n")
    else: 
        f.write("metCorrector = createJMECorrector(isMC="+str(isMC)+", dataYear=\""+str(sample.year)+"\", runPeriod='"+str(sample.runP)+"', applyHEMfix=True, jesUncert='All')\n")
        f.write("fatJetCorrector = createJMECorrector(isMC="+str(isMC)+", dataYear=\""+str(sample.year)+"\", runPeriod='"+str(sample.runP)+"', jesUncert='All', applyHEMfix=True, jetType = 'AK8PFPuppi')\n")
*/

RVec<float> ones(int n){
    RVec<float> result;
    for(int i; i<n; i++){
        result.emplace_back(1.);
    }
    return result;
}

#ifndef CMSJMECalculators_JMESystematicsCalculators_H
#define CMSJMECalculators_JMESystematicsCalculators_H

#include <map>
#include <ROOT/RVec.hxx>
#include "JetResolution.h"
#include "JetCorrectorParameters.h"
#include "SimpleJetCorrectionUncertainty.h"
#include "FactorizedJetCorrectorCalculator.h"
#include "FormulaEvaluator.h"

class TRandom3;

namespace rdfhelpers {

class ModifiedPtMCollection { // map of variation collections
public:
  using compv_t = ROOT::VecOps::RVec<float>;

  ModifiedPtMCollection() = default;
  ModifiedPtMCollection(std::size_t n, const compv_t& pt, const compv_t& mass)
    : m_pt(n, pt), m_mass(n, mass) {}

  std::size_t size() const { return m_pt.size(); }

  const compv_t& pt(std::size_t i) const { return m_pt[i]; }
  const compv_t& mass(std::size_t i) const { return m_mass[i]; }

  void set(std::size_t i, const compv_t& pt, const compv_t& mass) {
    m_pt[i] = pt;
    m_mass[i] = mass;
  }
  void set(std::size_t i, compv_t&& pt, compv_t&& mass) {
    m_pt[i] = std::move(pt);
    m_mass[i] = std::move(mass);
  }
private:
  std::vector<compv_t> m_pt;
  std::vector<compv_t> m_mass;
};

class ModifiedPtMMsdCollection { // for fat jets
public:
  using compv_t = ROOT::VecOps::RVec<float>;

  ModifiedPtMMsdCollection() = default;
  ModifiedPtMMsdCollection(std::size_t nPt, const compv_t& pt, const std::size_t nM, const compv_t& mass, const compv_t& msd)
    : m_pt(nPt, pt), m_mass(nM, mass), m_msd(nM, msd) {}

  std::size_t size() const { return m_pt.size(); }
  std::size_t sizeM() const { return m_mass.size(); }

  const compv_t& pt(std::size_t i) const { return m_pt[i]; }
  const compv_t& mass(std::size_t i) const { return m_mass[i]; }
  const compv_t& msoftdrop(std::size_t i) const { return m_msd[i]; }

  void set(std::size_t i, const compv_t& pt, const compv_t& mass, const compv_t& msd) {
    m_pt[i] = pt;
    m_mass[i] = mass;
    m_msd[i] = msd;
  }
  void setM(std::size_t i, const compv_t& mass, const compv_t& msd) {
    m_mass[i] = mass;
    m_msd[i] = msd;
  }
  void set(std::size_t i, compv_t&& pt, compv_t&& mass, compv_t&& msd) {
    m_pt[i] = std::move(pt);
    m_mass[i] = std::move(mass);
    m_msd[i] = std::move(msd);
  }
  void setM(std::size_t i, compv_t&& mass, compv_t&& msd) {
    m_mass[i] = std::move(mass);
    m_msd[i] = std::move(msd);
  }
private:
  std::vector<compv_t> m_pt;
  std::vector<compv_t> m_mass;
  std::vector<compv_t> m_msd;
};

class ModifiedMET {
public:
  using compv_t = ROOT::VecOps::RVec<double>;

  ModifiedMET() = default;
  // initialize with the nominal value for all variations
  ModifiedMET(std::size_t n, double px_nom, double py_nom)
    : m_px(n, px_nom), m_py(n, py_nom) {}

  std::size_t size() const { return m_px.size(); }
  const compv_t& px() const { return m_px; }
  const compv_t& py() const { return m_py; }
  double px (std::size_t i) const { return m_px[i]; }
  double py (std::size_t i) const { return m_py[i]; }
  double pt (std::size_t i) const { return std::sqrt(m_px[i]*m_px[i]+m_py[i]*m_py[i]); }
  double phi(std::size_t i) const { return std::atan2(m_py[i], m_px[i]); }

  void setXY(std::size_t i, double dpx, double dpy) {
    m_px[i] = dpx;
    m_py[i] = dpy;
  }
  void addR_proj(std::size_t i, double cosphi, double sinphi, double dp) {
    m_px[i] += dp*cosphi;
    m_py[i] += dp*sinphi;
  }
private:
  compv_t m_px;
  compv_t m_py;
};
}

class JetMETVariationsCalculatorBase {
public:
  using p4compv_t = ROOT::VecOps::RVec<float>;

  JetMETVariationsCalculatorBase() = default;

  // set up smearing (and JER systematics)
  void setSmearing(const std::string& ptResolution, const std::string& ptResolutionSF, bool splitJER, bool doGenMatch, float genMatch_maxDR=-1., float genMatch_maxDPT=-1.)
  {
    m_doSmearing = true;
    m_jetPtRes   = JME::JetResolution(ptResolution);
    m_jetEResSF = JME::JetResolutionScaleFactor(ptResolutionSF);
    m_splitJER = splitJER;
    m_smearDoGenMatch = doGenMatch;
    m_genMatch_dR2max = genMatch_maxDR*genMatch_maxDR;
    m_genMatch_dPtmax = genMatch_maxDPT;
  }

  void setJEC(const std::vector<JetCorrectorParameters>& jecParams);
  void setAddHEM2018Issue(bool enable) { m_addHEM2018Issue = enable; }

  void addJESUncertainty(const std::string& name, const JetCorrectorParameters& params)
  {
    m_jesUncSources.emplace(std::piecewise_construct,
        std::forward_as_tuple(name),
        std::forward_as_tuple(params));
  }
protected:
  std::size_t findGenMatch(const double pt, const float eta, const float phi, const ROOT::VecOps::RVec<float>& gen_pt, const ROOT::VecOps::RVec<float>& gen_eta, const ROOT::VecOps::RVec<float>& gen_phi, const double resolution ) const;

  // config options
  bool m_doSmearing{false}, m_smearDoGenMatch;      // default: yes, yes
  bool m_addHEM2018Issue{false}, m_splitJER{false}; // default: no, no
  float m_genMatch_dR2max, m_genMatch_dPtmax;       // default: R/2 (0.2) and 3
  // parameters and helpers
  JME::JetResolution m_jetPtRes;
  JME::JetResolutionScaleFactor m_jetEResSF;
  std::unique_ptr<FactorizedJetCorrectorCalculator> m_jetCorrector;
  std::unordered_map<std::string,SimpleJetCorrectionUncertainty> m_jesUncSources;
};

class JetVariationsCalculator : public JetMETVariationsCalculatorBase {
public:
  using result_t = rdfhelpers::ModifiedPtMCollection;

  JetVariationsCalculator() = default;

  static JetVariationsCalculator create(
      const std::vector<std::string>& jecParams,
      const std::vector<std::pair<std::string,std::string>>& jesUncertainties,
      bool addHEM2018Issue,
      const std::string& ptResolution, const std::string& ptResolutionSF, bool splitJER,
      bool doGenMatch, float genMatch_maxDR, float genMatch_maxDPT);

  std::vector<std::string> available(const std::string& attr = {}) const;
  // interface for NanoAOD
  result_t produce(
      const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
      const p4compv_t& jet_rawcorr, const p4compv_t& jet_area, const ROOT::VecOps::RVec<int>& jet_jetId, const float rho,
      // MC-only
      const std::uint32_t seed,
      const p4compv_t& genjet_pt, const p4compv_t& genjet_eta, const p4compv_t& genjet_phi, const p4compv_t& genjet_mass
      ) const;
};

class Type1METVariationsCalculator : public JetMETVariationsCalculatorBase {
public:
  using result_t = rdfhelpers::ModifiedMET;

  Type1METVariationsCalculator() = default;

  static Type1METVariationsCalculator create(
      const std::vector<std::string>& jecParams, const std::vector<std::string>& jecParamsL1,
      float unclEnThreshold,
      const std::vector<std::pair<std::string,std::string>>& jesUncertainties,
      bool addHEM2018Issue,
      bool isT1SmearedMET,
      const std::string& ptResolution, const std::string& ptResolutionSF, bool splitJER,
      bool doGenMatch, float genMatch_maxDR, float genMatch_maxDPT);

  // additional settings: L1-only JEC
  void setL1JEC(const std::vector<JetCorrectorParameters>& jecParams);
  void setUnclusteredEnergyTreshold(float threshold) { m_unclEnThreshold = threshold; }
  void setIsT1SmearedMET(bool isT1SmearedMET) { m_isT1SmearedMET = isT1SmearedMET; }

  std::vector<std::string> available(const std::string& attr = {}) const;
  // interface for NanoAOD
  result_t produce(
      const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
      const p4compv_t& jet_rawcorr, const p4compv_t& jet_area,
      const p4compv_t& jet_muonSubtrFactor, const p4compv_t& jet_neEmEF, const p4compv_t& jet_chEmEF, const ROOT::VecOps::RVec<int>& jet_jetId,
      const float rho,
      // MC-only
      const std::uint32_t seed,
      const p4compv_t& genjet_pt, const p4compv_t& genjet_eta, const p4compv_t& genjet_phi, const p4compv_t& genjet_mass,
      // MET-specific
      const float rawmet_phi, const float rawmet_pt,
      const float met_unclustenupdx, const float met_unclustenupdy,
      const p4compv_t& lowptjet_rawpt, const p4compv_t& lowptjet_eta, const p4compv_t& lowptjet_phi, const p4compv_t& lowptjet_area,
      const p4compv_t& lowptjet_muonSubtrFactor, const p4compv_t& lowptjet_neEmEF, const p4compv_t& lowptjet_chEmEF
      ) const;
protected:
  float m_unclEnThreshold = 15.;
  bool m_isT1SmearedMET = false;
  std::unique_ptr<FactorizedJetCorrectorCalculator> m_jetCorrectorL1;
  void addVariations(Type1METVariationsCalculator::result_t& out,
      const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
      const p4compv_t& jet_rawcorr, const p4compv_t& jet_area, const p4compv_t& jet_muonSubtrFactor,
      const p4compv_t& jet_neEmEF, const p4compv_t& jet_chEmEF,
      const ROOT::VecOps::RVec<int>& jet_jetId, const ROOT::VecOps::RVec<bool>& jet_mask,
      const float rho,
      const p4compv_t& genjet_pt, const p4compv_t& genjet_eta, const p4compv_t& genjet_phi, const p4compv_t& genjet_mass,
      TRandom3& rg) const;
};

class FixEE2017Type1METVariationsCalculator : public Type1METVariationsCalculator {
public:
  FixEE2017Type1METVariationsCalculator() = default;

  static FixEE2017Type1METVariationsCalculator create(
      const std::vector<std::string>& jecParams, const std::vector<std::string>& jecParamsL1,
      float unclEnThreshold,
      const std::vector<std::string>& jecParamsProd, const std::vector<std::string>& jecParamsL1Prod,
      const std::vector<std::pair<std::string,std::string>>& jesUncertainties,
      bool isT1SmearedMET,
      const std::string& ptResolution, const std::string& ptResolutionSF, bool splitJER,
      bool doGenMatch, float genMatch_maxDR, float genMatch_maxDPT);

  // additional settings: full and L1-only JEC used in production
  void setJECProd(const std::vector<JetCorrectorParameters>& jecParams);
  void setL1JECProd(const std::vector<JetCorrectorParameters>& jecParams);

  // interface for NanoAOD
  result_t produce(
      const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
      const p4compv_t& jet_rawcorr, const p4compv_t& jet_area,
      const p4compv_t& jet_muonSubtrFactor, const p4compv_t& jet_neEmEF, const p4compv_t& jet_chEmEF,
      const float rho,
      // MC-only
      const std::uint32_t seed,
      const p4compv_t& genjet_pt, const p4compv_t& genjet_eta, const p4compv_t& genjet_phi, const p4compv_t& genjet_mass,
      // MET-specific
      const float rawmet_phi, const float rawmet_pt, // "RawMET"
      const float met_unclustenupdx, const float met_unclustenupdy,
      const p4compv_t& lowptjet_rawpt, const p4compv_t& lowptjet_eta, const p4compv_t& lowptjet_phi, const p4compv_t& lowptjet_area,
      const p4compv_t& lowptjet_muonSubtrFactor, const p4compv_t& lowptjet_neEmEF, const p4compv_t& lowptjet_chEmEF,
      const float defmet_phi, const float defmet_pt, // "MET"
      const float t1met_phi, const float t1met_pt    // "METFixEE2017"
      ) const;
protected:
  std::unique_ptr<FactorizedJetCorrectorCalculator> m_jetCorrectorProd;
  std::unique_ptr<FactorizedJetCorrectorCalculator> m_jetCorrectorL1Prod;
  std::array<double,4> calculateFixEE2017Offset(ROOT::VecOps::RVec<bool>& jet_mask,
      const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
      const p4compv_t& jet_rawcorr, const p4compv_t& jet_area, const p4compv_t& jet_muonSubtrFactor,
      const float rho) const;
};

class FatJetVariationsCalculator : public JetMETVariationsCalculatorBase {
public:
  using result_t = rdfhelpers::ModifiedPtMMsdCollection;

  FatJetVariationsCalculator() = default;

  static FatJetVariationsCalculator create(
      const std::vector<std::string>& jecParams,
      const std::vector<std::pair<std::string,std::string>>& jesUncertainties,
      bool addHEM2018Issue,
      const std::string& ptResolution, const std::string& ptResolutionSF, bool splitJER,
      bool doGenMatch, float genMatch_maxDR, float genMatch_maxDPT,
      std::vector<double> jmsValues,
      std::vector<double> gmsValues,
      std::vector<double> jmrValues,
      std::vector<double> gmrValues,
      const std::string& puppiGenFormula,
      const std::array<double, 6>& puppi_reco_cen_params,
      const std::array<double, 6>& puppi_reco_for_params,
      const std::array<double, 6>& puppi_resol_cen_params,
      const std::array<double, 6>& puppi_resol_for_params);

  void setJMRValues(double nominal, double up=1., double down=1.) { m_jmrVals = {{ nominal, up, down }}; }
  void setGMRValues(double nominal, double up=1., double down=1.) { m_gmrVals = {{ nominal, up, down }}; }
  void setJMSValues(double nominal, double up=1., double down=1.) { m_jmsVals = {{ nominal, up/nominal, down/nominal }}; }
  void setGMSValues(double nominal, double up=1., double down=1.) { m_gmsVals = {{ nominal, up/nominal, down/nominal }}; }
  void setPuppiCorrections(const std::string& genFormula, const std::array<double, 6>& reco_cen_params, const std::array<double, 6>& reco_for_params, const std::array<double, 6>& resol_cen_params, const std::array<double, 6>& resol_for_params);

  std::vector<std::string> available(const std::string& attr = {}) const;
  // interface for NanoAOD
  result_t produce(
      const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
      const p4compv_t& jet_rawcorr, const p4compv_t& jet_area,
      const p4compv_t& jet_msoftdrop, const ROOT::VecOps::RVec<int>& jet_subJetIdx1, const ROOT::VecOps::RVec<int>& jet_subJetIdx2,
      const p4compv_t& subjet_pt, const p4compv_t& subjet_eta, const p4compv_t& subjet_phi, const p4compv_t& subjet_mass,
      const ROOT::VecOps::RVec<int>& jet_jetId, const float rho,
      // MC-only
      const std::uint32_t seed,
      const p4compv_t& genjet_pt, const p4compv_t& genjet_eta, const p4compv_t& genjet_phi, const p4compv_t& genjet_mass,
      const p4compv_t& gensubjet_pt, const p4compv_t& gensubjet_eta, const p4compv_t& gensubjet_phi, const p4compv_t& gensubjet_mass
      ) const;
private:
  std::unique_ptr<reco::FormulaEvaluator> m_puppiCorrGen, m_puppiPoly5;
  std::array<double,6> m_puppiCorrRecoParam_cen, m_puppiCorrRecoParam_for, m_puppiResolParam_cen, m_puppiResolParam_for;
  std::array<double,3> m_jmrVals = {{1., 1., 1.}}; // nominal, up, down
  std::array<double,3> m_gmrVals = {{1., 1., 1.}}; // nominal, up, down
  std::array<double,3> m_jmsVals = {{1., 1., 1.}}; // nominal, up/nominal, down/nominal
  std::array<double,3> m_gmsVals = {{1., 1., 1.}}; // nominal, up/nominal, down/nominal
};

#endif // CMSJMECalculators_JMESystematicsCalculators_H


#include "JMESystematicsCalculators.h"
#include "FactorizedJetCorrectorCalculator.h"

#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include "Math/VectorUtil.h"

#include "TRandom3.h"

#include <cassert>

// #define BAMBOO_JME_DEBUG // uncomment to debug

#ifdef BAMBOO_JME_DEBUG
#define LogDebug std::cout
#else
#define LogDebug if (false) std::cout
#endif

void JetMETVariationsCalculatorBase::setJEC(const std::vector<JetCorrectorParameters>& jecParams)
{
  if ( ! jecParams.empty() ) {
    m_jetCorrector = std::unique_ptr<FactorizedJetCorrectorCalculator>{new FactorizedJetCorrectorCalculator(jecParams)};
  }
}

namespace {
  TRandom3& getTRandom3(uint32_t seed) {
    static thread_local TRandom3 rg{};
    rg.SetSeed(seed);
    return rg;
  }

  double jetESmearFactor( double pt, double eOrig, float genPt, float ptRes, float sfUncert, double rand )
  {
    double smear = 1.;
    if ( genPt > 0. ) {
      smear = 1. + (sfUncert-1.)*(pt - genPt)/pt;
    } else if ( sfUncert > 1. ) {
      smear = 1. + rand*std::sqrt(sfUncert*sfUncert-1.);
    }
    if ( smear*eOrig < 1.e-2 ) {
      smear = 1.e-2/eOrig;
    }
    return smear;
  }

  // because something goes wrong with linking ROOT::Math::VectorUtil::Phi_mpi_pi
  template<typename T>
  T phi_mpi_pi(T angle) {
    if ( angle <= M_PI && angle > -M_PI ) {
      return angle;
    }
    if ( angle > 0 ) {
      const int n = static_cast<int>(.5*(angle*M_1_PI+1.));
      angle -= 2*n*M_PI;
    } else {
      const int n = static_cast<int>(-.5*(angle*M_1_PI-1.));
      angle += 2*n*M_PI;
    }
    return angle;
  }

  std::vector<float> fillVector(const std::vector<std::string>& names, const JME::JetParameters::value_type& jetParams)
  {
    static const std::unordered_map<std::string,JME::Binning> jmeBinningFromString = {
        {"JetEta", JME::Binning::JetEta},
        {"JetPt" , JME::Binning::JetPt},
        // TODO JetPhi, JetEMF, LepPx, LepPy, LepPz
        {"JetE"  , JME::Binning::JetE}
      };
    std::vector<float> result;
    result.reserve(names.size());
    for ( const auto& nm : names ) {
      const auto it_key = jmeBinningFromString.find(nm);
      if ( std::end(jmeBinningFromString) == it_key ) {
        throw std::runtime_error{"Unknown binning variable: "+nm};
      } else {
        const auto it_par = jetParams.find(it_key->second);
        if ( std::end(jetParams) == it_par ) {
          throw std::runtime_error{"Binning variable "+nm+" not found"};
        } else {
          result.push_back(it_par->second);
        }
      }
    }
    return result;
  }

  float getUncertainty(const SimpleJetCorrectionUncertainty& uncert, const JME::JetParameters::value_type& jetParams, bool direction)
  {
    const auto vx = fillVector(uncert.parameters().definitions().binVar(), jetParams);
    const auto vy = fillVector(uncert.parameters().definitions().parVar(), jetParams);
    return uncert.uncertainty(vx, vy[0], direction);
  }

  float deltaHEM2018Issue(float pt_nom, int jetId, float phi, float eta ) {
    float delta = 1.;
    if ( pt_nom > 15. && ( jetId & 0x2 ) && phi > -1.57 && phi < -0.87 ) {
      if ( eta > -2.5 && eta < -1.3 ) {
        delta = 0.8;
      } else if ( eta <= -2.5 && eta > -3. ) {
        delta = 0.65;
      }
    }
    return delta;
  }

  int jerSplitID(float pt, float eta) {
    const auto aEta = std::abs(eta);
    if ( aEta < 1.93 )
      return 0;
    else if ( aEta < 2.5 )
      return 1;
    else if ( aEta < 3. )
      if ( pt < 50. )
        return 2;
      else
        return 3;
    else
      if ( pt < 50. )
        return 4;
      else
        return 5;
  }

  std::vector<JetCorrectorParameters> makeJCPList(const std::vector<std::string>& paths) {
    std::vector<JetCorrectorParameters> params;
    std::transform(std::begin(paths), std::end(paths), std::back_inserter(params),
      [] (const std::string& path) { return JetCorrectorParameters{path}; });
    return params;
  }

  template<typename CALC>
  void configureBaseCalc(CALC& calc,
    const std::vector<std::string>& jecParams,
    const std::vector<std::pair<std::string,std::string>>& jesUncertainties,
    const std::string& ptResolution, const std::string& ptResolutionSF, bool splitJER,
    bool doGenMatch, float genMatch_maxDR, float genMatch_maxDPT)
  {
    calc.setJEC(makeJCPList(jecParams));
    for (const auto& entry : jesUncertainties) {
      calc.addJESUncertainty(entry.first, JetCorrectorParameters{entry.second, entry.first});
    }
    if (!ptResolution.empty()) {
      calc.setSmearing(ptResolution, ptResolutionSF, splitJER,
                       doGenMatch, genMatch_maxDR, genMatch_maxDPT);
    }
  }

  template<typename CALC>
  void configureMETCalc_common(CALC& calc,
    const std::vector<std::string>& jecParams, const std::vector<std::string>& jecParamsL1,
    float unclEnThreshold,
    const std::vector<std::pair<std::string,std::string>>& jesUncertainties,
    bool isT1SmearedMET,
    const std::string& ptResolution, const std::string& ptResolutionSF, bool splitJER,
    bool doGenMatch, float genMatch_maxDR, float genMatch_maxDPT)
  {
    configureBaseCalc(calc, jecParams, jesUncertainties,
      ptResolution, ptResolutionSF, splitJER, doGenMatch, genMatch_maxDR, genMatch_maxDPT);
    calc.setL1JEC(makeJCPList(jecParamsL1));
    calc.setUnclusteredEnergyTreshold(unclEnThreshold);
    if (isT1SmearedMET) {
      calc.setIsT1SmearedMET(isT1SmearedMET);
    }
  }
}

// TODO with orig MET and jets (sumpx,sumpy): calc modif MET(sig), produce bigger results type

std::size_t JetMETVariationsCalculatorBase::findGenMatch(const double pt, const float eta, const float phi, const ROOT::VecOps::RVec<float>& gen_pt, const ROOT::VecOps::RVec<float>& gen_eta, const ROOT::VecOps::RVec<float>& gen_phi, const double resolution ) const
{
  auto dr2Min = std::numeric_limits<float>::max();
  std::size_t igBest{gen_pt.size()};
  LogDebug << "(DRs: ";
  for ( std::size_t ig{0}; ig != gen_pt.size(); ++ig ) {
    const auto dphi = phi_mpi_pi(gen_phi[ig]-phi);
    const auto deta = (gen_eta[ig]-eta);
    const auto dr2 = dphi*dphi + deta*deta;
    LogDebug << "dr2=" << dr2;
    if ( ( dr2 < dr2Min ) && ( dr2 < m_genMatch_dR2max ) ) {
      LogDebug << "->dpt=" << std::abs(gen_pt[ig]-pt) << ",res=" << resolution;
      if ( std::abs(gen_pt[ig]-pt) < m_genMatch_dPtmax*resolution ) {
        LogDebug << "->best:" << ig;
        dr2Min = dr2;
        igBest = ig;
      }
    }
    LogDebug << ", ";
  }
  LogDebug << ")";
  return igBest;
}

JetVariationsCalculator JetVariationsCalculator::create(
    const std::vector<std::string>& jecParams,
    const std::vector<std::pair<std::string,std::string>>& jesUncertainties,
    bool addHEM2018Issue,
    const std::string& ptResolution, const std::string& ptResolutionSF, bool splitJER,
    bool doGenMatch, float genMatch_maxDR, float genMatch_maxDPT)
{
  JetVariationsCalculator inst{};
  configureBaseCalc(inst, jecParams, jesUncertainties,
    ptResolution, ptResolutionSF, splitJER, doGenMatch, genMatch_maxDR, genMatch_maxDPT);
  inst.setAddHEM2018Issue(addHEM2018Issue);
  return std::move(inst);
}

JetVariationsCalculator::result_t JetVariationsCalculator::produce(
    const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
    const p4compv_t& jet_rawcorr, const p4compv_t& jet_area, const ROOT::VecOps::RVec<int>& jet_jetId, const float rho,
    const std::uint32_t seed,
    const p4compv_t& genjet_pt, const p4compv_t& genjet_eta, const p4compv_t& genjet_phi, const p4compv_t& genjet_mass ) const
{
  const auto nVariations = 1+( m_doSmearing ? 2*( m_splitJER ? 6 : 1 ) : 0 )+2*m_jesUncSources.size()+( m_addHEM2018Issue ? 2 : 0 ); // 1(nom)+2(JER)+2*len(JES)[+2(HEM)]
  LogDebug << "JME:: hello from JetVariations produce. Got " << jet_pt.size() << " jets" << std::endl;
  const auto nJets = jet_pt.size();
  result_t out{nVariations, jet_pt, jet_mass};
  ROOT::VecOps::RVec<double> pt_nom{jet_pt}, mass_nom{jet_mass};
  if ( m_jetCorrector ) {
    LogDebug << "JME:: reapplying JEC" << std::endl;
    FactorizedJetCorrectorCalculator::VariableValues vals;
    for ( std::size_t i{0}; i != nJets; ++i ) {
      vals.setJetEta(jet_eta[i]);
      vals.setJetPt(jet_pt[i]*(1.-jet_rawcorr[i]));
      vals.setJetA(jet_area[i]);
      vals.setRho(rho);
      const auto corr = m_jetCorrector->getCorrection(vals);
      if ( corr > 0. ) {
        const double newc = (1.-jet_rawcorr[i])*corr;
        pt_nom[i]   *= newc;
        mass_nom[i] *= newc;
      }
    }
#ifdef BAMBOO_JME_DEBUG
    LogDebug << "JME:: with reapplied JEC: ";
    for ( std::size_t i{0}; i != nJets; ++i ) {
      LogDebug << "(PT=" << pt_nom[i] << ", ETA=" << jet_eta[i] << ", PHI=" << jet_phi[i] << ", M=" << mass_nom[i] << ") ";
    }
    LogDebug << std::endl;
#endif
  } else {
    LogDebug << "JME:: Not reapplying JEC" << std::endl;
  }
  // smearing and JER
  std::size_t iVar = 1; // after nominal
  if ( m_doSmearing ) {
    LogDebug << "JME:: Smearing (seed=" << seed << ")" << std::endl;
    auto& rg = getTRandom3(seed);
    p4compv_t pt_jerUp(pt_nom.size(), 0.), mass_jerUp(mass_nom.size(), 0.);
    p4compv_t pt_jerDown(pt_nom.size(), 0.), mass_jerDown(mass_nom.size(), 0.);
    for ( std::size_t i{0}; i != nJets; ++i ) {
      const auto eOrig = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>(pt_nom[i], jet_eta[i], jet_phi[i], mass_nom[i]).E();
      double smearFactor_nom{1.}, smearFactor_down{1.}, smearFactor_up{1.};
      if ( pt_nom[i] > 0. ) {
        JME::JetParameters jPar{
            {JME::Binning::JetPt , pt_nom[i]},
            {JME::Binning::JetEta, jet_eta[i]},
            {JME::Binning::Rho   , rho} };
        const auto ptRes  = m_jetPtRes.getResolution(jPar);
        LogDebug << "JME:: JetParameters: pt=" << pt_nom[i] << ", eta=" << jet_eta[i] << ", rho=" << rho << "; ptRes=" << ptRes << std::endl;
        LogDebug << "JME:: ";
        float genPt = -1;
        if ( m_smearDoGenMatch ) {
          const auto iGen = findGenMatch(pt_nom[i], jet_eta[i], jet_phi[i], genjet_pt, genjet_eta, genjet_phi, ptRes*pt_nom[i]);
          if ( iGen != genjet_pt.size() ) {
            genPt = genjet_pt[iGen];
            LogDebug << "genPt=" << genPt << " ";
          }
        }
        const auto rand = ( genPt < 0. ) ? rg.Gaus(0, ptRes) : -1.;
        LogDebug << "jet_pt_resolution: " << ptRes << ", rand: " << rand << std::endl;
        smearFactor_nom  = jetESmearFactor(pt_nom[i], eOrig, genPt, ptRes, m_jetEResSF.getScaleFactor(jPar, Variation::NOMINAL), rand);
        smearFactor_down = jetESmearFactor(pt_nom[i], eOrig, genPt, ptRes, m_jetEResSF.getScaleFactor(jPar, Variation::DOWN   ), rand);
        smearFactor_up   = jetESmearFactor(pt_nom[i], eOrig, genPt, ptRes, m_jetEResSF.getScaleFactor(jPar, Variation::UP     ), rand);
        // LogDebug << "  scalefactors are NOMINAL=" << m_jetEResSF.getScaleFactor(jPar, Variation::NOMINAL) << ", DOWN=" << m_jetEResSF.getScaleFactor(jPar, Variation::DOWN) << ", UP=" << m_jetEResSF.getScaleFactor(jPar, Variation::UP) << std::endl;
        // LogDebug << "  smearfactors are NOMINAL=" << smearFactor_nom << ", DOWN=" << smearFactor_down << ", UP=" << smearFactor_up << std::endl;
      }
      pt_jerDown[i]   = pt_nom[i]*smearFactor_down;
      mass_jerDown[i] = mass_nom[i]*smearFactor_down;
      pt_jerUp[i]     = pt_nom[i]*smearFactor_up;
      mass_jerUp[i]   = mass_nom[i]*smearFactor_up;
      pt_nom[i]       *= smearFactor_nom;
      mass_nom[i]     *= smearFactor_nom;
    }
    if ( m_splitJER ) {
      ROOT::VecOps::RVec<int> jerBin(pt_nom.size(), -1);
      for ( std::size_t j{0}; j != nJets; ++j ) {
        jerBin[j] = jerSplitID(pt_nom[j], jet_eta[j]);
      }
      for ( int i{0}; i != 6; ++i ) {
        p4compv_t pt_jeriUp{pt_nom}, mass_jeriUp{mass_nom};
        p4compv_t pt_jeriDown{pt_nom}, mass_jeriDown{mass_nom};
        for ( std::size_t j{0}; j != nJets; ++j ) {
          if ( jerBin[j] == i ) {
            pt_jeriUp[j] = pt_jerUp[j];
            pt_jeriDown[j] = pt_jerDown[j];
            mass_jeriUp[j] = mass_jerUp[j];
            mass_jeriDown[j] = mass_jerDown[j];
          }
        }
        out.set(iVar++, std::move(pt_jeriUp)  , std::move(mass_jeriUp)  );
        out.set(iVar++, std::move(pt_jeriDown), std::move(mass_jeriDown));
      }
    } else {
      out.set(iVar++, std::move(pt_jerUp)  , std::move(mass_jerUp)  );
      out.set(iVar++, std::move(pt_jerDown), std::move(mass_jerDown));
    }
    LogDebug << "JME:: Done with smearing" << std::endl;
  } else {
    LogDebug << "JME:: No smearing" << std::endl;
  }
  out.set(0, pt_nom, mass_nom);

  // HEM issue 2018, see https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/2000.html
  if ( m_addHEM2018Issue ) {
    p4compv_t pt_down(pt_nom.size(), 0.), mass_down(mass_nom.size(), 0.);
    for ( std::size_t j{0}; j != nJets; ++j ) {
      const auto delta = deltaHEM2018Issue(pt_nom[j], jet_jetId[j], jet_phi[j], jet_eta[j]);
      pt_down[j] = pt_nom[j]*delta;
      mass_down[j] = mass_nom[j]*delta;
    }
    out.set(iVar++, pt_nom, mass_nom);
    out.set(iVar++, std::move(pt_down), std::move(mass_down));
  }
  // JES uncertainties
  for ( auto& jesUnc : m_jesUncSources ) {
    LogDebug << "JME:: evaluating JES uncertainty: " << jesUnc.first << std::endl;
    p4compv_t pt_jesDown(pt_nom.size(), 0.), mass_jesDown(mass_nom.size(), 0.);
    p4compv_t pt_jesUp(pt_nom.size(), 0.), mass_jesUp(mass_nom.size(), 0.);
    for ( std::size_t i{0}; i != nJets; ++i ) {
      const auto delta = getUncertainty(jesUnc.second, { {JME::Binning::JetPt, pt_nom[i]}, {JME::Binning::JetEta, jet_eta[i]} }, true);
      pt_jesDown[i]   = pt_nom[i]*(1.-delta);
      mass_jesDown[i] = mass_nom[i]*(1.-delta);
      pt_jesUp[i]     = pt_nom[i]*(1.+delta);
      mass_jesUp[i]   = mass_nom[i]*(1.+delta);
    }
    out.set(iVar++, std::move(pt_jesUp)  , std::move(mass_jesUp)  );
    out.set(iVar++, std::move(pt_jesDown), std::move(mass_jesDown));
  }

#ifdef BAMBOO_JME_DEBUG
  assert(iVar == out.size());
  LogDebug << "JME:: returning " << out.size() << " modified jet collections" << std::endl;
  const auto varNames = available();
  assert(varNames.size() == nVariations);
  for ( std::size_t i{0}; i != nVariations; ++i ) {
    LogDebug << "JME:: Jet_" << varNames[i] << ": ";
    for ( std::size_t j{0}; j != nJets; ++j ) {
      LogDebug << "(PT=" << out.pt(i)[j] << ", ETA=" << jet_eta[j] << ", PHI=" << jet_phi[j] << ", M=" << out.mass(i)[j] << ") ";
    }
    LogDebug << std::endl;
  }
#endif
  return out;
}

std::vector<std::string> JetVariationsCalculator::available(const std::string&) const
{
  std::vector<std::string> products = { "nominal" };
  if ( m_doSmearing ) {
    if ( m_splitJER ) {
      for ( int i = 0; i != 6; ++i ) {
        products.emplace_back("jer"+std::to_string(i)+"up");
        products.emplace_back("jer"+std::to_string(i)+"down");
      }
    } else {
      products.emplace_back("jerup");
      products.emplace_back("jerdown");
    }
  }
  if ( m_addHEM2018Issue ) {
    products.emplace_back("jesHEMIssueup");
    products.emplace_back("jesHEMIssuedown");
  }
  for ( const auto& src : m_jesUncSources ) {
    products.emplace_back("jes"+src.first+"up");
    products.emplace_back("jes"+src.first+"down");
  }
  return products;
}

FatJetVariationsCalculator FatJetVariationsCalculator::create(
    const std::vector<std::string>& jecParams,
    const std::vector<std::pair<std::string,std::string>>& jesUncertainties,
    bool addHEM2018Issue,
    const std::string& ptResolution, const std::string& ptResolutionSF, bool splitJER,
    bool doGenMatch, float genMatch_maxDR, float genMatch_maxDPT,
    std::vector<double> jmsValues, std::vector<double> gmsValues,
    std::vector<double> jmrValues, std::vector<double> gmrValues,
    const std::string& puppiGenFormula,
    const std::array<double, 6>& puppi_reco_cen_params, const std::array<double, 6>& puppi_reco_for_params,
    const std::array<double, 6>& puppi_resol_cen_params, const std::array<double, 6>& puppi_resol_for_params)
{
  FatJetVariationsCalculator inst{};
  configureBaseCalc(inst, jecParams, jesUncertainties,
    ptResolution, ptResolutionSF, splitJER, doGenMatch, genMatch_maxDR, genMatch_maxDPT);
  inst.setAddHEM2018Issue(addHEM2018Issue);
  if (jmsValues.size() == 3) {
    inst.setJMSValues(jmsValues[0], jmsValues[1], jmsValues[2]);
  } else if (jmsValues.size() == 1) {
    inst.setJMSValues(jmsValues[0]);
  }
  if (gmsValues.size() == 3) {
    inst.setGMSValues(gmsValues[0], gmsValues[1], gmsValues[2]);
  } else if (gmsValues.size() == 1) {
    inst.setGMSValues(gmsValues[0]);
  }
  if (jmrValues.size() == 3) {
    inst.setJMRValues(jmrValues[0], jmrValues[1], jmrValues[2]);
  } else if (jmrValues.size() == 1) {
    inst.setJMRValues(jmrValues[0]);
  }
  if (gmrValues.size() == 3) {
    inst.setGMRValues(gmrValues[0], gmrValues[1], gmrValues[2]);
  } else if (gmrValues.size() == 1) {
    inst.setGMRValues(gmrValues[0]);
  }
  if (!puppiGenFormula.empty()) {
    inst.setPuppiCorrections(puppiGenFormula,
        puppi_reco_cen_params, puppi_reco_for_params,
        puppi_resol_cen_params, puppi_resol_for_params);
  }
  return inst;
}

void FatJetVariationsCalculator::setPuppiCorrections(const std::string& genFormula, const std::array<double, 6>& reco_cen_params, const std::array<double, 6>& reco_for_params, const std::array<double, 6>& resol_cen_params, const std::array<double, 6>& resol_for_params) {
  m_puppiCorrGen = std::make_unique<reco::FormulaEvaluator>(genFormula);
  m_puppiPoly5 = std::make_unique<reco::FormulaEvaluator>("[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)");
  m_puppiCorrRecoParam_cen = reco_cen_params;
  m_puppiCorrRecoParam_for = reco_for_params;
  m_puppiResolParam_cen = resol_cen_params;
  m_puppiResolParam_for = resol_for_params;
}

FatJetVariationsCalculator::result_t FatJetVariationsCalculator::produce(
    const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
    const p4compv_t& jet_rawcorr, const p4compv_t& jet_area, const p4compv_t& jet_msoftdrop, const ROOT::VecOps::RVec<int>& jet_subJetIdx1, const ROOT::VecOps::RVec<int>& jet_subJetIdx2,
    const p4compv_t& subjet_pt, const p4compv_t& subjet_eta, const p4compv_t& subjet_phi, const p4compv_t& subjet_mass,
    const ROOT::VecOps::RVec<int>& jet_jetId, const float rho,
    const std::uint32_t seed,
    const p4compv_t& genjet_pt, const p4compv_t& genjet_eta, const p4compv_t& genjet_phi, const p4compv_t& genjet_mass,
    const p4compv_t& gensubjet_pt, const p4compv_t& gensubjet_eta, const p4compv_t& gensubjet_phi, const p4compv_t& gensubjet_mass) const
{
  const auto nVariations = 1+( m_doSmearing ? 2*( m_splitJER ? 6 : 1 ) : 0 )+2*m_jesUncSources.size()+( m_addHEM2018Issue ? 2 : 0 ); // 1(nom)+2(JER)+2*len(JES)[+2(HEM)]
  const bool doJMSVars = m_jmsVals[1] != 1.;
  const auto nVariationsM = nVariations + ( m_doSmearing ? 2 : 0 ) + ( doJMSVars ? 2 : 0 ); // 2(JMR)+2(JMS), both optional
  const auto nJets = jet_pt.size();
  LogDebug << "JME:: hello from FatJetVariations produce. Got " << nJets << " jets" << std::endl;
  LogDebug << "JME:: variations for PT: " << nVariations << " and for M, Msd: " << nVariationsM << std::endl;
  result_t out{nVariations, jet_pt, nVariationsM, jet_mass, jet_msoftdrop};

  ROOT::VecOps::RVec<double> pt_nom{jet_pt}, mass_raw{jet_mass};
  if ( m_jetCorrector ) {
    LogDebug << "JME:: reapplying JEC" << std::endl;
    FactorizedJetCorrectorCalculator::VariableValues vals;
    for ( std::size_t i{0}; i != nJets; ++i ) {
      vals.setJetEta(jet_eta[i]);
      vals.setJetPt(jet_pt[i]*(1.-jet_rawcorr[i]));
      vals.setJetA(jet_area[i]);
      vals.setRho(rho);
      const auto corr = m_jetCorrector->getCorrection(vals);
      if ( corr > 0. ) {
        const double newc = (1.-jet_rawcorr[i])*corr;
        pt_nom[i]   *= newc;
        mass_raw[i] *= newc;
      }
    }
#ifdef BAMBOO_JME_DEBUG
    LogDebug << "JME:: with reapplied JEC: ";
    for ( std::size_t i{0}; i != nJets; ++i ) {
      LogDebug << "(PT=" << pt_nom[i] << ", ETA=" << jet_eta[i] << ", PHI=" << jet_phi[i] << ", M=" << mass_raw[i] << ") ";
    }
    LogDebug << std::endl;
#endif
  } else {
    LogDebug << "JME:: Not reapplying JEC" << std::endl;
  }
  // calculate groomed P4 (and mass)
  using LVectorM = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>;
  auto jet_groomedP4 = std::vector<LVectorM>{nJets};
  ROOT::VecOps::RVec<double> msd_nom(nJets, 0.); // ?
  for ( std::size_t j{0}; j != nJets; ++j ) {
    if ( jet_subJetIdx1[j] >= 0 && jet_subJetIdx2[j] >= 0 ) {
      jet_groomedP4[j] = (
          LVectorM(subjet_pt[jet_subJetIdx1[j]], subjet_eta[jet_subJetIdx1[j]], subjet_phi[jet_subJetIdx1[j]], subjet_mass[jet_subJetIdx1[j]])
        + LVectorM(subjet_pt[jet_subJetIdx2[j]], subjet_eta[jet_subJetIdx2[j]], subjet_phi[jet_subJetIdx2[j]], subjet_mass[jet_subJetIdx2[j]]));
      // PUPPI SD mass correction https://github.com/cms-jet/PuppiSoftdropMassCorr/
      if ( m_puppiCorrGen->numberOfVariables() ) {
        const auto puppisd_corr = (
              m_puppiCorrGen->evaluate(std::array<double,1>{{ pt_nom[j] }}, std::array<double,0>{{}})
            * m_puppiPoly5->evaluate(std::array<double,1>{{ pt_nom[j] }},
               ( std::abs(jet_eta[j]) <= 1.3 ? m_puppiCorrRecoParam_cen : m_puppiCorrRecoParam_for ) ));
        LogDebug << "JME:: PUPPI gen mass correction: " << puppisd_corr << std::endl;
        jet_groomedP4[j].SetM(puppisd_corr*jet_groomedP4[j].M());
      }
      msd_nom[j] = jet_groomedP4[j].M();
      if ( msd_nom[j] < 0.0) {
        msd_nom[j] *= -1.;
      }
    }
  }
#ifdef BAMBOO_JME_DEBUG
  LogDebug << "JME:: Groomed momenta: ";
  for ( std::size_t i{0}; i != nJets; ++i ) {
    const auto& p4_g = jet_groomedP4[i];
    LogDebug << "(PT=" << p4_g.Pt() << ", ETA=" << p4_g.Eta() << ", PHI=" << p4_g.Phi() << ", M=" << msd_nom[i] << ") ";
  }
  LogDebug << std::endl;
#endif
  // apply nominal JMS (JMS variations by storing (up|down)/nom)
  auto mass_nom = mass_raw*m_jmsVals[0];
  msd_nom *= m_gmsVals[0];
  std::size_t iVar = 1; // after nominal
  // smearing and JER
  if ( m_doSmearing ) {
    p4compv_t pt_jerUp{pt_nom}, pt_jerDown{pt_nom};
    p4compv_t mass_jerUp{mass_nom}, mass_jerDown{mass_nom}, mass_jmrUp{mass_nom}, mass_jmrDown{mass_nom};
    p4compv_t msd_jerUp{msd_nom}, msd_jerDown{msd_nom}, msd_jmrUp{msd_nom}, msd_jmrDown{msd_nom};
    LogDebug << "JME:: Smearing (seed=" << seed << ")" << std::endl;
    auto& rg = getTRandom3(seed);
    for ( std::size_t i{0}; i != nJets; ++i ) {
      double jer_nom{1.}, jer_up{1.}, jer_down{1.}, jmr_nom{1.}, jmr_up{1.}, jmr_down{1.}, gmr_nom{1.}, gmr_up{1.}, gmr_down{1.};
      if ( pt_nom[i] > 0. || mass_raw[i] > 0. ) {
        const auto eOrig = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>(pt_nom[i], jet_phi[i], jet_eta[i], mass_raw[i]).E();
        JME::JetParameters jPar{
            {JME::Binning::JetPt , pt_nom[i]},
            {JME::Binning::JetEta, jet_eta[i]},
            {JME::Binning::Rho   , rho} };
        const auto ptRes  = m_jetPtRes.getResolution(jPar);
        LogDebug << "JME:: JetParameters: pt=" << pt_nom[i] << ", eta=" << jet_eta[i] << ", rho=" << rho << "; ptRes=" << ptRes << std::endl;
        LogDebug << "JME:: ";
        float genPt{-1}, genM{-1.};
        if ( m_smearDoGenMatch ) {
          const auto iGen = findGenMatch(pt_nom[i], jet_eta[i], jet_phi[i], genjet_pt, genjet_eta, genjet_phi, ptRes*pt_nom[i]);
          if ( iGen != genjet_pt.size() ) {
            genPt = genjet_pt[iGen];
            genM  = genjet_mass[iGen];
            LogDebug << "genPt=" << genPt << " genM=" << genM << " ";
          }
        }
        if ( pt_nom[i] > 0. ) {
          const auto rand = ( genPt < 0. ) ? rg.Gaus(0, ptRes) : -1.;
          LogDebug << "jet_pt_resolution: " << ptRes << ", rand: " << rand << std::endl;
          jer_nom  = jetESmearFactor(pt_nom[i], eOrig, genPt, ptRes, m_jetEResSF.getScaleFactor(jPar, Variation::NOMINAL), rand);
          jer_down = jetESmearFactor(pt_nom[i], eOrig, genPt, ptRes, m_jetEResSF.getScaleFactor(jPar, Variation::DOWN   ), rand);
          jer_up   = jetESmearFactor(pt_nom[i], eOrig, genPt, ptRes, m_jetEResSF.getScaleFactor(jPar, Variation::UP     ), rand);
        }
        // LogDebug << "  scalefactors are NOMINAL=" << m_jetEResSF.getScaleFactor(jPar, Variation::NOMINAL) << ", DOWN=" << m_jetEResSF.getScaleFactor(jPar, Variation::DOWN) << ", UP=" << m_jetEResSF.getScaleFactor(jPar, Variation::UP) << std::endl;
        // LogDebug << "  smearfactors are NOMINAL=" << jer_nom[i] << ", DOWN=" << jer_down[i] << ", UP=" << jer_up[i] << std::endl;
        // JMR
        if ( mass_raw[i] > 0. ) {
          if ( genM != -1. ) {
            LogDebug << "JME:: JMR with genM=" << genM << " and raw " << mass_raw[i];
            const auto dMoM = 1.-(genM/mass_raw[i]);
            jmr_nom  = 1.+(m_jmrVals[0]-1.)*dMoM;
            jmr_up   = 1.+(m_jmrVals[1]-1.)*dMoM;
            jmr_down = 1.+(m_jmrVals[2]-1.)*dMoM;
          } else {
            const auto mRes = m_puppiPoly5->evaluate(std::array<double,1>{{ pt_nom[i] }},
               ( std::abs(jet_eta[i]) <= 1.3 ? m_puppiResolParam_cen : m_puppiResolParam_for ) );
            LogDebug << "JME:: JMR parametric with resolution " << mRes;
            const auto rand = rg.Gaus(0, mRes);
            jmr_nom  = ( m_jmrVals[0] <= 1. ? 1. : rand*std::sqrt(m_jmrVals[0]*m_jmrVals[0]-1));
            jmr_up   = ( m_jmrVals[1] <= 1. ? 1. : rand*std::sqrt(m_jmrVals[1]*m_jmrVals[1]-1));
            jmr_down = ( m_jmrVals[2] <= 1. ? 1. : rand*std::sqrt(m_jmrVals[2]*m_jmrVals[2]-1));
          }
          if ( jmr_nom *mass_raw[i] < 1.e-2 ) { jmr_nom  = 1.e-2; }
          if ( jmr_up  *mass_raw[i] < 1.e-2 ) { jmr_up   = 1.e-2; }
          if ( jmr_down*mass_raw[i] < 1.e-2 ) { jmr_down = 1.e-2; }
        }
        LogDebug << "  mass smearfactors are NOMINAL=" << jmr_nom << ", DOWN=" << jmr_down << ", UP=" << jmr_up << std::endl;
      }
      if ( msd_nom[i] > 0. ) { // JMR for groomed
        // genGroomedJet
        int igsj1{-1}, igsj2{-1};
        for ( int ig{0}; ig != int(gensubjet_eta.size()); ++ig ) {
          const auto dphi = phi_mpi_pi(gensubjet_phi[ig]-jet_groomedP4[i].Phi());
          const auto deta = (gensubjet_eta[ig]-jet_groomedP4[i].Eta());
          const auto dr2 = dphi*dphi + deta*deta;
          if ( dr2 < 0.64 ) { // dR < 0.8
            if ( igsj1 == -1 ) {
              igsj1 = ig;
            } else {
              igsj2 = ig;
              break; // only need the first two
            }
          }
        }
        if ( igsj2 != -1 ) {
          LogDebug << "JME:: Matched gen-level subjets " << igsj1 << "," << igsj2 << std::endl;
          const auto genM = (
              LVectorM{gensubjet_pt[igsj1], gensubjet_eta[igsj1], gensubjet_phi[igsj1], gensubjet_mass[igsj1]}
            + LVectorM{gensubjet_pt[igsj2], gensubjet_eta[igsj2], gensubjet_phi[igsj2], gensubjet_mass[igsj2]}).M();
          const auto dMoM = 1.-(genM/jet_groomedP4[i].M()); // raw
          gmr_nom  = 1.+(m_gmrVals[0]-1.)*dMoM;
          gmr_up   = 1.+(m_gmrVals[1]-1.)*dMoM;
          gmr_down = 1.+(m_gmrVals[2]-1.)*dMoM;
        } else {
          const auto mRes = m_puppiPoly5->evaluate(std::array<double,1>{{ jet_groomedP4[i].Pt() }},
             ( std::abs(jet_eta[i]) <= 1.3 ? m_puppiResolParam_cen : m_puppiResolParam_for ) );
          const auto rand = rg.Gaus(0, mRes);
          gmr_nom  = ( m_gmrVals[0] <= 1. ? 1. : rand*std::sqrt(m_gmrVals[0]*m_gmrVals[0]-1));
          gmr_up   = ( m_gmrVals[1] <= 1. ? 1. : rand*std::sqrt(m_gmrVals[1]*m_gmrVals[1]-1));
          gmr_down = ( m_gmrVals[2] <= 1. ? 1. : rand*std::sqrt(m_gmrVals[2]*m_gmrVals[2]-1));
        }
        if ( gmr_nom *msd_nom[i] < 1.e-2 ) { gmr_nom  = 1.e-2; }
        if ( gmr_up  *msd_nom[i] < 1.e-2 ) { gmr_up   = 1.e-2; }
        if ( gmr_down*msd_nom[i] < 1.e-2 ) { gmr_down = 1.e-2; }
      }
      LogDebug << "  groomed mass smearfactors are NOMINAL=" << gmr_nom << ", DOWN=" << gmr_down << ", UP=" << gmr_up << std::endl;
      // fill variation arrays
      pt_jerDown  [i] *= jer_down;
      pt_jerUp    [i] *= jer_up;
      mass_jerUp  [i] *= jer_up  *jmr_nom;
      mass_jerDown[i] *= jer_down*jmr_nom;
      mass_jmrUp  [i] *= jer_nom *jmr_up  ;
      mass_jmrDown[i] *= jer_nom *jmr_down;
      msd_jerUp   [i] *= jer_up  *gmr_nom;
      msd_jerDown [i] *= jer_down*gmr_nom;
      msd_jmrUp   [i] *= jer_nom *gmr_up  ;
      msd_jmrDown [i] *= jer_nom *gmr_down;
      // finally apply JER and JMR to nominal
      pt_nom[i]   *= jer_nom;
      mass_nom[i] *= jer_nom*jmr_nom;
      msd_nom[i]  *= jer_nom*gmr_nom;
    }
    if ( m_splitJER ) {
      ROOT::VecOps::RVec<int> jerBin(pt_nom.size(), -1);
      for ( std::size_t j{0}; j != nJets; ++j ) {
        jerBin[j] = jerSplitID(pt_nom[j], jet_eta[j]);
      }
      for ( int i{0}; i != 6; ++i ) {
        p4compv_t pt_jeriUp{pt_nom}, mass_jeriUp{mass_nom}, msd_jeriUp{msd_nom};
        p4compv_t pt_jeriDown{pt_nom}, mass_jeriDown{mass_nom}, msd_jeriDown{msd_nom};
        for ( std::size_t j{0}; j != nJets; ++j ) {
          if ( jerBin[j] == i ) {
            pt_jeriUp[j] = pt_jerUp[j];
            pt_jeriDown[j] = pt_jerDown[j];
            mass_jeriUp[j] = mass_jerUp[j];
            mass_jeriDown[j] = mass_jerDown[j];
            msd_jeriUp[j] = msd_jerUp[j];
            msd_jeriDown[j] = msd_jerDown[j];
          }
        }
        out.set(iVar++, std::move(pt_jeriUp)  , std::move(mass_jeriUp)  , std::move(msd_jeriUp));
        out.set(iVar++, std::move(pt_jeriDown), std::move(mass_jeriDown), std::move(msd_jeriDown));
      }
    } else {
      out.set(iVar++, std::move(pt_jerUp)  , std::move(mass_jerUp)  , std::move(msd_jerUp));
      out.set(iVar++, std::move(pt_jerDown), std::move(mass_jerDown), std::move(msd_jerDown));
    }
    out.setM(nVariationsM-2, std::move(mass_jmrUp)  , std::move(msd_jmrUp));
    out.setM(nVariationsM-1, std::move(mass_jmrDown), std::move(msd_jmrDown));
    LogDebug << "JME:: Done with smearing" << std::endl;
  } else {
    LogDebug << "JME:: No smearing" << std::endl;
  }
  out.set(0, pt_nom, mass_nom, msd_nom);
  if ( doJMSVars ) { // mass_nom has nominal JMS, jmsVals/gmsVals are divided by that
    out.setM(nVariations  , mass_nom*m_jmsVals[1], msd_nom*m_gmsVals[1]); // UP
    out.setM(nVariations+1, mass_nom*m_jmsVals[2], msd_nom*m_gmsVals[2]); // DOWN
  }

  // HEM issue 2018, see https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/2000.html
  if ( m_addHEM2018Issue ) {
    p4compv_t pt_down(pt_nom.size(), 0.), mass_down(mass_nom.size(), 0.), msd_down(msd_nom.size(), 0.);
    for ( std::size_t j{0}; j != nJets; ++j ) {
      const auto delta = deltaHEM2018Issue(pt_nom[j], jet_jetId[j], jet_phi[j], jet_eta[j]);
      pt_down[j] = pt_nom[j]*delta;
      mass_down[j] = mass_nom[j]*delta;
      msd_down[j] = msd_nom[j]*delta;
    }
    out.set(iVar++, pt_nom, mass_nom, msd_nom);
    out.set(iVar++, std::move(pt_down), std::move(mass_down), std::move(msd_down));
  }
  // JES uncertainties
  for ( auto& jesUnc : m_jesUncSources ) {
    LogDebug << "JME:: evaluating JES uncertainty: " << jesUnc.first << std::endl;
    p4compv_t pt_jesDown(pt_nom.size(), 0.), mass_jesDown(mass_nom.size(), 0.), msd_jesDown(msd_nom.size(), 0.);
    p4compv_t pt_jesUp(pt_nom.size(), 0.), mass_jesUp(mass_nom.size(), 0.), msd_jesUp(msd_nom.size(), 0.);
    for ( std::size_t i{0}; i != nJets; ++i ) {
      const auto delta = getUncertainty(jesUnc.second, { {JME::Binning::JetPt, pt_nom[i]}, {JME::Binning::JetEta, jet_eta[i]} }, true);
      pt_jesDown  [i] = pt_nom  [i]*(1.-delta);
      mass_jesDown[i] = mass_nom[i]*(1.-delta);
      msd_jesDown [i] = msd_nom [i]*(1.-delta);
      pt_jesUp    [i] = pt_nom  [i]*(1.+delta);
      mass_jesUp  [i] = mass_nom[i]*(1.+delta);
      msd_jesUp   [i] = msd_nom [i]*(1.+delta);
    }
    out.set(iVar++, std::move(pt_jesUp)  , std::move(mass_jesUp)  , std::move(msd_jesUp)  );
    out.set(iVar++, std::move(pt_jesDown), std::move(mass_jesDown), std::move(msd_jesDown));
  }

#ifdef BAMBOO_JME_DEBUG
  assert(iVar == out.size());
  LogDebug << "JME:: returning " << out.size() << " modified jet collections" << std::endl;
  const auto varNames = available("mass");
  assert(varNames.size() == nVariationsM);
  for ( std::size_t i{0}; i != nVariations; ++i ) {
    LogDebug << "JME:: Jet_" << varNames[i] << ": ";
    for ( std::size_t j{0}; j != nJets; ++j ) {
      LogDebug << "(PT=" << out.pt(i)[j] << ", ETA=" << jet_eta[j] << ", PHI=" << jet_phi[j] << ", M=" << out.mass(i)[j] << ", Msd=" << out.msoftdrop(i)[j] << ") ";
    }
    LogDebug << std::endl;
  }
  for ( std::size_t i{nVariations}; i != nVariationsM; ++i ) {
    LogDebug << "JME:: Jet_" << varNames[i] << ": ";
    for ( std::size_t j{0}; j != nJets; ++j ) {
      LogDebug << "(PT=" << out.pt(0)[j] << ", ETA=" << jet_eta[j] << ", PHI=" << jet_phi[j] << ", M=" << out.mass(i)[j] << ", Msd=" << out.msoftdrop(i)[j] << ") ";
    }
    LogDebug << std::endl;
  }
#endif
  return out;
}

std::vector<std::string> FatJetVariationsCalculator::available(const std::string& attr) const
{
  std::vector<std::string> products = { "nominal" };
  if ( m_doSmearing ) {
    if ( m_splitJER ) {
      for ( int i = 0; i != 6; ++i ) {
        products.emplace_back("jer"+std::to_string(i)+"up");
        products.emplace_back("jer"+std::to_string(i)+"down");
      }
    } else {
      products.emplace_back("jerup");
      products.emplace_back("jerdown");
    }
  }
  if ( m_addHEM2018Issue ) {
    products.emplace_back("jesHEMIssueup");
    products.emplace_back("jesHEMIssuedown");
  }
  for ( const auto& src : m_jesUncSources ) {
    products.emplace_back("jes"+src.first+"up");
    products.emplace_back("jes"+src.first+"down");
  }
  if ( attr == "mass" || attr == "msoftdrop" ) {
    if ( m_jmsVals[1] != 1. ) {
      products.emplace_back("jmsup");
      products.emplace_back("jmsdown");
    }
    if ( m_doSmearing ) {
      products.emplace_back("jmrup");
      products.emplace_back("jmrdown");
    }
  }
  return products;
}

Type1METVariationsCalculator Type1METVariationsCalculator::create(
  const std::vector<std::string>& jecParams, const std::vector<std::string>& jecParamsL1,
  float unclEnThreshold,
  const std::vector<std::pair<std::string,std::string>>& jesUncertainties,
  bool addHEM2018Issue,
  bool isT1SmearedMET,
  const std::string& ptResolution, const std::string& ptResolutionSF, bool splitJER,
  bool doGenMatch, float genMatch_maxDR, float genMatch_maxDPT)
{
  auto inst = Type1METVariationsCalculator();
  configureMETCalc_common(
    inst, jecParams, jecParamsL1, unclEnThreshold, jesUncertainties,
    isT1SmearedMET, ptResolution, ptResolutionSF, splitJER,
    doGenMatch, genMatch_maxDR, genMatch_maxDPT
  );
  inst.setAddHEM2018Issue(addHEM2018Issue);
  return inst;
}

void Type1METVariationsCalculator::setL1JEC(const std::vector<JetCorrectorParameters>& jecParams)
{
  if ( ! jecParams.empty() ) {
    m_jetCorrectorL1 = std::unique_ptr<FactorizedJetCorrectorCalculator>{new FactorizedJetCorrectorCalculator(jecParams)};
  }
}

Type1METVariationsCalculator::result_t Type1METVariationsCalculator::produce(
    const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
    const p4compv_t& jet_rawcorr, const p4compv_t& jet_area,
    const p4compv_t& jet_muonSubtrFactor, const p4compv_t& jet_neEmEF, const p4compv_t& jet_chEmEF, const ROOT::VecOps::RVec<int>& jet_jetId,
    const float rho,
    const std::uint32_t seed,
    const p4compv_t& genjet_pt, const p4compv_t& genjet_eta, const p4compv_t& genjet_phi, const p4compv_t& genjet_mass,
    const float rawmet_phi, const float rawmet_pt,
    const float met_unclustenupdx, const float met_unclustenupdy,
    const p4compv_t& lowptjet_rawpt, const p4compv_t& lowptjet_eta, const p4compv_t& lowptjet_phi, const p4compv_t& lowptjet_area,
    const p4compv_t& lowptjet_muonSubtrFactor, const p4compv_t& lowptjet_neEmEF, const p4compv_t& lowptjet_chEmEF
    ) const
{
  const auto nVariations = 3+( m_doSmearing ? 2*( m_splitJER ? 6 : 1 ) : 0 )+2*m_jesUncSources.size()+( m_addHEM2018Issue ? 2 : 0 ); // 1(nom)+2(unclust)+2(JER[*6])+2*len(JES)[+2(HEM)]
  result_t out{nVariations, rawmet_pt*std::cos(rawmet_phi), rawmet_pt*std::sin(rawmet_phi)};
  auto& rg = getTRandom3(seed);
  LogDebug << "JME:: hello from Type1METVariations produce. Got " << jet_pt.size() << " jets and " << lowptjet_rawpt.size() << " low-PT jets" << std::endl;
  // normal jets
  addVariations(out, jet_pt, jet_eta, jet_phi, jet_mass,
      jet_rawcorr, jet_area, jet_muonSubtrFactor, jet_neEmEF, jet_chEmEF, jet_jetId,
      ROOT::VecOps::RVec<bool>(), rho, genjet_pt, genjet_eta, genjet_phi, genjet_mass, rg);
  // low-PT jets
  p4compv_t lowptjet_zero(lowptjet_rawpt.size(), 0.);
  addVariations(out, lowptjet_rawpt, lowptjet_eta, lowptjet_phi, lowptjet_zero,
      lowptjet_zero, lowptjet_area, lowptjet_muonSubtrFactor,
      ( lowptjet_neEmEF.empty() ? lowptjet_zero : lowptjet_neEmEF  ),
      ( lowptjet_chEmEF.empty() ? lowptjet_zero : lowptjet_chEmEF  ),
      ROOT::VecOps::RVec<int>(lowptjet_rawpt.size(), 0), ROOT::VecOps::RVec<bool>(), rho,
      genjet_pt, genjet_eta, genjet_phi, genjet_mass, rg);
  // unclustered energy, based on nominal (0)
  out.setXY(nVariations-2, out.px(0)+met_unclustenupdx, out.py(0)+met_unclustenupdy);
  out.setXY(nVariations-1, out.px(0)-met_unclustenupdx, out.py(0)-met_unclustenupdy);

#ifdef BAMBOO_JME_DEBUG
  LogDebug << "JME:: returning " << out.size() << " modified METs" << std::endl;
  const auto varNames = available();
  assert(varNames.size() == nVariations);
  for ( std::size_t i{0}; i != nVariations; ++i ) {
    LogDebug << "JME:: MET_" << varNames[i] << ": PT=" << out.pt(i) << ", PHI=" << out.phi(i) << std::endl;
  }
#endif
  return out;
}

// for a single jet collection
void Type1METVariationsCalculator::addVariations(Type1METVariationsCalculator::result_t& out,
    const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
    const p4compv_t& jet_rawcorr, const p4compv_t& jet_area, const p4compv_t& jet_muonSubtrFactor,
    const p4compv_t& jet_neEmEF, const p4compv_t& jet_chEmEF, const ROOT::VecOps::RVec<int>& jet_jetId, const ROOT::VecOps::RVec<bool>& jet_mask, const float rho,
    const p4compv_t& genjet_pt, const p4compv_t& genjet_eta, const p4compv_t& genjet_phi, const p4compv_t& genjet_mass,
    TRandom3& rg) const
{
  FactorizedJetCorrectorCalculator::VariableValues vals, valsL1;
  const auto nJets = jet_pt.size();
  for ( std::size_t i{0}; i != nJets; ++i ) {
    // L1 and full (L1L2L3) JEC for muon-subtracted jet
    double jet_pt_nom = jet_pt[i];
    double jet_mass_nom = jet_mass[i];
    const auto jet_pt_raw = jet_pt_nom*(1-jet_rawcorr[i]);
    vals.setJetEta(jet_eta[i]);
    vals.setJetPt(jet_pt_raw);
    vals.setJetA(jet_area[i]);
    vals.setRho(rho);
    LogDebug << "Jet #" << i << " ETA=" << jet_eta[i] << ", PT_raw=" << jet_pt_raw << ", area=" << jet_area[i] << std::endl;
    auto jecL1L2L3 = m_jetCorrector->getCorrection(vals);
    if ( jecL1L2L3 <= 0. ) {
      jecL1L2L3 = 1.;
    } else {
      jet_pt_nom = jet_pt_raw*jecL1L2L3;
      jet_mass_nom = jet_mass_nom*(1-jet_rawcorr[i])*jecL1L2L3;
    }
    valsL1.setJetEta(jet_eta[i]);
    valsL1.setJetPt(jet_pt_raw);
    valsL1.setJetA(jet_area[i]);
    valsL1.setRho(rho);
    auto jecL1 = m_jetCorrectorL1->getCorrection(valsL1);
    if ( jecL1     <= 0. ) { jecL1     = 1.; }
    const double jet_pt_raw_nomu = jet_pt_raw*(1-jet_muonSubtrFactor[i]);
    const double muon_pt = jet_pt_raw*jet_muonSubtrFactor[i];
    const auto jet_pt_nomuL1L2L3 = jet_pt_raw_nomu*jecL1L2L3;
    const auto jet_pt_nomuL1     = jet_pt_raw_nomu*jecL1;
    const auto jet_pt_L1L2L3 = jet_pt_nomuL1L2L3 + muon_pt;
    const auto jet_pt_L1     = jet_pt_nomuL1     + muon_pt;
    LogDebug << "jecL1L2L3=" << jecL1L2L3 << ", jecL1=" << jecL1 << "; PT_L1L2L3=" << jet_pt_L1L2L3 << ", PT_L1=" << jet_pt_L1 << ", PT_mu=" << muon_pt << std::endl;
    // JER / smearing
    double jerF_nom{1.}, jerF_up{1.}, jerF_down{1.};
    if ( m_doSmearing ) {
      const auto eOrig = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>(jet_pt_nom, jet_phi[i], jet_eta[i], jet_mass_nom).E();
      if ( jet_pt_nom > 0. ) {
        JME::JetParameters jPar{
            {JME::Binning::JetPt , jet_pt_nom},
            {JME::Binning::JetEta, jet_eta[i]},
            {JME::Binning::Rho   , rho} };
        const auto ptRes  = m_jetPtRes.getResolution(jPar);
        LogDebug << "JME:: JetParameters: pt=" << jet_pt_nom << ", eta=" << jet_eta[i] << ", rho=" << rho << "; ptRes=" << ptRes << std::endl;
        LogDebug << "JME:: ";
        float genPt = -1;
        if ( m_smearDoGenMatch ) {
          const auto iGen = findGenMatch(jet_pt_nom, jet_eta[i], jet_phi[i], genjet_pt, genjet_eta, genjet_phi, ptRes*jet_pt_nom);
          if ( iGen != genjet_pt.size() ) {
            genPt = genjet_pt[iGen];
            LogDebug << "genPt=" << genPt << " ";
          }
        }
        const auto rand = ( genPt < 0. ) ? rg.Gaus(0, ptRes) : -1.;
        LogDebug << "jet_pt_resolution: " << ptRes << ", rand: " << rand << std::endl;
        jerF_nom  = jetESmearFactor(jet_pt_nom, eOrig, genPt, ptRes, m_jetEResSF.getScaleFactor(jPar, Variation::NOMINAL), rand);
        jerF_down = jetESmearFactor(jet_pt_nom, eOrig, genPt, ptRes, m_jetEResSF.getScaleFactor(jPar, Variation::DOWN   ), rand);
        jerF_up   = jetESmearFactor(jet_pt_nom, eOrig, genPt, ptRes, m_jetEResSF.getScaleFactor(jPar, Variation::UP     ), rand);
        // LogDebug << "  scalefactors are NOMINAL=" << m_jetEResSF.getScaleFactor(jPar, Variation::NOMINAL) << ", DOWN=" << m_jetEResSF.getScaleFactor(jPar, Variation::DOWN) << ", UP=" << m_jetEResSF.getScaleFactor(jPar, Variation::UP) << std::endl;
        // LogDebug << "  smearfactors are NOMINAL=" << jerF_nom << ", DOWN=" << jerF_down << ", UP=" << jerF_up << std::endl;
        jet_pt_nom *= jerF_nom; // for consistency with jet case (used for split JER and HEM issue below)
      }
    }
    if ( ( jet_mask.empty() || jet_mask[i] ) && ( jet_pt_nomuL1L2L3 > m_unclEnThreshold ) && ( (jet_neEmEF[i]+jet_chEmEF[i]) < 0.9 ) ) {
      std::size_t iVar = 0;
      const auto jet_cosPhi = std::cos(jet_phi[i]);
      const auto jet_sinPhi = std::sin(jet_phi[i]);
      if ( ! ( m_doSmearing && m_isT1SmearedMET ) ) {
        out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, jet_pt_L1 - jet_pt_L1L2L3);           // nominal
      }
      auto jet_pt_L1p = jet_pt_L1; // with optional offset for JES uncertainty calculation if nominal is smeared
      if ( m_doSmearing ) {
        const auto dr_jernom  = jet_pt_L1 - jet_pt_L1L2L3*jerF_nom;
        if ( m_isT1SmearedMET ) {
          const auto dr_jerup   = jet_pt_L1 - jet_pt_L1L2L3*jerF_up;
          const auto dr_jerdown = jet_pt_L1 - jet_pt_L1L2L3*jerF_down;
          out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, dr_jernom);                         // smeared nominal
          if ( m_splitJER ) {
            const auto jerBin = jerSplitID(jet_pt_nom, jet_eta[i]);
            for ( int k{0}; k != 6; ++k ) {
              if ( jerBin == k ) { // vary
                out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, dr_jerup);                    // JER[k]-up
                out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, dr_jerdown);                  // JER[k]-down
              } else { // keep nominal
                out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, dr_jernom);                   // JER[k]-up
                out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, dr_jernom);                   // JER[k]-down
              }
            }
          } else {
            out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, dr_jerup);                        // JER-up
            out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, dr_jerdown);                      // JER-down
          }
          jet_pt_L1p += jet_pt_L1L2L3*(1.-jerF_nom); // offset for JES uncertainties, since the nominal is smeared
        } else {
          for ( std::size_t k{0}; k != ( m_splitJER ? 6 : 1 ); ++k ) {
            out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, dr_jernom);                     // JER[k]-up
            out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, dr_jernom);                     // JER[k]-down
          }
        }
      }
      if ( m_addHEM2018Issue ) {
        const auto delta = deltaHEM2018Issue(jet_pt_nom, jet_jetId[i], jet_phi[i], jet_eta[i]);
        out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, jet_pt_L1p - jet_pt_L1L2L3);           // up = nominal
        out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, jet_pt_L1p - jet_pt_L1L2L3*delta);     // down
      }
      for ( auto& jesUnc : m_jesUncSources ) {
        const auto delta = getUncertainty(jesUnc.second, { {JME::Binning::JetPt, jet_pt_L1L2L3}, {JME::Binning::JetEta, jet_eta[i]} }, true);
        out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, jet_pt_L1p - jet_pt_L1L2L3*(1+delta)); // JES_i-up
        out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, jet_pt_L1p - jet_pt_L1L2L3*(1-delta)); // JES_i-down
      }
#ifdef BAMBOO_JME_DEBUG
      assert(iVar+2 == out.size()); // last two are unclustered energy up and down
#endif
    }
  }
}

std::vector<std::string> Type1METVariationsCalculator::available(const std::string&) const
{
  if ((!m_jetCorrector) || (!m_jetCorrectorL1)) {
    throw std::runtime_error("The calculator is not fully configured (for MET variations both setJEC and setL1JEC need to be called)");
  }
  std::vector<std::string> products = { "nominal" };
  if ( m_doSmearing ) {
    if ( m_splitJER ) {
      for ( int i = 0; i != 6; ++i ) {
        products.emplace_back("jer"+std::to_string(i)+"up");
        products.emplace_back("jer"+std::to_string(i)+"down");
      }
    } else {
      products.emplace_back("jerup");
      products.emplace_back("jerdown");
    }
  }
  if ( m_addHEM2018Issue ) {
    products.emplace_back("jesHEMIssueup");
    products.emplace_back("jesHEMIssuedown");
  }
  for ( const auto& src : m_jesUncSources ) {
    products.emplace_back("jes"+src.first+"up");
    products.emplace_back("jes"+src.first+"down");
  }
  products.emplace_back("unclustEnup");
  products.emplace_back("unclustEndown");
  return products;
}

FixEE2017Type1METVariationsCalculator FixEE2017Type1METVariationsCalculator::create(
  const std::vector<std::string>& jecParams, const std::vector<std::string>& jecParamsL1,
  float unclEnThreshold,
  const std::vector<std::string>& jecParamsProd, const std::vector<std::string>& jecParamsL1Prod,
  const std::vector<std::pair<std::string,std::string>>& jesUncertainties,
  bool isT1SmearedMET,
  const std::string& ptResolution, const std::string& ptResolutionSF, bool splitJER,
  bool doGenMatch, float genMatch_maxDR, float genMatch_maxDPT)
{
  auto inst = FixEE2017Type1METVariationsCalculator();
  configureMETCalc_common(
    inst, jecParams, jecParamsL1, unclEnThreshold, jesUncertainties,
    isT1SmearedMET, ptResolution, ptResolutionSF, splitJER,
    doGenMatch, genMatch_maxDR, genMatch_maxDPT
  );
  inst.setJECProd(makeJCPList(jecParamsProd));
  inst.setL1JECProd(makeJCPList(jecParamsL1Prod));
  return inst;
}

void FixEE2017Type1METVariationsCalculator::setJECProd(const std::vector<JetCorrectorParameters>& jecParams)
{
  if ( ! jecParams.empty() ) {
    m_jetCorrectorProd = std::unique_ptr<FactorizedJetCorrectorCalculator>{new FactorizedJetCorrectorCalculator(jecParams)};
  }
}
void FixEE2017Type1METVariationsCalculator::setL1JECProd(const std::vector<JetCorrectorParameters>& jecParams)
{
  if ( ! jecParams.empty() ) {
    m_jetCorrectorL1Prod = std::unique_ptr<FactorizedJetCorrectorCalculator>{new FactorizedJetCorrectorCalculator(jecParams)};
  }
}

FixEE2017Type1METVariationsCalculator::result_t FixEE2017Type1METVariationsCalculator::produce(
    const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
    const p4compv_t& jet_rawcorr, const p4compv_t& jet_area,
    const p4compv_t& jet_muonSubtrFactor, const p4compv_t& jet_neEmEF, const p4compv_t& jet_chEmEF,
    const float rho,
    const std::uint32_t seed,
    const p4compv_t& genjet_pt, const p4compv_t& genjet_eta, const p4compv_t& genjet_phi, const p4compv_t& genjet_mass,
    const float rawmet_phi, const float rawmet_pt,
    const float met_unclustenupdx, const float met_unclustenupdy,
    const p4compv_t& lowptjet_rawpt, const p4compv_t& lowptjet_eta, const p4compv_t& lowptjet_phi, const p4compv_t& lowptjet_area,
    const p4compv_t& lowptjet_muonSubtrFactor, const p4compv_t& lowptjet_neEmEF, const p4compv_t& lowptjet_chEmEF,
    const float defmet_phi, const float defmet_pt, const float t1met_phi, const float t1met_pt
    ) const
{
  LogDebug << "JME:: hello from Type1METVariations produce with 2017 EE Fix. Got " << jet_pt.size() << " jets and " << lowptjet_rawpt.size() << " low-PT jets" << std::endl;
  assert(!m_addHEM2018Issue);
  const auto nVariations = 3+( m_doSmearing ? 2*( m_splitJER ? 6 : 1 ) : 0 )+2*m_jesUncSources.size(); // 1(nom)+2(unclust)+3(JER)+2*len(JES)
  p4compv_t lowptjet_zero(lowptjet_rawpt.size(), 0.);
  // the actual MET fix
  auto jet_mask = ROOT::VecOps::RVec<bool>(jet_pt.size(), true);
  auto lowptjet_mask = ROOT::VecOps::RVec<bool>(lowptjet_rawpt.size(), true);
  LogDebug << "JME:: First the (vetoed) jets in the noisy region" << std::endl;
  const auto offset_jets = calculateFixEE2017Offset(jet_mask,
      jet_pt, jet_eta, jet_phi, jet_mass,
      jet_rawcorr, jet_area, jet_muonSubtrFactor,
      rho);
  const auto offset_lowptjets = calculateFixEE2017Offset(lowptjet_mask,
      lowptjet_rawpt, lowptjet_eta, lowptjet_phi, lowptjet_zero,
      lowptjet_zero, lowptjet_area, lowptjet_muonSubtrFactor,
      rho);
  const auto delta_x_T1Jet  = offset_jets[0]+offset_lowptjets[0];
  const auto delta_y_T1Jet  = offset_jets[1]+offset_lowptjets[1];
  const auto delta_x_RawJet = offset_jets[2]+offset_lowptjets[2];
  const auto delta_y_RawJet = offset_jets[3]+offset_lowptjets[3];
  // default MET : add delta_T1Jet
  // unclustered EE : default MET (with delta T1Jet) - T1MET
  // correction for all MET variations: add delta_RawJet - unclustered EE
  const auto dx = delta_x_RawJet - ( defmet_pt*std::cos(defmet_phi) + delta_x_T1Jet - t1met_pt*std::cos(t1met_phi) );
  const auto dy = delta_y_RawJet - ( defmet_pt*std::sin(defmet_phi) + delta_y_T1Jet - t1met_pt*std::sin(t1met_phi) );
#ifdef BAMBOO_JME_DEBUG
  LogDebug << "JME:: T1      MET px=" << t1met_pt*std::cos(t1met_phi) << " py=" << t1met_pt*std::sin(t1met_phi) << std::endl;
  LogDebug << "JME:: default MET px=" << defmet_pt*std::cos(defmet_phi) << " py=" << defmet_pt*std::sin(defmet_phi) << std::endl;
  LogDebug << "JME:: raw     MET px=" << rawmet_pt*std::cos(rawmet_phi) << " py=" << rawmet_pt*std::sin(rawmet_phi) << std::endl;
  LogDebug << "JME:: deltas T1Jet x=" << delta_x_T1Jet << " y=" << delta_y_T1Jet << " RawJet x=" << delta_x_RawJet << " y= " << delta_y_RawJet << std::endl;
  LogDebug << "JME:: MET offset from jets in the noisy region: dx=" << dx << " and dy=" << dy << std::endl; // NO these are just minus the unclustered
#endif
  result_t out{nVariations, rawmet_pt*std::cos(rawmet_phi)+dx, rawmet_pt*std::sin(rawmet_phi)+dy};
  // usual variations, with jets that are in "unclustered EE" now vetoed
  auto& rg = getTRandom3(seed);
  // normal jets
  addVariations(out, jet_pt, jet_eta, jet_phi, jet_mass,
      jet_rawcorr, jet_area, jet_muonSubtrFactor, jet_neEmEF, jet_chEmEF, ROOT::VecOps::RVec<int>{},
      jet_mask, rho, genjet_pt, genjet_eta, genjet_phi, genjet_mass, rg);
  // low-PT jets
  addVariations(out, lowptjet_rawpt, lowptjet_eta, lowptjet_phi, lowptjet_zero,
      lowptjet_zero, lowptjet_area, lowptjet_muonSubtrFactor,
      ( lowptjet_neEmEF.empty() ? lowptjet_zero : lowptjet_neEmEF  ),
      ( lowptjet_chEmEF.empty() ? lowptjet_zero : lowptjet_chEmEF  ),
      ROOT::VecOps::RVec<int>{}, lowptjet_mask, rho,
      genjet_pt, genjet_eta, genjet_phi, genjet_mass, rg);
  // unclustered energy, based on nominal (0)
  out.setXY(nVariations-2, out.px(0)+met_unclustenupdx, out.py(0)+met_unclustenupdy);
  out.setXY(nVariations-1, out.px(0)-met_unclustenupdx, out.py(0)-met_unclustenupdy);

#ifdef BAMBOO_JME_DEBUG
  LogDebug << "JME:: returning " << out.size() << " modified METs" << std::endl;
  const auto varNames = available();
  assert(varNames.size() == nVariations);
  for ( std::size_t i{0}; i != nVariations; ++i ) {
    LogDebug << "JME:: MET_" << varNames[i] << ": PT=" << out.pt(i) << ", PHI=" << out.phi(i) << std::endl;
  }
#endif
  return out;
}

std::array<double,4> FixEE2017Type1METVariationsCalculator::calculateFixEE2017Offset(ROOT::VecOps::RVec<bool>& jet_mask,
    const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
    const p4compv_t& jet_rawcorr, const p4compv_t& jet_area, const p4compv_t& jet_muonSubtrFactor,
    const float rho
    ) const
{
  double delta_x_T1Jet{0.}, delta_y_T1Jet{0.};
  double delta_x_rawJet{0.}, delta_y_rawJet{0.};
  FactorizedJetCorrectorCalculator::VariableValues vals, valsL1;
  const auto nJets = jet_pt.size();
  for ( std::size_t i{0}; i != nJets; ++i ) {
    if ( ( 2.65 < std::abs(jet_eta[i]) ) && ( std::abs(jet_eta[i]) < 3.14 ) ) {
      const double jet_pt_raw = jet_pt[i]*(1-jet_rawcorr[i]);
      if ( jet_pt_raw < 50. ) {
        jet_mask[i] = false; // these are the jets to veto for the nominal variations
        // L1 and full (L1L2L3) JEC for muon-subtracted jet
        vals.setJetEta(jet_eta[i]);
        vals.setJetPt(jet_pt_raw);
        vals.setJetA(jet_area[i]);
        vals.setRho(rho);
        LogDebug << "Jet #" << i << " ETA=" << jet_eta[i] << ", PT_raw=" << jet_pt_raw << ", area=" << jet_area[i] << std::endl;
        auto jecL1L2L3 = ( m_jetCorrectorProd ? m_jetCorrectorProd : m_jetCorrector )->getCorrection(vals);
        if ( jecL1L2L3 <= 0. ) { jecL1L2L3 = 1.; }
        valsL1.setJetEta(jet_eta[i]);
        valsL1.setJetPt(jet_pt_raw);
        valsL1.setJetA(jet_area[i]);
        valsL1.setRho(rho);
        auto jecL1 = ( m_jetCorrectorL1Prod ? m_jetCorrectorL1Prod : m_jetCorrectorL1 )->getCorrection(valsL1);
        if ( jecL1     <= 0. ) { jecL1     = 1.; }
        const auto jet_pt_raw_nomu = jet_pt_raw*(1-jet_muonSubtrFactor[i]);
        const auto muon_pt = jet_pt_raw*jet_muonSubtrFactor[i];
        const auto jet_pt_nomuL1L2L3 = jet_pt_raw_nomu*jecL1L2L3;
        const auto jet_pt_nomuL1     = jet_pt_raw_nomu*jecL1;
        const auto jet_pt_L1L2L3 = jet_pt_nomuL1L2L3 + muon_pt;
        const auto jet_pt_L1     = jet_pt_nomuL1     + muon_pt;
        LogDebug << "jecL1L2L3=" << jecL1L2L3 << ", jecL1=" << jecL1 << "; PT_L1L2L3=" << jet_pt_L1L2L3 << ", PT_L1=" << jet_pt_L1 << ", PT_mu=" << muon_pt << std::endl;
        if ( jet_pt_nomuL1L2L3 > m_unclEnThreshold ) {
          const auto cosphi = std::cos(jet_phi[i]);
          const auto sinphi = std::sin(jet_phi[i]);
          delta_x_T1Jet += (jet_pt_L1L2L3-jet_pt_L1+jet_pt_raw)*cosphi;
          delta_y_T1Jet += (jet_pt_L1L2L3-jet_pt_L1+jet_pt_raw)*sinphi;
          delta_x_rawJet += jet_pt_raw*cosphi;
          delta_y_rawJet += jet_pt_raw*sinphi;
        }
      }
    }
  }
  return {{ delta_x_T1Jet, delta_y_T1Jet, delta_x_rawJet, delta_y_rawJet }};
}
