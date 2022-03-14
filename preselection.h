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


//int download_file(const char* url, const char* file_name)
//{
//  CURL* easyhandle = curl_easy_init();
//  curl_easy_setopt( easyhandle, CURLOPT_URL, url ) ;
//  FILE* file = fopen( file_name, "w");
//  curl_easy_setopt( easyhandle, CURLOPT_WRITEDATA, file) ;
//  curl_easy_perform( easyhandle );
//  curl_easy_cleanup( easyhandle );
//  fclose(file);
//  return 0;
//}

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

string path = "https://ttedesch.web.cern.ch/ttedesch/nanoAOD-tools/python/postprocessing/data/leptonSF/";

std::vector<std::string> mu_f_2016{path + "Mu_RunBCDEFGH_SF_ID_2016_syst.root"};
std::vector<std::string> el_f_2016{path + "EGM2D_RECO_SF_2016.root", path + "2016LegacyReReco_ElectronMVA90noiso_Fall17V2.root"};
std::vector<std::string> mu_h_2016{"NUM_TightID_DEN_genTracks_eta_pt"};
std::vector<std::string> el_h_2016{"EGamma_SF2D", "EGamma_SF2D"};
std::vector<std::string> mu_f_2017{path + "Muon_RunBCDEF_SF_ID_2017.root", path + "Muon_RunBCDEF_SF_ISO_2017.root"};
std::vector<std::string> el_f_2017{path + "EGM2D_2017_passingRECO_highEt.root", path + "Electron_MVA90_2017.root"};
std::vector<std::string> mu_h_2017{"NUM_TightID_DEN_genTracks_pt_abseta", "NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta"};
std::vector<std::string> el_h_2017{"EGamma_SF2D", "EGamma_SF2D"};
std::vector<std::string> mu_f_2018{path + "Muon_RunBCDEF_SF_ID_2018.root", path + "Muon_RunBCDEF_SF_ISO_2018.root"};
std::vector<std::string> el_f_2018{path + "EGM2D_passingRECO_2018All.root", path + "2018_ElectronMVA90Iso.root"};
std::vector<std::string> mu_h_2018{"NUM_TightID_DEN_TrackerMuons_pt_abseta", "NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta"};
std::vector<std::string> el_h_2018{"EGamma_SF2D", "EGamma_SF2D"};

/*
{} 
mu_h_2016, el_f_2016, el_h_2016;
std::vector<std::string> mu_f_2017, mu_h_2017, el_f_2017, el_h_2017;
std::vector<std::string> mu_f_2018, mu_h_2018, el_f_2018, el_h_2018;



mu_f_2016.push_back(path + "Mu_RunBCDEFGH_SF_ID_2016_syst.root");
mu_h_2016.push_back("NUM_TightID_DEN_genTracks_eta_pt");

el_f_2016.push_back(path + "EGM2D_RECO_SF_2016.root");
el_f_2016.push_back(path + "2016LegacyReReco_ElectronMVA90noiso_Fall17V2.root");
el_h_2016.push_back("EGamma_SF2D");
el_h_2016.push_back("EGamma_SF2D");

mu_f_2017.push_back(path + "Muon_RunBCDEF_SF_ID_2017.root");
mu_f_2017.push_back(path + "Muon_RunBCDEF_SF_ISO_2017.root");
mu_h_2017.push_back("NUM_TightID_DEN_genTracks_pt_abseta");
mu_h_2017.push_back("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta");

el_f_2017.push_back(path + "EGM2D_2017_passingRECO_highEt.root");
el_f_2017.push_back(path + "Electron_MVA90_2017.root");
el_h_2017.push_back("EGamma_SF2D");
el_h_2017.push_back("EGamma_SF2D");

mu_f_2018.push_back(path + "Muon_RunABCD_SF_ID_2018.root");
mu_f_2018.push_back(path + "Muon_RunABCD_SF_ISO_2018.root");
mu_h_2018.push_back("NUM_TightID_DEN_TrackerMuons_pt_abseta");
mu_h_2018.push_back("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta");

el_f_2018.push_back(path + "EGM2D_passingRECO_2018All.root");
el_f_2018.push_back(path + "2018_ElectronMVA90Iso.root");
el_h_2018.push_back("EGamma_SF2D");
el_h_2018.push_back("EGamma_SF2D");

*/
//LeptonEfficiencyCorrector worker_mu;
//LeptonEfficiencyCorrector worker_el;

LeptonEfficiencyCorrector worker_mu_2016(mu_f_2016, mu_h_2016);
LeptonEfficiencyCorrector worker_el_2016(el_f_2016, el_h_2016);
LeptonEfficiencyCorrector worker_mu_2017(mu_f_2017, mu_h_2017);
LeptonEfficiencyCorrector worker_el_2017(el_f_2017, el_h_2017);
LeptonEfficiencyCorrector worker_mu_2018(mu_f_2018, mu_h_2018);
LeptonEfficiencyCorrector worker_el_2018(el_f_2018, el_h_2018);

/*
RVec<RVec<float>> LepSF(rvec_f Electron_pt, rvec_f Electron_eta, rvec_i Electron_pdgId, rvec_f Muon_pt, rvec_f Muon_eta, rvec_i Muon_pdgId, string Year){
    
    LeptonEfficiencyCorrector worker_mu, worker_el;

    if (Year == "2016") {
        worker_mu = worker_mu_2016;
        worker_el = worker_mu_2016;
    }
    if (Year == "2017"){
        worker_mu = worker_mu_2017;
        worker_el = worker_mu_2017;
    }
    
    if (Year == "2018"){
        worker_mu = worker_mu_2018;
        worker_el = worker_mu_2018;
    }
    
    RVec<float> sf_mu(Muon_pt.size());
    RVec<float> sferr_mu(Muon_pt.size()); 
    RVec<float> sf_el(Electron_pt.size());
    RVec<float> sferr_el(Electron_pt.size());
    
    for (size_t j = 0; j < Electron_pt.size(); j++) sf_el[j] =  worker_el.getSF(Electron_pdgId[j], Electron_pt[j], Electron_eta[j]);
    for (size_t j = 0; j < Electron_pt.size(); j++) sferr_el[j] =  worker_el.getSFErr(Electron_pdgId[j], Electron_pt[j], Electron_eta[j]);
     
    for (size_t j = 0; j < Muon_pt.size(); j++) sf_mu[j] =  worker_el.getSF(Muon_pdgId[j], Muon_pt[j], Muon_eta[j]);
    for (size_t j = 0; j < Muon_pt.size(); j++) sferr_mu[j] =  worker_el.getSFErr(Muon_pdgId[j], Muon_pt[j], Muon_eta[j]);
        
    RVec<RVec<float>> result;
    
    result.emplace_back(sf_mu);
    result.emplace_back(sf_el);
    
    RVec<float> Muon_effSF_errUp(Muon_pt.size());
    for (size_t j = 0; j < Muon_pt.size(); j++) Muon_effSF_errUp[j] = sferr_mu[j] + sf_mu[j];
    result.emplace_back(Muon_effSF_errUp);
    
    RVec<float> Electron_effSF_errUp(Electron_pt.size());
    for (size_t j = 0; j < Electron_pt.size(); j++) Electron_effSF_errUp[j] = sferr_el[j] + sf_el[j];
    result.emplace_back(Electron_effSF_errUp);
    
    RVec<float> Muon_effSF_errDown(Muon_pt.size());
    for (size_t j = 0; j < Muon_pt.size(); j++) Muon_effSF_errDown[j] = sferr_mu[j] - sf_mu[j];
    result.emplace_back(Muon_effSF_errDown);
    
    RVec<float> Electron_effSF_errDown(Electron_pt.size());
    for (size_t j = 0; j < Electron_pt.size(); j++) Electron_effSF_errUp[j] = sferr_el[j] - sf_el[j];
    result.emplace_back(Electron_effSF_errDown);
    
    return result;
}
*/

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

string remote_storage = "https://ttedesch.web.cern.ch/ttedesch/nanoAOD-tools/python/postprocessing/data/";
void set_null_directory(TH2F *histo){
    histo->SetDirectory(NULL);
}

/*
TFile *pufile_data2016 = TFile::Open(TString(remote_storage) + TString("pileup/PileupData_GoldenJSON_Full2016.root"));
TH1 *histo_target_2016 = (TH1*)pufile_data2016->Get("pileup");
TH1 *histo_target_2016_plus = (TH1*)pufile_data2016->Get("pileup_plus");
TH1 *histo_target_2016_minus = (TH1*)pufile_data2016->Get("pileup_minus");
//histo_target_2016->SetDirectory(NULL);
//histo_target_2016_plus->SetDirectory(NULL);
//histo_target_2016_minus->SetDirectory(NULL);
TFile *pufile_mc2016 = TFile::Open(TString(remote_storage) + TString("pileup/pileup_profile_Summer16.root"));
TH1 *histo_2016 = (TH1*)pufile_mc2016->Get("pileup");
TH1 *histo_2016_plus = (TH1*)pufile_mc2016->Get("pileup_plus");
TH1 *histo_2016_minus = (TH1*)pufile_mc2016->Get("pileup_minus");
//histo_2016->SetDirectory(NULL);
//histo_2016_plus->SetDirectory(NULL);
//histo_2016_minus->SetDirectory(NULL);

TFile *pufile_data2017 = TFile::Open(TString(remote_storage) + TString("pileup/PileupHistogram-goldenJSON-13tev-2017-99bins_withVar.root"));
TH1 *histo_target_2017 = (TH1*)pufile_data2017->Get("pileup");
TH1 *histo_target_2017_plus = (TH1*)pufile_data2017->Get("pileup_plus");
TH1 *histo_target_2017_minus = (TH1*)pufile_data2017->Get("pileup_minus");
//histo_target_2017->SetDirectory(NULL);
//histo_target_2017_plus->SetDirectory(NULL);
//histo_target_2017_minus->SetDirectory(NULL);
TFile *pufile_mc2017 = TFile::Open(TString(remote_storage) + TString("pileup/mcPileup2017.root"));
TH1 *histo_2017 = (TH1*)pufile_mc2017->Get("pileup");
TH1 *histo_2017_plus = (TH1*)pufile_mc2017->Get("pileup_plus");
TH1 *histo_2017_minus = (TH1*)pufile_mc2017->Get("pileup_minus");
//histo_2017->SetDirectory(NULL);
//histo_2017_plus->SetDirectory(NULL);
//histo_2017_minus->SetDirectory(NULL);

TFile *pufile_data2018 = TFile::Open(TString(remote_storage) + TString("pileup/PileupHistogram-goldenJSON-13tev-2018-100bins_withVar.root"));
TH1 *histo_target_2018 = (TH1*)pufile_data2018->Get("pileup");
TH1 *histo_target_2018_plus = (TH1*)pufile_data2018->Get("pileup_plus");
TH1 *histo_target_2018_minus = (TH1*)pufile_data2018->Get("pileup_minus");
//histo_target_2018->SetDirectory(NULL);
//histo_target_2018_plus->SetDirectory(NULL);
//histo_target_2018_minus->SetDirectory(NULL);
TFile *pufile_mc2018 = TFile::Open(TString(remote_storage) + TString("pileup/mcPileup2018.root"));
TH1 *histo_2018 = (TH1*)pufile_mc2018->Get("pileup");
TH1 *histo_2018_plus = (TH1*)pufile_mc2018->Get("pileup_plus");
TH1 *histo_2018_minus = (TH1*)pufile_mc2018->Get("pileup_minus");
//histo_2018->SetDirectory(NULL);
//histo_2018_plus->SetDirectory(NULL);
//histo_2018_minus->SetDirectory(NULL);
//close files??

RVec<float> puWeight(string Year, int Pileup_nTrueInt){
    RVec<float> result(3);
    if (Year == "2016"){ 
        histo_2016->SetDirectory(NULL);
        histo_2016_plus->SetDirectory(NULL);
        histo_2016_minus->SetDirectory(NULL);
        histo_target_2016->SetDirectory(NULL);
        histo_target_2016_plus->SetDirectory(NULL);
        histo_target_2016_minus->SetDirectory(NULL);
        WeightCalculatorFromHistogram worker_2016(histo_2016, histo_target_2016, true, true, false);
        WeightCalculatorFromHistogram worker_2016_plus(histo_2016, histo_target_2016, true, true, false);
        WeightCalculatorFromHistogram worker_2016_minus(histo_2016, histo_target_2016, true, true, false);
        if(Pileup_nTrueInt < histo_2016->GetNbinsX()) result[0] = worker_2016.getWeight(Pileup_nTrueInt);
        else result[0] = 1;
        if(Pileup_nTrueInt < histo_2016_plus->GetNbinsX()) result[1] = worker_2016_plus.getWeight(Pileup_nTrueInt);
        else result[1] = 1;
        if(Pileup_nTrueInt < histo_2016_minus->GetNbinsX()) result[2] = worker_2016_minus.getWeight(Pileup_nTrueInt);
        else result[2] = 1;
    }
    else if (Year == "2017"){
        histo_2017->SetDirectory(NULL);
        histo_2017_plus->SetDirectory(NULL);
        histo_2017_minus->SetDirectory(NULL);
        histo_target_2017->SetDirectory(NULL);
        histo_target_2017_plus->SetDirectory(NULL);
        histo_target_2017_minus->SetDirectory(NULL);
        WeightCalculatorFromHistogram worker_2017(histo_2017, histo_target_2017, true, true, false);
        WeightCalculatorFromHistogram worker_2017_plus(histo_2017, histo_target_2017, true, true, false);
        WeightCalculatorFromHistogram worker_2017_minus(histo_2017, histo_target_2017, true, true, false);
        if(Pileup_nTrueInt < histo_2017->GetNbinsX()) result[0] = worker_2017.getWeight(Pileup_nTrueInt);
        else result[0] = 1;
        if(Pileup_nTrueInt < histo_2017_plus->GetNbinsX()) result[1] = worker_2017_plus.getWeight(Pileup_nTrueInt);
        else result[1] = 1;
        if(Pileup_nTrueInt < histo_2017_minus->GetNbinsX()) result[2] = worker_2017_minus.getWeight(Pileup_nTrueInt);
        else result[2] = 1;
    }
    else if (Year == "2018"){
        histo_2018->SetDirectory(NULL);
        histo_2018_plus->SetDirectory(NULL);
        histo_2018_minus->SetDirectory(NULL);
        histo_target_2018->SetDirectory(NULL);
        histo_target_2018_plus->SetDirectory(NULL);
        histo_target_2018_minus->SetDirectory(NULL);
        WeightCalculatorFromHistogram worker_2018(histo_2018, histo_target_2018, true, true, false);
        WeightCalculatorFromHistogram worker_2018_plus(histo_2018, histo_target_2018, true, true, false);
        WeightCalculatorFromHistogram worker_2018_minus(histo_2018, histo_target_2018, true, true, false);
        if(Pileup_nTrueInt < histo_2018->GetNbinsX()) result[0] = worker_2018.getWeight(Pileup_nTrueInt);
        else result[0] = 1;
        if(Pileup_nTrueInt < histo_2018_plus->GetNbinsX()) result[1] = worker_2018_plus.getWeight(Pileup_nTrueInt);
        else result[1] = 1;
        if(Pileup_nTrueInt < histo_2018_minus->GetNbinsX()) result[2] = worker_2018_minus.getWeight(Pileup_nTrueInt);
        else result[2] = 1;
    }
    return result;    
}
*/


TFile *pufile_data2016 = TFile::Open(TString(remote_storage) + TString("pileup/PileupData_GoldenJSON_Full2016.root"));
TH1 *histo_target_2016 = (TH1*)pufile_data2016->Get("pileup");
TH1 *histo_target_2016_plus = (TH1*)pufile_data2016->Get("pileup_plus");
TH1 *histo_target_2016_minus = (TH1*)pufile_data2016->Get("pileup_minus");
TFile *pufile_mc2016 = TFile::Open(TString(remote_storage) + TString("pileup/pileup_profile_Summer16.root"));
TH1 *histo_2016 = (TH1*)pufile_mc2016->Get("pu_mc");

TFile *pufile_data2017 = TFile::Open(TString(remote_storage) + TString("pileup/PileupHistogram-goldenJSON-13tev-2017-99bins_withVar.root"));
TH1 *histo_target_2017 = (TH1*)pufile_data2017->Get("pileup");
TH1 *histo_target_2017_plus = (TH1*)pufile_data2017->Get("pileup_plus");
TH1 *histo_target_2017_minus = (TH1*)pufile_data2017->Get("pileup_minus");
TFile *pufile_mc2017 = TFile::Open(TString(remote_storage) + TString("pileup/mcPileup2017.root"));
TH1 *histo_2017 = (TH1*)pufile_mc2017->Get("pu_mc");

TFile *pufile_data2018 = TFile::Open(TString(remote_storage) + TString("pileup/PileupHistogram-goldenJSON-13tev-2018-100bins_withVar.root"));
TH1 *histo_target_2018 = (TH1*)pufile_data2018->Get("pileup");
TH1 *histo_target_2018_plus = (TH1*)pufile_data2018->Get("pileup_plus");
TH1 *histo_target_2018_minus = (TH1*)pufile_data2018->Get("pileup_minus");
TFile *pufile_mc2018 = TFile::Open(TString(remote_storage) + TString("pileup/mcPileup2018.root"));
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
string remote_storage_ = "https://ttedesch.web.cern.ch/ttedesch/nanoAOD-tools/";
TFile *L1PrefiringMaps = TFile::Open(TString(remote_storage_) + TString("data/prefire_maps/L1PrefiringMaps.root"));
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

    vector<int> v{0,-1,1};
    for (int j = 0; j < v.size(); j++){
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

//int a = download_file("https://ttedesch.web.cern.ch/ttedesch/nanoAOD-tools/python/postprocessing/data/roccor.Run2.v3/RoccoR2016.txt", "RoccoR2016.txt"); 
//int b = download_file("https://ttedesch.web.cern.ch/ttedesch/nanoAOD-tools/python/postprocessing/data/roccor.Run2.v3/RoccoR2017.txt", "RoccoR2017.txt"); 
//int c = download_file("https://ttedesch.web.cern.ch/ttedesch/nanoAOD-tools/python/postprocessing/data/roccor.Run2.v3/RoccoR2018.txt", "RoccoR2018.txt"); 
//a = download_file("https://ttedesch.web.cern.ch/ttedesch/nanoAOD-tools/python/postprocessing/data/roccor.Run2.v3/RoccoR2016aUL.txt", "RoccoR2016aUL.txt"); 
//a = download_file("https://ttedesch.web.cern.ch/ttedesch/nanoAOD-tools/python/postprocessing/data/roccor.Run2.v3/RoccoR2016bUL.txt", "RoccoR2016bUL.txt"); 
//a = download_file("https://ttedesch.web.cern.ch/ttedesch/nanoAOD-tools/python/postprocessing/data/roccor.Run2.v3/RoccoR2017UL.txt", "RoccoR2017UL.txt"); 
//a = download_file("https://ttedesch.web.cern.ch/ttedesch/nanoAOD-tools/python/postprocessing/data/roccor.Run2.v3/RoccoR2018UL.txt", "RoccoR2018UL.txt");  
//int d = download_file("https://ttedesch.web.cern.ch/ttedesch/nanoAOD-tools/python/postprocessing/data/roccor.Run2.v3/RoccoR.cc", "RoccoR.cc");
//int f = download_file("https://ttedesch.web.cern.ch/ttedesch/nanoAOD-tools/python/postprocessing/data/roccor.Run2.v3/RoccoR.h", "RoccoR.h");
//R__LOAD_LIBRARY(./RoccoR.cc);

RoccoR roccor_2016("RoccoR2016.txt");
RoccoR roccor_2017("RoccoR2017.txt");
RoccoR roccor_2018("RoccoR2018.txt");
//RoccoR roccor_2016aUL("RoccoR2016aUL.txt");
//RoccoR roccor_2016bUL("RoccoR2016bUL.txt");
//RoccoR roccor_2017UL("RoccoR2017UL.txt");
//RoccoR roccor_2018UL("RoccoR2018UL.txt");

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0, 1);//uniform distribution between 0 and 1

/*
RVec<RVec<float>> muonScaleRes(rvec_f Muon_pt, rvec_f Muon_eta, rvec_f Muon_phi, rvec_i Muon_charge, rvec_i Muon_nTrackerLayers, rvec_i Muon_genPartIdx, rvec_f GenPart_pt, string era){
    RVec<RVec<float>> result;
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
    
    result.emplace_back(pt_corr);
    result.emplace_back(pt_corr_up);
    result.emplace_back(pt_corr_down);
    
    return result;
}
*/

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



//download_file("https://ttedesch.web.cern.ch/ttedesch/nanoAOD-tools/data/btagSF/DeepJet_2016LegacySF_V1.csv", "DeepJet_2016LegacySF_V1.csv"); 
//download_file("https://ttedesch.web.cern.ch/ttedesch/nanoAOD-tools/data/btagSF/DeepFlavour_94XSF_V3_B_F.csv", "DeepFlavour_94XSF_V3_B_F.csv"); download_file("https://ttedesch.web.cern.ch/ttedesch/nanoAOD-tools/data/btagSF/DeepJet_102XSF_V1.csv", "DeepJet_102XSF_V1.csv"); download_file("https://ttedesch.web.cern.ch/ttedesch/nanoAOD-tools/data/btagSF/DeepJet_106XUL16preVFPSF_v1_new.csv", "DeepJet_106XUL16preVFPSF_v1_new.csv"); download_file("https://ttedesch.web.cern.ch/ttedesch/nanoAOD-tools/data/btagSF/DeepJet_106XUL16postVFPSF_v2_new.csv", "DeepJet_106XUL16postVFPSF_v2_new.csv"); download_file("https://ttedesch.web.cern.ch/ttedesch/nanoAOD-tools/data/btagSF/DeepJet_106XUL17_v3_new.csv", "DeepJet_106XUL17_v3_new.csv"); download_file("https://ttedesch.web.cern.ch/ttedesch/nanoAOD-tools/data/btagSF/DeepJet_106XUL18_v2_new.csv", "DeepJet_106XUL18_v2_new.csv"); 



//BTagCalibration calibration_Legacy2016("deepjet", "DeepJet_2016LegacySF_V1.csv");
BTagCalibration calibration_2017("deepjet", "DeepFlavour_94XSF_V3_B_F.csv");
//BTagCalibration calibration_2018("deepjet", "DeepJet_102XSF_V1.csv");
//BTagCalibration calibration_UL2016APV("deepjet", "DeepJet_106XUL16preVFPSF_v1_new.csv");
//BTagCalibration calibration_UL2016("deepjet", "DeepJet_106XUL16postVFPSF_v2_new.csv");
//BTagCalibration calibration_UL2017("deepjet", "DeepJet_106XUL17_v3_new.csv");
//BTagCalibration calibration_UL2018("deepjet", "DeepJet_106XUL18_v2_new.csv");



/*
vector<string> v_systs{"up", "down"};
BTagCalibrationReader reader_0_Legacy2016(0,"central",v_systs); //0 is wp_btv
BTagCalibrationReader reader_1_Legacy2016(1,"central",v_systs);
BTagCalibrationReader reader_2_Legacy2016(2,"central",v_systs); 
BTagCalibrationReader reader_0_2017(0,"central",v_systs); //0 is wp_btv
BTagCalibrationReader reader_1_2017(1,"central",v_systs);
BTagCalibrationReader reader_2_2017(2,"central",v_systs); 
BTagCalibrationReader reader_0_2018(0,"central",v_systs); //0 is wp_btv
BTagCalibrationReader reader_1_2018(1,"central",v_systs);
BTagCalibrationReader reader_2_2018(2,"central",v_systs); 
BTagCalibrationReader reader_0_UL2016APV(0,"central",v_systs); //0 is wp_btv
BTagCalibrationReader reader_1_UL2016APV(1,"central",v_systs);
BTagCalibrationReader reader_2_UL2016APV(2,"central",v_systs); 
BTagCalibrationReader reader_0_UL2016(0,"central",v_systs); //0 is wp_btv
BTagCalibrationReader reader_1_UL2016(1,"central",v_systs);
BTagCalibrationReader reader_2_UL2016(2,"central",v_systs); 
BTagCalibrationReader reader_0_UL2017(0,"central",v_systs); //0 is wp_btv
BTagCalibrationReader reader_1_UL2017(1,"central",v_systs);
BTagCalibrationReader reader_2_UL2017(2,"central",v_systs); 
BTagCalibrationReader reader_0_UL2018(0,"central",v_systs); //0 is wp_btv
BTagCalibrationReader reader_1_UL2018(1,"central",v_systs);
BTagCalibrationReader reader_2_UL2018(2,"central",v_systs); 
*/

//BTagEntry btag_dummy;
//BTagEntry::OperatingPoint a(OP_LOOSE);


vector<string> v_systs{"up", "down"};
//BTagCalibrationReader reader_0_Legacy2016(BTagEntry::OP_LOOSE,"central",v_systs); //0 is wp_btv
//BTagCalibrationReader reader_1_Legacy2016(BTagEntry::OP_MEDIUM,"central",v_systs);
//BTagCalibrationReader reader_2_Legacy2016(BTagEntry::OP_TIGHT,"central",v_systs); 
BTagCalibrationReader reader_0_2017(BTagEntry::OP_LOOSE,"central",v_systs); //0 is wp_btv
BTagCalibrationReader reader_1_2017(BTagEntry::OP_MEDIUM,"central",v_systs);
BTagCalibrationReader reader_2_2017(BTagEntry::OP_TIGHT,"central",v_systs); 
//BTagCalibrationReader reader_0_2018(BTagEntry::OP_LOOSE,"central",v_systs); //0 is wp_btv
//BTagCalibrationReader reader_1_2018(BTagEntry::OP_MEDIUM,"central",v_systs);
//BTagCalibrationReader reader_2_2018(BTagEntry::OP_TIGHT,"central",v_systs); 
//BTagCalibrationReader reader_0_UL2016APV(BTagEntry::OP_LOOSE,"central",v_systs); //0 is wp_btv
//BTagCalibrationReader reader_1_UL2016APV(BTagEntry::OP_MEDIUM,"central",v_systs);
//BTagCalibrationReader reader_2_UL2016APV(BTagEntry::OP_TIGHT,"central",v_systs); 
//BTagCalibrationReader reader_0_UL2016(BTagEntry::OP_LOOSE,"central",v_systs); //0 is wp_btv
//BTagCalibrationReader reader_1_UL2016(BTagEntry::OP_MEDIUM,"central",v_systs);
//BTagCalibrationReader reader_2_UL2016(BTagEntry::OP_TIGHT,"central",v_systs); 
//BTagCalibrationReader reader_0_UL2017(BTagEntry::OP_LOOSE,"central",v_systs); //0 is wp_btv
//BTagCalibrationReader reader_1_UL2017(BTagEntry::OP_MEDIUM,"central",v_systs);
//BTagCalibrationReader reader_2_UL2017(BTagEntry::OP_TIGHT,"central",v_systs); 
//BTagCalibrationReader reader_0_UL2018(BTagEntry::OP_LOOSE,"central",v_systs); //0 is wp_btv
//BTagCalibrationReader reader_1_UL2018(BTagEntry::OP_MEDIUM,"central",v_systs);
//BTagCalibrationReader reader_2_UL2018(BTagEntry::OP_TIGHT,"central",v_systs); 


/*
reader_0_Legacy2016.load(calibration_Legacy2016, 0, "comb"); //0 is flavor_btv
reader_0_Legacy2016.load(calibration_Legacy2016, 1, "comb");
reader_0_Legacy2016.load(calibration_Legacy2016, 2, "incl");
reader_1_Legacy2016.load(calibration_Legacy2016, 0, "comb"); //0 is flavor_btv
reader_1_Legacy2016.load(calibration_Legacy2016, 1, "comb");
reader_1_Legacy2016.load(calibration_Legacy2016, 2, "incl");
reader_2_Legacy2016.load(calibration_Legacy2016, 0, "comb"); //0 is flavor_btv
reader_2_Legacy2016.load(calibration_Legacy2016, 1, "comb");
reader_2_Legacy2016.load(calibration_Legacy2016, 2, "incl");

reader_0_2017.load(calibration_2017, 0, "comb"); //0 is flavor_btv
reader_0_2017.load(calibration_2017, 1, "comb");
reader_0_2017.load(calibration_2017, 2, "incl");
reader_1_2017.load(calibration_2017, 0, "comb"); //0 is flavor_btv
reader_1_2017.load(calibration_2017, 1, "comb");
reader_1_2017.load(calibration_2017, 2, "incl");
reader_2_2017.load(calibration_2017, 0, "comb"); //0 is flavor_btv
reader_2_2017.load(calibration_2017, 1, "comb");
reader_2_2017.load(calibration_2017, 2, "incl");
    
reader_0_2018.load(calibration_2018, 0, "comb"); //0 is flavor_btv
reader_0_2018.load(calibration_2018, 1, "comb");
reader_0_2018.load(calibration_2018, 2, "incl");
reader_1_2018.load(calibration_2018, 0, "comb"); //0 is flavor_btv
reader_1_2018.load(calibration_2018, 1, "comb");
reader_1_2018.load(calibration_2018, 2, "incl");
reader_2_2018.load(calibration_2018, 0, "comb"); //0 is flavor_btv
reader_2_2018.load(calibration_2018, 1, "comb");
reader_2_2018.load(calibration_2018, 2, "incl");
    
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