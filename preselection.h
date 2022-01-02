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
  float x = abs(pdgid)==13 ? pt : eta;
  float y = abs(pdgid)==13 ? fabs(eta) : pt;
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


RVec<RVec<float>> LepSF(rvec_f Electron_pt, rvec_f Electron_eta, rvec_i Electron_pdgId, rvec_f Muon_pt, rvec_f Muon_eta, rvec_i Muon_pdgId, Int_t Year){
    
    LeptonEfficiencyCorrector worker_mu, worker_el;

    if (Year == 2016) {
        worker_mu = worker_mu_2016;
        worker_el = worker_mu_2016;
    }
    if (Year == 2017){
        worker_mu = worker_mu_2017;
        worker_el = worker_mu_2017;
    }
    
    if (Year == 2018){
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

string remote_storage = "https://ttedesch.web.cern.ch/ttedesch/nanoAOD-tools/python/postprocessing/data/";
void set_null_directory(TH2F *histo){
    histo->SetDirectory(NULL);
}



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

RVec<float> puWeight(int Year, int Pileup_nTrueInt){
    RVec<float> result(3);
    if (Year == 2016){ 
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
    else if (Year == 2017){
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
    else if (Year == 2018){
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