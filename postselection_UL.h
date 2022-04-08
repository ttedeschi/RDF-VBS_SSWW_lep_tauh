/// Header file with functions needed to execute the Python version of
/// postselection step of the analysis. The header is declared to the
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

//values for cuts and constant 

const size_t ONLYELE=1;
const size_t ONLYMU=0;

const float PT_CUT_MU=  30;
const float ETA_CUT_MU= 2.4;
const float ISO_CUT_MU= 0.15;

//const size_t PT_CUT_ELE=  35;
const float PT_CUT_ELE_UL2016=  30;
const float PT_CUT_ELE_UL2017=  38;
const float PT_CUT_ELE_UL2018=  35;

const float ETA_CUT_ELE= 2.4;
const float ISO_CUT_ELE= 0.08;

const float REL_ISO_CUT_LEP_VETO_ELE=   0.2;
const float PT_CUT_LEP_VETO_ELE=        15;
const float ETA_CUT_LEP_VETO_ELE=       2.4;

const float REL_ISO_CUT_LEP_VETO_MU=    0.4;
const float PT_CUT_LEP_VETO_MU=         10;
const float ETA_CUT_LEP_VETO_MU=        2.4;

const float DR_OVERLAP_CONE_TAU=        0.5;
const float DR_OVERLAP_CONE_OTHER=      0.4;

const float PT_CUT_JET= 30;
const float ETA_CUT_JET=5;

const float DELTAETA_JJ_CUT=2.5;

const float BTAG_PT_CUT =   30;
const float BTAG_ETA_CUT=   5;
string BTAG_ALGO   =   "DeepFlv";
string BTAG_WP     =   "M";
const float BTAG_WP_VALUE = 0.3033;
const float BTAG_WP_VALUE_LOOSE = 0.0521;

const size_t ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_ELE_UL2016APV = 8;
const size_t ID_TAU_RECO_DEEPTAU_VSJET_VETO_ELE_UL2016APV = 8;
const size_t ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_MU_UL2016APV = 4;
const size_t ID_TAU_RECO_DEEPTAU_VSJET_VETO_MU_UL2016APV = 4;
const size_t ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_ELE_UL2016 = 4 ;
const size_t ID_TAU_RECO_DEEPTAU_VSJET_VETO_ELE_UL2016 = 4;
const size_t ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_MU_UL2016 = 8;
const size_t ID_TAU_RECO_DEEPTAU_VSJET_VETO_MU_UL2016 = 8;
const size_t ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_ELE_UL2017 = 16;
const size_t ID_TAU_RECO_DEEPTAU_VSJET_VETO_ELE_UL2017 = 16;
const size_t ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_MU_UL2017 = 8;
const size_t ID_TAU_RECO_DEEPTAU_VSJET_VETO_MU_UL2017 = 8;
const size_t ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_ELE_UL2018 = 16;
const size_t ID_TAU_RECO_DEEPTAU_VSJET_VETO_ELE_UL2018 = 16;
const size_t ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_MU_UL2018 = 16;
const size_t ID_TAU_RECO_DEEPTAU_VSJET_VETO_MU_UL2018 = 16;

//const size_t ID_TAU_RECO_DEEPTAU_VSJET_VETO_ELE = 16;  //DUMMY!!!!!!!!
//const size_t ID_TAU_RECO_DEEPTAU_VSJET_VETO_MU = 8; 

const size_t ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_ELE = 16; //byDeepTau2017v2p1VSjet ID working points (deepTau2017v2p1): bitmask 1 = VVVLoose, 2 = VVLoose, 4 = VLoose, 8 = Loose, 16 = Medium, 32 = Tight, 64 = VTight, 128 = VVTight
const size_t ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_MU = 8; //byDeepTau2017v2p1VSjet ID working points (deepTau2017v2p1): bitmask 1 = VVVLoose, 2 = VVLoose, 4 = VLoose, 8 = Loose, 16 = Medium, 32 = Tight, 64 = VTight, 128 = VVTight

const size_t ID_TAU_RECO_DEEPTAU_VSJET=  64; //byDeepTau2017v2p1VSjet ID working points (deepTau2017v2p1): bitmask 1 = VVVLoose, 2 = VVLoose, 4 = VLoose, 8 = Loose, 16 = Medium, 32 = Tight, 64 = VTight, 128 = VVTight
const size_t ID_TAU_RECO_DEEPTAU_VSELE=  4;  //byDeepTau2017v2p1VSe ID working points (deepTau2017v2p1): bitmask 1 = VVVLoose, 2 = VVLoose, 4 = VLoose, 8 = Loose, 16 = Medium, 32 = Tight, 64 = VTight, 128 = VVTight
const size_t ID_TAU_RECO_DEEPTAU_VSMU=   8;  //byDeepTau2017v2p1VSmu ID working points (deepTau2017v2p1): bitmask 1 = VLoose, 2 = Loose, 4 = Medium, 8 = Tight

const size_t ID_TAU_RECO_MVA=            8; //IsolationMVArun2v1DBoldDMwLT ID working point (2017v1): bitmask 1 = VVLoose, 2 = VLoose, 4 = Loose, 8 = Medium, 16 = Tight, 32 = VTight, 64 = VVTight
const size_t ID_TAU_ANTIMU=              1; //Anti-muon discriminator V3: : bitmask 1 = Loose, 2 = Tight
const size_t ID_TAU_ANTIELE=             2; //Anti-electron MVA discriminator V6 (2015): bitmask 1 = VLoose, 2 = Loose, 4 = Medium, 8 = Tight, 16 = VTight
const float PT_CUT_TAU=30;
const float ETA_CUT_TAU=2.3;
const float M_JJ_CUT=   500;
const float MET_CUT=    40;

const string vsJetwp = "VTight";
const string vsElewp = "VLoose";
const string vsMuwp = "Tight";

//const string remote_storage = "https://ttedesch.web.cern.ch/ttedesch/nanoAOD-tools_UL/python/postprocessing/";
const string remote_storage = "https://vbs-pg-support.web.cern.ch/nanoAOD-tools/python/postprocessing/";


TFile *TauID_SF_pt_DeepTau2017v2p1VSjet_UL2017 = TFile::Open(TString(remote_storage) + TString("data/tauSF/TauID_SF_pt_") + TString("DeepTau2017v2p1VSjet") + TString("_") + TString("2017") + TString("ReReco") + TString(".root"));
TString path_down = TString(TString(vsJetwp) + TString("_down"));
TString path_cent = TString(TString(vsJetwp) + TString("_cent"));
TString path_up = TString(TString(vsJetwp) + TString("_up"));
TF1 * TauID_SF_pt_DeepTau2017v2p1VSjet_UL2017_h_down = (TF1*)TauID_SF_pt_DeepTau2017v2p1VSjet_UL2017->Get(path_down);
TF1 * TauID_SF_pt_DeepTau2017v2p1VSjet_UL2017_h_cent =  (TF1*)TauID_SF_pt_DeepTau2017v2p1VSjet_UL2017->Get(path_cent);
TF1 * TauID_SF_pt_DeepTau2017v2p1VSjet_UL2017_h_up =  (TF1*)TauID_SF_pt_DeepTau2017v2p1VSjet_UL2017->Get(path_up);

TFile *TauID_SF_eta_DeepTau2017v2p1VSe_UL2017 = TFile::Open(TString(remote_storage) + TString("data/tauSF/TauID_SF_eta_") + TString("DeepTau2017v2p1VSe") + TString("_") + TString("2017") + TString("ReReco") + TString(".root"));
TString histoname_ele = TString(vsElewp);
TH1F * TauID_SF_eta_DeepTau2017v2p1VSe_UL2017_hist = (TH1F *) TauID_SF_eta_DeepTau2017v2p1VSe_UL2017->Get(histoname_ele);

TFile *TauID_SF_eta_DeepTau2017v2p1VSmu_UL2017 = TFile::Open(TString(remote_storage) + TString("data/tauSF/TauID_SF_eta_") + TString("DeepTau2017v2p1VSmu") + TString("_") + TString("2017") + TString("ReReco") + TString(".root"));
TString histoname_mu = TString(vsMuwp);
TH1F * TauID_SF_eta_DeepTau2017v2p1VSmu_UL2017_hist = (TH1F *) TauID_SF_eta_DeepTau2017v2p1VSmu_UL2017->Get(histoname_mu);

TFile *TauES_dm_DeepTau2017v2p1VSjet_UL2017 = TFile::Open(TString(remote_storage) + TString("data/tauSF/TauES_dm_") + TString("DeepTau2017v2p1VSjet") + TString("_") + TString("2017") + TString("ReReco") +  TString(".root"));        
TH1F * TauES_dm_DeepTau2017v2p1VSjet_UL2017_hist_low = (TH1F *) TauES_dm_DeepTau2017v2p1VSjet_UL2017->Get("tes");

TFile *TauES_dm_DeepTau2017v2p1VSjet_UL2017_ptgt100 = TFile::Open(TString(remote_storage) + TString("data/tauSF/TauES_dm_") + TString("DeepTau2017v2p1VSjet") + TString("_") + TString("2017") + TString("ReReco") + TString("_ptgt100.root"));
TH1F * TauES_dm_DeepTau2017v2p1VSjet_UL2017_ptgt100_hist_high = (TH1F *) TauES_dm_DeepTau2017v2p1VSjet_UL2017_ptgt100->Get("tes");

TFile *TauFES_eta_dm_DeepTau2017v2p1VSe_UL2017 = TFile::Open(TString(remote_storage) + TString("data/tauSF/TauFES_eta-dm_") + TString("DeepTau2017v2p1VSe") + TString("_") + TString("2017")  + TString("ReReco") + TString(".root"));
TGraphAsymmErrors * TauFES_eta_dm_DeepTau2017v2p1VSe_UL2017_graph = (TGraphAsymmErrors *) TauFES_eta_dm_DeepTau2017v2p1VSe_UL2017->Get("fes");

TFile *Btag_eff_UL2017 = TFile::Open(TString(remote_storage) + TString("Btag_eff_") + TString("UL2017") + TString(".root"));
TEfficiency *eff_b = (TEfficiency *) Btag_eff_UL2017->Get("h2_BTaggingEff_b");
TH2F *Btag_eff_UL2017_h_b = (TH2F *) eff_b->CreateHistogram();
TEfficiency *eff_c = (TEfficiency *) Btag_eff_UL2017->Get("h2_BTaggingEff_c");
TH2F *Btag_eff_UL2017_h_c = (TH2F *) eff_c->CreateHistogram();
TEfficiency *eff_udsg = (TEfficiency *) Btag_eff_UL2017->Get("h2_BTaggingEff_udsg");
TH2F *Btag_eff_UL2017_h_udsg = (TH2F *) eff_udsg->CreateHistogram();

//TFile *FR_vsjet2_vsmuT_ZZ = TFile::Open(TString(remote_storage) + TString("FR_vsjet2_vsmuT_ZZ.root"));
TFile *FR_vsjet2_vsmuT_ZZ = TFile::Open(TString(remote_storage) + TString("FR_vsjet2_UL2017.root"));
TH2F *FR_vsjet2_vsmuT_ZZ_histo_ele = (TH2F*)FR_vsjet2_vsmuT_ZZ->Get("FakeRatio_Electron");
TH2F *FR_vsjet2_vsmuT_ZZ_histo_mu = (TH2F*)FR_vsjet2_vsmuT_ZZ->Get("FakeRatio_Muon");
TH2F *FR_vsjet2_vsmuT_ZZ_histo_tau = (TH2F*)FR_vsjet2_vsmuT_ZZ->Get("FakeRatio_Tau");
                                        
//TFile *FR_vsjet4_vsmuT_ZZ = TFile::Open(TString(remote_storage) + TString("FR_vsjet4_vsmuT_ZZ.root"));
TFile *FR_vsjet4_vsmuT_ZZ = TFile::Open(TString(remote_storage) + TString("FR_vsjet4_UL2017.root"));
TH2F *FR_vsjet4_vsmuT_ZZ_histo_ele = (TH2F*)FR_vsjet4_vsmuT_ZZ->Get("FakeRatio_Electron");
TH2F *FR_vsjet4_vsmuT_ZZ_histo_mu = (TH2F*)FR_vsjet4_vsmuT_ZZ->Get("FakeRatio_Muon");
TH2F *FR_vsjet4_vsmuT_ZZ_histo_tau = (TH2F*)FR_vsjet4_vsmuT_ZZ->Get("FakeRatio_Tau");

//TFile *FR_vsjet8_vsmuT_ZZ = TFile::Open(TString(remote_storage) + TString("FR_vsjet4_vsmuT_ZZ.root"));
TFile *FR_vsjet8_vsmuT_ZZ = TFile::Open(TString(remote_storage) + TString("FR_vsjet8_UL2017.root"));
TH2F *FR_vsjet8_vsmuT_ZZ_histo_ele = (TH2F*)FR_vsjet8_vsmuT_ZZ->Get("FakeRatio_Electron");
TH2F *FR_vsjet8_vsmuT_ZZ_histo_mu = (TH2F*)FR_vsjet8_vsmuT_ZZ->Get("FakeRatio_Muon");
TH2F *FR_vsjet8_vsmuT_ZZ_histo_tau = (TH2F*)FR_vsjet8_vsmuT_ZZ->Get("FakeRatio_Tau");

//std::map<string, std::map<string, float> > WP_btagger = {
//    { "CSVv2", {{"L",0.5803},{"M", 0.8838},{"T", 0.9693}}},
//    { "DeepCSV", {{"L", 0.1522},{"M", 0.4941},{"T", 0.8001}}},
//    { "DeepFlv", {{"L", 0.0521},{"M", 0.3033},{"T", 0.7489}}},
//};

float deltaPhi (float phi1, float phi2){
    float dphi = (phi1-phi2);
    while(dphi >  M_PI) dphi -= 2*M_PI;
    while(dphi < -M_PI) dphi += 2*M_PI;
    return dphi;
}

float deltaR(float eta1, float phi1, float eta2, float phi2){
    return hypot(eta1 - eta2, deltaPhi(phi1, phi2)); 
}

RVec<float> getFlattenedMatrixColumn(rvec_f flattened_matrix, int nColumns, int column_index){
    RVec<float> result;
    for (int i = 0; i < flattened_matrix.size()/nColumns; i++) result.emplace_back(flattened_matrix[column_index + i*nColumns]);
    return result;
}


bool LepVetoEle(rvec_i Electron_idx, rvec_f Electron_pt, rvec_f Electron_eta, rvec_b Iso_WPL, rvec_f jetRelIso, rvec_f Muon_pt, rvec_f Muon_eta, rvec_f Iso04_all, rvec_b Muon_looseId )
{
    bool IsEleVetoPassed = true;
    bool IsMuVetoPassed = true;
    for (size_t i = 0; i < Electron_pt.size(); i++) {
        if(i != Electron_idx[0] && Iso_WPL[i] && Electron_pt[i] > PT_CUT_LEP_VETO_ELE && abs(Electron_eta[i]) < ETA_CUT_LEP_VETO_ELE && !(abs(Electron_eta[i])>1.4442 && abs(Electron_eta[i])<1.566) && jetRelIso[i] < REL_ISO_CUT_LEP_VETO_ELE) IsEleVetoPassed = false;
    }
    if(IsEleVetoPassed == true){
        for (size_t i = 0; i < Muon_pt.size(); i++) {
            if(Muon_looseId[i] && Muon_pt[i] > PT_CUT_LEP_VETO_MU && abs(Muon_eta[i]) < ETA_CUT_LEP_VETO_MU && Iso04_all[i] < REL_ISO_CUT_LEP_VETO_MU) IsMuVetoPassed = false;
        }
    }
    return IsEleVetoPassed && IsMuVetoPassed;
}



bool LepVetoMu(rvec_i Muon_idx, rvec_f Electron_pt, rvec_f Electron_eta, rvec_b Electron_mvaFall17V2Iso_WPL, rvec_f Electron_jetRelIso, rvec_f Muon_pt, rvec_f Muon_eta, rvec_f Muon_pfRelIso04_all, rvec_b Muon_looseId)
{
    bool IsEleVetoPassed = true;
    bool IsMuVetoPassed = true;
    for (size_t i = 0; i < Electron_pt.size(); i++) {
        if(Electron_mvaFall17V2Iso_WPL[i] && Electron_pt[i] > PT_CUT_LEP_VETO_ELE && abs(Electron_eta[i]) < ETA_CUT_LEP_VETO_ELE && !(abs(Electron_eta[i])>1.4442 && abs(Electron_eta[i])<1.566) && Electron_jetRelIso[i] < REL_ISO_CUT_LEP_VETO_ELE) IsEleVetoPassed = false;
    }
    if(IsEleVetoPassed == true){
        for (size_t i = 0; i < Muon_pt.size(); i++) {
            if(i != Muon_idx[0] && Muon_looseId[i] && Muon_pt[i] > PT_CUT_LEP_VETO_MU && abs(Muon_eta[i]) < ETA_CUT_LEP_VETO_MU && Muon_pfRelIso04_all[i] < REL_ISO_CUT_LEP_VETO_MU) IsMuVetoPassed = false;
        }
    }
    return IsEleVetoPassed && IsMuVetoPassed;
}

RVec<size_t> GoodJets(rvec_i jetId, rvec_f eta, rvec_f pt, rvec_i puId){
   RVec<int> idx;
   for (size_t i = 0; i < pt.size(); i++) {
      if (jetId[i] >= 2 && abs(eta[i]) < 5. && pt[i] > PT_CUT_JET && (pt[i] > 50. || (pt[i] <= 50. && puId[i] >= 7))) idx.emplace_back(i);
   }
   return idx;
}

bool atleast2GoodJets(rvec_i GoodJets_idx){
    if (GoodJets_idx.size() >= 2) return true;
    else return false;
}


RVec<size_t> SelectVBSJets_invmass(rvec_f pt, rvec_f eta, rvec_f phi, rvec_f mass, rvec_i GoodJets_idx)
{
    
    RVec<size_t> idx;
    // Find first lepton pair with invariant mass closest to Z mass
    auto idx_cmb = Combinations(GoodJets_idx, 2);
    //auto best_mass = -1.;
    float best_mass = -1.;
    size_t best_i1 = 0; size_t best_i2 = 0;
    for (size_t i = 0; i < idx_cmb[0].size(); i++) {
        const auto i1 = idx_cmb[0][i];
        const auto i2 = idx_cmb[1][i];
        //cout<<i1<<i2<<endl;
        if (abs(eta[GoodJets_idx[i1]] - eta[GoodJets_idx[i2]]) >= DELTAETA_JJ_CUT) {
            ROOT::Math::PtEtaPhiMVector p1(pt[GoodJets_idx[i1]], eta[GoodJets_idx[i1]], phi[GoodJets_idx[i1]], mass[GoodJets_idx[i1]]);
            ROOT::Math::PtEtaPhiMVector p2(pt[GoodJets_idx[i2]], eta[GoodJets_idx[i2]], phi[GoodJets_idx[i2]], mass[GoodJets_idx[i2]]);
            //const auto this_mass = (p1 + p2).M();
            float this_mass = (p1 + p2).M();
            //cout<<this_mass<<endl;
            if (this_mass > best_mass) {
                best_mass = this_mass;
                best_i1 = GoodJets_idx[i1];
                best_i2 = GoodJets_idx[i2];
            }
        }
    } 
    //cout<<"best "<<best_i1<<best_i2<<best_mass;
    idx.emplace_back(best_i1);
    idx.emplace_back(best_i2);
    return idx;
}

float GetInvMass(rvec_f pt, rvec_f eta, rvec_f phi, rvec_f mass, rvec_i VBSJets_idx)
{
    ROOT::Math::PtEtaPhiMVector p1(pt[VBSJets_idx[0]], eta[VBSJets_idx[0]], phi[VBSJets_idx[0]], mass[VBSJets_idx[0]]);
    ROOT::Math::PtEtaPhiMVector p2(pt[VBSJets_idx[1]], eta[VBSJets_idx[1]], phi[VBSJets_idx[1]], mass[VBSJets_idx[1]]);
    return (p1 + p2).M();
}

float GetInvMassNoIndex(float pt1, float eta1, float phi1, float mass1, float pt2, float eta2, float phi2, float mass2)
{
    ROOT::Math::PtEtaPhiMVector p1(pt1, eta1, phi1, mass1);
    ROOT::Math::PtEtaPhiMVector p2(pt2, eta2, phi2, mass2);
    return (p1 + p2).M();
}

float GetInvMassNoIndex3(float pt1, float eta1, float phi1, float mass1, float pt2, float eta2, float phi2, float mass2, float pt3, float eta3, float phi3, float mass3)
{
    ROOT::Math::PtEtaPhiMVector p1(pt1, eta1, phi1, mass1);
    ROOT::Math::PtEtaPhiMVector p2(pt2, eta2, phi2, mass2);
    ROOT::Math::PtEtaPhiMVector p3(pt3, eta3, phi3, mass3);
    return (p1 + p2 + p3).M();
}

float GetInvMassNoIndex4(float pt1, float eta1, float phi1, float mass1, float pt2, float eta2, float phi2, float mass2, float pt3, float eta3, float phi3, float mass3, float pt4, float eta4, float phi4, float mass4)
{
    ROOT::Math::PtEtaPhiMVector p1(pt1, eta1, phi1, mass1);
    ROOT::Math::PtEtaPhiMVector p2(pt2, eta2, phi2, mass2);
    ROOT::Math::PtEtaPhiMVector p3(pt3, eta3, phi3, mass3);
    ROOT::Math::PtEtaPhiMVector p4(pt4, eta4, phi3, mass4);
    return (p1 + p2 + p3 + p4).M();
}

//float ComputeTransverseMass(float met_et, float met_phi, float lep_pt, float lep_eta, float lep_phi, float lep_e)
//{
//    ROOT::Math::PtEtaPhiEVector met(met_et, 0, met_phi, met_et);
//    ROOT::Math::PtEtaPhiEVector lep(lep_pt, lep_eta, lep_phi, lep_e);
//    return (met + lep).Mt() / 1000.0;
//}


float GetLeading(rvec_f Jet_pt, rvec_i VBSJet_idx){
    return Jet_pt[VBSJet_idx[0]];
}

float GetSubLeading(rvec_f Jet_pt, rvec_i VBSJet_idx){
    return Jet_pt[VBSJet_idx[1]];
}

float GetLepton(rvec_f Electron_pt, rvec_i Electron_idx, rvec_f Muon_pt, rvec_i Muon_idx, int GoodLeptonFamily){
    if (GoodLeptonFamily == 0) return Electron_pt[Electron_idx[0]];
    else return Muon_pt[Muon_idx[0]];
}

float GetLeptonSF(rvec_f Electron_pt, rvec_i Electron_idx, rvec_f Muon_pt, rvec_i Muon_idx, int GoodLeptonFamily, bool IsMC){
    if (IsMC == false) return 1;
    else {
        if (GoodLeptonFamily == 0) return Electron_pt[Electron_idx[0]];
        else return Muon_pt[Muon_idx[0]];
    }
}

//float GetLeptonTightFlag(rvec_i Electron_idx, rvec_i Muon_idx, int GoodLeptonFamily){
int GetLeptonTightFlag(rvec_i Electron_idx, rvec_i Muon_idx, int GoodLeptonFamily){
    if (GoodLeptonFamily == 0) return Electron_idx[1];
    else return Muon_idx[1];
}

float GetTau(rvec_f pt, rvec_i idx){
    return pt[idx[0]];
}

RVec<int> SelectElectron(rvec_f lepton_pt, rvec_f lepton_eta, rvec_f lepton_phi, rvec_f lepton_jetRelIso, rvec_b lepton_mvaFall17V2Iso_WPL, rvec_f lepton_mvaFall17V2Iso_WP90, rvec_f jet_eta, rvec_f jet_phi, rvec_i VBSJets_idx, string Year){
    //setting jet-related quantities if isolation from them is needed
    //const auto jet1_idx = VBSJets_idx[0];
    //const auto jet2_idx = VBSJets_idx[1];

    //const auto jet1eta = jet_eta[jet1_idx];
    //const auto jet2eta = jet_eta[jet2_idx];
    //const auto jet1phi = jet_phi[jet1_idx];
    //const auto jet2phi = jet_phi[jet2_idx];
    
    float PT_CUT_ELE;
    
    if(Year == "UL2016" || Year == "UL2016APV") PT_CUT_ELE = PT_CUT_ELE_UL2016;
    else if (Year == "UL2017") PT_CUT_ELE = PT_CUT_ELE_UL2017;
    else PT_CUT_ELE = PT_CUT_ELE_UL2018;
    
    float jet1_idx = VBSJets_idx[0];
    float jet2_idx = VBSJets_idx[1];
    //const auto jet1eta = lepton_pt[0]
    float jet1eta = jet_eta[jet1_idx];
    float jet2eta = jet_eta[jet2_idx];
    float jet1phi = jet_phi[jet1_idx];
    float jet2phi = jet_phi[jet2_idx];
    
    const float isocone = DR_OVERLAP_CONE_OTHER;

    RVec<size_t> Tleptons_idx;
    RVec<size_t> LnTleptons_idx;
    
    bool IsLooseID, IsTightID, IsTightIso, IsLooseIso, IsInEtaRegion, IsInPtRegion;
    
    for (size_t i = 0; i < lepton_pt.size(); i++) {
        //setting loose and tight, eta, and pt criteria for leptons depending on lepton flavour
        IsLooseID = lepton_mvaFall17V2Iso_WPL[i];
        IsTightID = lepton_mvaFall17V2Iso_WP90[i];
        IsTightIso = lepton_jetRelIso[i]<ISO_CUT_ELE && lepton_jetRelIso[i]>=0.;
        IsLooseIso = lepton_jetRelIso[i]<1. && lepton_jetRelIso[i]>=0.;
        IsInEtaRegion = abs(lepton_eta[i])<ETA_CUT_ELE && !(abs(lepton_eta[i])>1.4442 && abs(lepton_eta[i])<1.566);
        IsInPtRegion = lepton_pt[i] > PT_CUT_ELE;

        //find tight and loose-not-tight leptons filtering with jet-lep isolation criteria
        if (IsInEtaRegion && IsInPtRegion){
            if(IsLooseID && IsLooseIso){
                if(deltaR(lepton_eta[i], lepton_phi[i], jet1eta,  jet1phi) > isocone && deltaR(lepton_eta[i], lepton_phi[i], jet2eta,  jet2phi) > isocone){
                    if(IsTightID && IsTightIso) Tleptons_idx.emplace_back(i);
                    else LnTleptons_idx.emplace_back(i);
                }
            }
        }
    }
 
    RVec<int> idx(2);
    //select leading tight/loose-not-tight lepton
    if (Tleptons_idx.size() > 0){
        idx[0] = Tleptons_idx[0];
        idx[1] = 1;
    }
    else if (LnTleptons_idx.size() > 0){
        idx[0] = LnTleptons_idx[0];
        idx[1] = 0;
    }
    else{
        idx[0] = -1;
        idx[1] = -1;  
    }
    return idx;
}


RVec<int> SelectMuon(rvec_f lepton_pt, rvec_f lepton_eta, rvec_f lepton_phi, rvec_b lepton_tightId, rvec_b lepton_looseId, rvec_f Iso04_all, rvec_f jet_eta, rvec_f jet_phi, rvec_i VBSJets_idx){
    //setting jet-related quantities if isolation from them is needed
    //const auto jet1_idx = VBSJets_idx[0];
    //const auto jet2_idx = VBSJets_idx[1];
    //const auto jet1eta = lepton_pt[0]
    //const auto jet1eta = jet_eta[jet1_idx];
    //const auto jet2eta = jet_eta[jet2_idx];
    //const auto jet1phi = jet_phi[jet1_idx];
    //const auto jet2phi = jet_phi[jet2_idx];
    
    float jet1_idx = VBSJets_idx[0];
    float jet2_idx = VBSJets_idx[1];
    //const auto jet1eta = lepton_pt[0]
    float jet1eta = jet_eta[jet1_idx];
    float jet2eta = jet_eta[jet2_idx];
    float jet1phi = jet_phi[jet1_idx];
    float jet2phi = jet_phi[jet2_idx];
    
    const float isocone = DR_OVERLAP_CONE_OTHER;

    RVec<size_t> Tleptons_idx;
    RVec<size_t> LnTleptons_idx;
    
    bool IsLooseID, IsTightID, IsTightIso, IsLooseIso, IsInEtaRegion, IsInPtRegion;
    
    for (size_t i = 0; i < lepton_pt.size(); i++) {
        //setting loose and tight, eta, and pt criteria for leptons depending on lepton flavour
        IsTightID = lepton_tightId[i];
        IsLooseID = lepton_looseId[i];
        IsTightIso = Iso04_all[i]<ISO_CUT_MU && Iso04_all[i]>=0.;
        IsLooseIso = Iso04_all[i]<1. && Iso04_all[i]>=0.;
        IsInEtaRegion = abs(lepton_eta[i]) < ETA_CUT_MU;
        IsInPtRegion = lepton_pt[i] > PT_CUT_MU;
        

        //find tight and loose-not-tight leptons filtering with jet-lep isolation criteria
        if (IsInEtaRegion && IsInPtRegion){
            if(IsLooseID && IsLooseIso){
                if(deltaR(lepton_eta[i], lepton_phi[i], jet1eta,  jet1phi) > isocone && deltaR(lepton_eta[i], lepton_phi[i], jet2eta,  jet2phi) > isocone){
                    if(IsTightID && IsTightIso) Tleptons_idx.emplace_back(i);
                    else LnTleptons_idx.emplace_back(i);
                }
            }
        }
    }
    RVec<int> idx(2);
    //select leading tight/loose-not-tight lepton
    if (Tleptons_idx.size() > 0){
        idx[0] = Tleptons_idx[0];
        idx[1] = 1;
    }
    else if (LnTleptons_idx.size() > 0){
        idx[0] = LnTleptons_idx[0];
        idx[1] = 0;
    }
    else{
        idx[0] = -1;
        idx[1] = -1;  
    }
    
    return idx;
}


int DetermineGoodLepton(bool HLT_IsoMu27, bool HLT_Mu50, bool HLT_Ele35_WPTight_Gsf, bool HLT_Ele32_WPTight_Gsf_L1DoubleEG, bool HLT_Photon200, bool HLT_PFHT250, bool HLT_PFHT350, rvec_i Electron_idx, rvec_f Electron_pt, rvec_f Electron_eta, rvec_b Electron_mvaFall17V2Iso_WPL, rvec_f Electron_jetRelIso, rvec_i Muon_idx, rvec_f Muon_pt, rvec_f Muon_eta, rvec_f Muon_pfRelIso04_all, rvec_b Muon_looseId){
    bool passMu = false;
    bool passEle = false;
    bool passHT = false;
    int GoodLeptonFamily;
    
    if(HLT_IsoMu27 || HLT_Mu50) passMu = true;
    if(HLT_Ele35_WPTight_Gsf || HLT_Ele32_WPTight_Gsf_L1DoubleEG || HLT_Photon200) passEle = true;
    if(HLT_PFHT250 || HLT_PFHT350) passHT = true;
    
    bool ele_lepton_veto = false;
    bool mu_lepton_veto = false;
    
    if(Electron_idx[1] != -1) ele_lepton_veto = LepVetoEle(Electron_idx, Electron_pt, Electron_eta, Electron_mvaFall17V2Iso_WPL, Electron_jetRelIso, Muon_pt, Muon_eta, Muon_pfRelIso04_all, Muon_looseId);
    if(Muon_idx[1] != -1) mu_lepton_veto = LepVetoMu(Muon_idx, Electron_pt, Electron_eta, Electron_mvaFall17V2Iso_WPL, Electron_jetRelIso, Muon_pt, Muon_eta, Muon_pfRelIso04_all, Muon_looseId);

    bool SingleEle=false;
    bool SingleMu=false;
    if(passEle && !passMu){
        if(Electron_idx[1] != -1 && ele_lepton_veto){
            GoodLeptonFamily = 0;
            //lepton_TightRegion[0] = copy.deepcopy(ele_TightRegion)
            SingleEle = true;
            SingleMu = false;
        }
        else GoodLeptonFamily = -1;
    }
    else if(!passEle && passMu){
        if(Muon_idx[1] != -1 && mu_lepton_veto){
            GoodLeptonFamily = 1;
            //lepton_TightRegion[0] = copy.deepcopy(mu_TightRegion)
            SingleEle = false;
            SingleMu = true;
        }
        else GoodLeptonFamily = -1;
    }
    else if(passEle && passMu){ 
        if(Muon_idx[1] == -1 && Electron_idx[1] != -1 && ele_lepton_veto){
            GoodLeptonFamily = 0;
            //lepton_TightRegion[0] = copy.deepcopy(ele_TightRegion)
            SingleEle = true;
            SingleMu = false;
        }
        else if(Electron_idx[1] == -1 && Muon_idx[1] != -1 && mu_lepton_veto){
            GoodLeptonFamily = 1;
            //lepton_TightRegion[0] = copy.deepcopy(mu_TightRegion)
            SingleEle = false;
            SingleMu = true;
        }
                
        else if(Muon_idx[1] != -1 && Electron_idx[1] != -1){
            if(ele_lepton_veto && !mu_lepton_veto){
                GoodLeptonFamily = 0;
                //lepton_TightRegion[0] = copy.deepcopy(ele_TightRegion)
                SingleEle = false;
                SingleMu = true;
            }
            else if(!ele_lepton_veto && mu_lepton_veto){           
                GoodLeptonFamily = 1;
                //lepton_TightRegion[0] = copy.deepcopy(mu_TightRegion)
                SingleEle = false;
                SingleMu = true;
            }

            else if(ele_lepton_veto && mu_lepton_veto){
                if(Electron_pt[Electron_idx[0]] > Muon_pt[Muon_idx[0]]){
                    GoodLeptonFamily = 0;
                    //lepton_TightRegion[0] = copy.deepcopy(ele_TightRegion)
                    SingleEle = true;
                    SingleMu = false;
                }
                else{
                    GoodLeptonFamily = 1;
                    //lepton_TightRegion[0] = copy.deepcopy(mu_TightRegion)
                    SingleEle = false;
                    SingleMu = true;
                }
            }
            else GoodLeptonFamily = -1;
        }
        else GoodLeptonFamily = -1;
    }
    if(!(SingleEle || SingleMu)) GoodLeptonFamily = -1;
    
    return GoodLeptonFamily;
}


RVec<int> SelectAndVetoTaus(rvec_f Tau_pt, rvec_f Tau_eta, rvec_f Tau_phi, RVec<UChar_t> &Tau_idDeepTau2017v2p1VSjet, RVec<UChar_t> &Tau_idDeepTau2017v2p1VSe, RVec<UChar_t> &Tau_idDeepTau2017v2p1VSmu, int GoodLeptonFamily, rvec_i Electron_idx, rvec_f Electron_eta, rvec_f Electron_phi, rvec_i Muon_idx, rvec_f Muon_eta, rvec_f Muon_phi, rvec_f Jet_eta, rvec_f Jet_phi, rvec_i VBSJet_idx, string Year)
{
    size_t ID_TAU_RECO_DEEPTAU_VSJET_VETO_ELE,  ID_TAU_RECO_DEEPTAU_VSJET_VETO_MU;
    
    if(Year == "UL2016" || Year == "UL2016APV") ID_TAU_RECO_DEEPTAU_VSJET_VETO_ELE = ID_TAU_RECO_DEEPTAU_VSJET_VETO_ELE_UL2016;
    else if (Year == "UL2017") ID_TAU_RECO_DEEPTAU_VSJET_VETO_ELE = ID_TAU_RECO_DEEPTAU_VSJET_VETO_ELE_UL2017;
    else ID_TAU_RECO_DEEPTAU_VSJET_VETO_ELE = ID_TAU_RECO_DEEPTAU_VSJET_VETO_ELE_UL2018;
    
    if(Year == "UL2016" || Year == "UL2016APV") ID_TAU_RECO_DEEPTAU_VSJET_VETO_MU = ID_TAU_RECO_DEEPTAU_VSJET_VETO_MU_UL2016;
    else if (Year == "UL2017") ID_TAU_RECO_DEEPTAU_VSJET_VETO_MU = ID_TAU_RECO_DEEPTAU_VSJET_VETO_MU_UL2017;
    else ID_TAU_RECO_DEEPTAU_VSJET_VETO_MU = ID_TAU_RECO_DEEPTAU_VSJET_VETO_MU_UL2018;
    
    //setting jet-related quantities if isolation from them is needed
    
    float jet1eta = Jet_eta[VBSJet_idx[0]];
    float jet2eta = Jet_eta[VBSJet_idx[1]];
    float jet1phi = Jet_phi[VBSJet_idx[0]];
    float jet2phi = Jet_phi[VBSJet_idx[1]];
    float isocone = DR_OVERLAP_CONE_OTHER;

    size_t nTau=0;
    RVec<int> idx(2);
    
    if (Tau_eta.size()==0){
        idx[0] = -1;
        idx[1] = -1;
        return idx;
    } 
    float cutloose_vsjet;
    for (size_t i = 0; i < Tau_eta.size(); i++) {
        
        if (GoodLeptonFamily == 0) {
            cutloose_vsjet = ID_TAU_RECO_DEEPTAU_VSJET_VETO_ELE;
            if ((Tau_idDeepTau2017v2p1VSjet[i]>=cutloose_vsjet && Tau_idDeepTau2017v2p1VSe[i]>=ID_TAU_RECO_DEEPTAU_VSELE && Tau_idDeepTau2017v2p1VSmu[i]>=ID_TAU_RECO_DEEPTAU_VSMU) && deltaR(Tau_eta[i], Tau_phi[i], Electron_eta[Electron_idx[0]], Electron_phi[Electron_idx[0]])>DR_OVERLAP_CONE_TAU && deltaR(Tau_eta[i], Tau_phi[i], jet1eta, jet1phi)>isocone && deltaR(Tau_eta[i], Tau_phi[i], jet2eta, jet2phi)>isocone && Tau_pt[i]>=PT_CUT_TAU && abs(Tau_eta[i])<=ETA_CUT_TAU){
                nTau++;
                if(Tau_idDeepTau2017v2p1VSjet[i]>=ID_TAU_RECO_DEEPTAU_VSJET){
                    idx[0] = i;
                    //idx[1] = 0;
                    idx[1] = 1;
                } 
                else{
                    idx[0] = i;
                    //idx[1] = 1;
                    idx[1] = 0;
                }
            }   
        }

        else if (GoodLeptonFamily == 1) {
            cutloose_vsjet = ID_TAU_RECO_DEEPTAU_VSJET_VETO_MU;
            if ((Tau_idDeepTau2017v2p1VSjet[i]>=cutloose_vsjet && Tau_idDeepTau2017v2p1VSe[i]>=ID_TAU_RECO_DEEPTAU_VSELE && Tau_idDeepTau2017v2p1VSmu[i]>=ID_TAU_RECO_DEEPTAU_VSMU) && deltaR(Tau_eta[i], Tau_phi[i], Muon_eta[Muon_idx[0]], Muon_phi[Muon_idx[0]])>DR_OVERLAP_CONE_TAU && deltaR(Tau_eta[i], Tau_phi[i], jet1eta, jet1phi)>isocone && deltaR(Tau_eta[i], Tau_phi[i], jet2eta, jet2phi)>isocone && Tau_pt[i]>=PT_CUT_TAU && abs(Tau_eta[i])<=ETA_CUT_TAU){
                nTau++;
                if(Tau_idDeepTau2017v2p1VSjet[i]>=ID_TAU_RECO_DEEPTAU_VSJET){
                    idx[0] = i;
                    idx[1] = 1;
                } 
                else{
                    idx[0] = i;
                    idx[1] = 0;
                }
            }
        }
    }
    if(nTau!=1) idx[1] = -1;                                                                                               
    return idx;
}

bool SameCharge(int GoodLeptonFamily, rvec_i Electron_idx, rvec_i Electron_charge, rvec_i Muon_idx, rvec_i Muon_charge, rvec_i Tau_idx, rvec_i Tau_charge){
    if(GoodLeptonFamily == 0){
        if(Electron_charge[Electron_idx[0]] == Tau_charge[Tau_idx[0]]) return true;
        else return false;
    }
    else if(GoodLeptonFamily == 1){
        if(Muon_charge[Muon_idx[0]] == Tau_charge[Tau_idx[0]]) return true;
        else return false;
    }
    return false;
}

bool BVeto(rvec_f Jet_pt, rvec_f Jet_eta, rvec_f Jet_btagDeepFlavB, rvec_i GoodJet_idx)
{
    bool veto = false;
    for (size_t i = 0; i < GoodJet_idx.size(); i++) {
        //if (Jet_btagDeepFlavB[i]>=WP_btagger[BTAG_ALGO][BTAG_WP])*(Jet_pt[i]>BTAG_PT_CUT)*(abs(Jet_eta[i])<BTAG_ETA_CUT) return true;
        //if ((Jet_btagDeepFlavB[i]>=BTAG_WP_VALUE) && (Jet_pt[i]>BTAG_PT_CUT) && (abs(Jet_eta[i])<BTAG_ETA_CUT)) return false;
        if ((Jet_btagDeepFlavB[GoodJet_idx[i]]>=BTAG_WP_VALUE) && (Jet_pt[GoodJet_idx[i]]>BTAG_PT_CUT) && (abs(Jet_eta[GoodJet_idx[i]])<BTAG_ETA_CUT)) return false;
    }
    return true;
}

bool BVeto_loose(rvec_f Jet_pt, rvec_f Jet_eta, rvec_f Jet_btagDeepFlavB, rvec_i GoodJet_idx)
{
    bool veto = false;
    for (size_t i = 0; i < GoodJet_idx.size(); i++) {
        //if (Jet_btagDeepFlavB[i]>=WP_btagger[BTAG_ALGO][BTAG_WP])*(Jet_pt[i]>BTAG_PT_CUT)*(abs(Jet_eta[i])<BTAG_ETA_CUT) return true;
        //if ((Jet_btagDeepFlavB[i]>=BTAG_WP_VALUE) && (Jet_pt[i]>BTAG_PT_CUT) && (abs(Jet_eta[i])<BTAG_ETA_CUT)) return false;
        if ((Jet_btagDeepFlavB[GoodJet_idx[i]]>=BTAG_WP_VALUE_LOOSE) && (Jet_pt[GoodJet_idx[i]]>BTAG_PT_CUT) && (abs(Jet_eta[GoodJet_idx[i]])<BTAG_ETA_CUT)) return false;
    }
    return true;
}


float GetLog2(float x){
    if(x > 0.) return log2(x);
    else return 0;
}

float deltaTheta(float pt1, float eta1, float phi1, float mass1, float pt2, float eta2, float phi2, float mass2){
    ROOT::Math::PtEtaPhiMVector p1(pt1, eta1, phi1, mass1);
    ROOT::Math::PtEtaPhiMVector p2(pt2, eta2, phi2, mass2);
    return cos((p1 - p2).Theta());
}

float Zeppenfeld(float lep_eta, float ljet_eta, float sljet_eta){
    float zepp_lepjj = lep_eta - 0.5*(ljet_eta+sljet_eta);
    return zepp_lepjj;
}

float M1T(float Lepton_pt, float Lepton_eta, float Lepton_phi, float Lepton_mass, float SelectedTau_pt, float SelectedTau_eta, float SelectedTau_phi, float SelectedTau_mass, float MET_pt, float MET_phi){
    ROOT::Math::PtEtaPhiMVector lep_p4(Lepton_pt, Lepton_eta, Lepton_phi, Lepton_mass);
    ROOT::Math::PtEtaPhiMVector tau_p4(SelectedTau_pt, SelectedTau_eta, SelectedTau_phi, SelectedTau_mass); //*taucorr
    auto leptau_p4 = lep_p4 + tau_p4;
    auto leptau_pt2 = leptau_p4.Perp2();
    auto leptau_px = leptau_p4.Px();
    auto leptau_py = leptau_p4.Py();
    auto leptau_mass2 = leptau_p4.M2();
    auto leptau_et = sqrt(leptau_mass2 + leptau_pt2);
    auto MET_px = MET_pt*cos(MET_phi);
    auto MET_py = MET_pt*sin(MET_phi);
    auto sys_et2 = pow(leptau_et + MET_pt,2.);
    auto sys_pt2 = pow(leptau_px + MET_px, 2.) + pow(leptau_py + MET_py, 2.);
    auto M1T2 = sys_et2 - sys_pt2;
    auto sign_M1T2 = M1T2/abs(M1T2);

    return sign_M1T2*sqrt(sign_M1T2*M1T2);
}

float Mo1(float Lepton_pt, float Lepton_eta, float Lepton_phi, float Lepton_mass, float SelectedTau_pt, float SelectedTau_eta, float SelectedTau_phi, float SelectedTau_mass, float MET_pt, float MET_phi){
    ROOT::Math::PtEtaPhiMVector lep_p4(Lepton_pt, Lepton_eta, Lepton_phi, Lepton_mass);
    ROOT::Math::PtEtaPhiMVector tau_p4(SelectedTau_pt, SelectedTau_eta, SelectedTau_phi, SelectedTau_mass); //*taucorr
    auto leptau_p4 = lep_p4 + tau_p4;
    auto lep_pt = Lepton_pt;
    auto tau_pt = SelectedTau_pt; //*taucorr
    auto leptau_px = leptau_p4.Px();
    auto leptau_py = leptau_p4.Py();
    auto MET_px = MET_pt*cos(MET_phi);
    auto MET_py = MET_pt*sin(MET_phi);
    auto sys_eo2 = pow(lep_pt + tau_pt + MET_pt, 2.);
    auto sys_pt2 = pow(leptau_px + MET_px, 2.) + pow(leptau_py + MET_py, 2.);
    auto Mo12 = sys_eo2 - sys_pt2;
    auto sign_Mo12 = Mo12/abs(Mo12);

    return sign_Mo12*sqrt(sign_Mo12*Mo12);
}


float SFFakeRatio_lep_calc_vsjet2(float pT, float eta, int pdgId){
    
    //TFile inFile("FR_vsjet2_vsmuT_ZZ.root"); 
    //TFile *inFile = FR_vsjet2_vsmuT_ZZ;
    TH2F *histo;
    //if (abs(pdgId) == 11) histo = (TH2F*)inFile.Get("hFRDataeledif");
    //else if (abs(pdgId) == 13) histo = (TH2F*)inFile.Get("hFRDatamudif");
    //if (abs(pdgId) == 11) histo = (TH2F*)inFile->Get("hFRDataeledif");
    //else if (abs(pdgId) == 13) histo = (TH2F*)inFile->Get("hFRDatamudif");
    if (abs(pdgId) == 11) histo = FR_vsjet2_vsmuT_ZZ_histo_ele;
    else if (abs(pdgId) == 13) histo = FR_vsjet2_vsmuT_ZZ_histo_mu;

    auto binx = histo->GetXaxis()->FindBin(pT);
    auto biny = histo->GetYaxis()->FindBin(eta);
    auto nxbins = histo->GetXaxis()->GetNbins();
    auto nybins = histo->GetYaxis()->GetNbins();
        
    if(binx > nxbins) binx = nxbins;
    else if (binx <= 0) binx = 1;
    if (biny > nybins) biny = nybins;
    else if (biny <= 0) biny = 1;

    auto FR = histo->GetBinContent(binx, biny);

    return FR/(1-FR);
}
        
float SFFakeRatio_lep_calc_vsjet4(float pT, float eta, int pdgId){
    
    //TFile inFile("FR_vsjet4_vsmuT_ZZ.root"); 
    //TFile *inFile = FR_vsjet4_vsmuT_ZZ;
    TH2F *histo;
    //if (abs(pdgId) == 11) histo = (TH2F*)inFile.Get("hFRDataeledif");
    //else if (abs(pdgId) == 13) histo = (TH2F*)inFile.Get("hFRDatamudif");
    //if (abs(pdgId) == 11) histo = (TH2F*)inFile->Get("hFRDataeledif");
    //else if (abs(pdgId) == 13) histo = (TH2F*)inFile->Get("hFRDatamudif");
    if (abs(pdgId) == 11) histo = FR_vsjet4_vsmuT_ZZ_histo_ele;
    else if (abs(pdgId) == 13) histo = FR_vsjet4_vsmuT_ZZ_histo_mu;

    auto binx = histo->GetXaxis()->FindBin(pT);
    auto biny = histo->GetYaxis()->FindBin(eta);
    auto nxbins = histo->GetXaxis()->GetNbins();
    auto nybins = histo->GetYaxis()->GetNbins();
        
    if(binx > nxbins) binx = nxbins;
    else if (binx <= 0) binx = 1;
    if (biny > nybins) biny = nybins;
    else if (biny <= 0) biny = 1;

    auto FR = histo->GetBinContent(binx, biny);

    return FR/(1-FR);
}

        
float SFFakeRatio_tau_calc_vsjet2(float pT, float eta){
    //TFile inFile("FR_vsjet2_vsmuT_ZZ.root"); 
    //TFile *inFile = FR_vsjet2_vsmuT_ZZ;

    //TH2F *histo = (TH2F*)inFile.Get("hFRDatataudif");
    //TH2F *histo = (TH2F*)inFile->Get("hFRDatataudif");
    
    TH2F *histo = FR_vsjet2_vsmuT_ZZ_histo_tau;

    auto binx = histo->GetXaxis()->FindBin(pT);
    //auto biny = histo->GetYaxis()->FindBin(eta);
    auto biny = histo->GetYaxis()->FindBin(abs(eta));
    auto nxbins = histo->GetXaxis()->GetNbins();
    auto nybins = histo->GetYaxis()->GetNbins();
        
    if(binx > nxbins) binx = nxbins;
    else if (binx <= 0) binx = 1;
    if (biny > nybins) biny = nybins;
    else if (biny <= 0) biny = 1;

    auto FR = histo->GetBinContent(binx, biny);

    return FR/(1-FR);
}
        
float SFFakeRatio_tau_calc_vsjet4(float pT, float eta){

    //TFile inFile("FR_vsjet4_vsmuT_ZZ.root"); 
    //TFile *inFile = FR_vsjet4_vsmuT_ZZ;

    //TH2F *histo = (TH2F*)inFile.Get("hFRDatataudif");
    //TH2F *histo = (TH2F*)inFile->Get("hFRDatataudif");
    
    TH2F *histo = FR_vsjet4_vsmuT_ZZ_histo_tau;

    auto binx = histo->GetXaxis()->FindBin(pT);
    //auto biny = histo->GetYaxis()->FindBin(eta);
    auto biny = histo->GetYaxis()->FindBin(abs(eta));
    auto nxbins = histo->GetXaxis()->GetNbins();
    auto nybins = histo->GetYaxis()->GetNbins();
        
    if(binx > nxbins) binx = nxbins;
    else if (binx <= 0) binx = 1;
    if (biny > nybins) biny = nybins;
    else if (biny <= 0) biny = 1;

    auto FR = histo->GetBinContent(binx, biny);

    return FR/(1-FR);
}

float GetEventSFFake(float lepton_SFFake, float tau_SFFake, int lepton_LnTRegion, int tau_LnTRegion){
    if(lepton_LnTRegion==1 && tau_LnTRegion==0) return lepton_SFFake;
    else if (lepton_LnTRegion==0 && tau_LnTRegion==1) return tau_SFFake;
    else if (lepton_LnTRegion==1 && tau_LnTRegion==1) return lepton_SFFake*tau_SFFake;
    else if (lepton_LnTRegion==0 && tau_LnTRegion==0) return 0.;
}

/*
RVec<int> SelectVBSQGenJet(rvec_i GenPart_pdgId, rvec_i GenPart_genPartIdxMother, rvec_f GenPart_pt, rvec_f GenPart_eta, rvec_i GenJet_partonFlavour, rvec_f GenJet_pt, rvec_f GenJet_eta){

    RVec<int> GenPart_idx;
    for (int i = 0; i < GenPart_pdgId.size(); i++) {
        if(GenPart_genPartIdxMother[i]==0 and abs(GenPart_pdgId[i])>0 and abs(GenPart_pdgId[i])<10) GenPart_idx.emplace_back(i);
    }
    
    int GenPart_idx1 = GenPart_idx[0];
    int GenPart_idx2 = GenPart_idx[1];
    
    int qflav1 = GenPart_pdgId[GenPart_idx1];
    int qflav2 = GenPart_pdgId[GenPart_idx2];
    
    RVec<int> GenJet_idx;

    for (int i = 0; i < GenJet_partonFlavour.size(); i++) {
        if(abs(GenJet_partonFlavour[i])>0 && abs(GenJet_partonFlavour[i])<10 && (GenJet_partonFlavour[i]==qflav1 || GenJet_partonFlavour[i]==qflav2)) GenJet_idx.emplace_back(i);
    }
    
    float discrim1 = 1000000.;
    float discrim2 = 1000000.;
    int idx_genjet1 = -1;
    int idx_genjet2 = -1;
    float tmpdiscr1;
    float tmpdiscr2;
           
    for (int i = 0; i < GenJet_idx.size(); i++) {
        tmpdiscr1 = abs(GenJet_eta[GenJet_idx[i]] -  GenPart_eta[GenPart_idx1]) + abs(GenJet_pt[GenJet_idx[i]] -  GenPart_pt[GenPart_idx1]);
        tmpdiscr2 = abs(GenJet_eta[GenJet_idx[i]] -  GenPart_eta[GenPart_idx1]) + abs(GenJet_pt[GenJet_idx[i]] -  GenPart_pt[GenPart_idx1]);
        if(tmpdiscr1 < discrim1){
            discrim1 = tmpdiscr1;
            idx_genjet1 = GenJet_idx[i];
        }
        if(tmpdiscr2 < discrim2){
            discrim2 = tmpdiscr2;
            idx_genjet2 = GenJet_idx[i];
        }
    }
        
    RVec<int> finalgenjets_idx(2);
           
    if(GenJet_pt[idx_genjet1] > GenJet_pt[idx_genjet2]){ 
           finalgenjets_idx[0] = idx_genjet1;
           finalgenjets_idx[1] = idx_genjet2;
    }
    else{ 
           finalgenjets_idx[0] = idx_genjet2;
           finalgenjets_idx[1] = idx_genjet1;
    }
    return finalgenjets_idx;
}
*/

RVec<int> SelectVBSQGenJet(rvec_i GenPart_pdgId, rvec_i GenPart_genPartIdxMother, rvec_f GenPart_pt, rvec_f GenPart_eta, rvec_i GenJet_partonFlavour, rvec_f GenJet_pt, rvec_f GenJet_eta){

    RVec<int> GenPart_idx;
    
    for (int i = 0; i < GenPart_pdgId.size(); i++) {
        if(GenPart_genPartIdxMother[i]==0 && abs(GenPart_pdgId[i])>0 && abs(GenPart_pdgId[i])<10) GenPart_idx.emplace_back(i);
    }
    
    RVec<int> dummy_idx;
    dummy_idx.emplace_back(-9999);
    dummy_idx.emplace_back(-9999);
    
    if (GenPart_idx.size() < 1) return dummy_idx;
    
    int GenPart_idx1 = GenPart_idx[0];
    int GenPart_idx2 = GenPart_idx[1];
    
    if(GenPart_pt[GenPart_idx1] == GenPart_pt[GenPart_idx2] && GenPart_eta[GenPart_idx1] == GenPart_eta[GenPart_idx2]) return dummy_idx;
    
    int qflav1 = GenPart_pdgId[GenPart_idx1];
    int qflav2 = GenPart_pdgId[GenPart_idx2];
    
    RVec<int> LightGenJet_idx;

    for (int i = 0; i < GenJet_partonFlavour.size(); i++) {
        if(abs(GenJet_partonFlavour[i])>0 && abs(GenJet_partonFlavour[i])<10 && (GenJet_partonFlavour[i]==qflav1 || GenJet_partonFlavour[i]==qflav2)) LightGenJet_idx.emplace_back(i);
    }
    
    if(LightGenJet_idx.size() < 2){
        LightGenJet_idx.clear();
        for (int i = 0; i < GenJet_partonFlavour.size(); i++) {
            if(abs(GenJet_partonFlavour[i])>0 && abs(GenJet_partonFlavour[i])<10) LightGenJet_idx.emplace_back(i);
        }
    }
    
    float discrim1 = 1000000.;
    float discrim2 = 1000000.;
    int idx_genjet1 = -1;
    int idx_genjet2 = -1;
    float tmpdiscr1;
    float tmpdiscr2;
           
    for (int i = 0; i < LightGenJet_idx.size(); i++) {
        tmpdiscr1 = abs(GenJet_eta[LightGenJet_idx[i]] -  GenPart_eta[GenPart_idx1]) + abs(GenJet_pt[LightGenJet_idx[i]] -  GenPart_pt[GenPart_idx1]);
        tmpdiscr2 = abs(GenJet_eta[LightGenJet_idx[i]] -  GenPart_eta[GenPart_idx2]) + abs(GenJet_pt[LightGenJet_idx[i]] -  GenPart_pt[GenPart_idx2]);
        if(tmpdiscr1 < discrim1){
            discrim1 = tmpdiscr1;
            idx_genjet1 = LightGenJet_idx[i];
        }
        if(tmpdiscr2 < discrim2){
            discrim2 = tmpdiscr2;
            idx_genjet2 = LightGenJet_idx[i];
        }
    }
    
    //cout<<idx_genjet1<<": "<<discrim1<<", "<<idx_genjet2<<": "<<discrim2<<endl;
    
    if((idx_genjet1 == -1 || idx_genjet2 == -1) || (idx_genjet1 == idx_genjet2)) return dummy_idx;
    
    RVec<int> finalgenjets_idx(2);
           
    if(GenJet_pt[idx_genjet1] > GenJet_pt[idx_genjet2]){ 
           finalgenjets_idx[0] = idx_genjet1;
           finalgenjets_idx[1] = idx_genjet2;
    }
    else{ 
           finalgenjets_idx[0] = idx_genjet2;
           finalgenjets_idx[1] = idx_genjet1;
    }
    return finalgenjets_idx;
}

bool atleast2Jets(rvec_i GenJet_idx){
    if(GenJet_idx.size() < 2) return false;
    else return true;
}

int IsGenMatched(int i, int j){
    if(i == j) return 1;
    else return 0;
}

RVec<RVec<float>> getTauSF(float SelectedTau_pt, float SelectedTau_eta, int SelectedTau_genPartFlav, bool IsMC, string year){
    RVec<float> vsJet, vsEle, vsMu;
    if (IsMC == false){
        vsJet.emplace_back(1.0);
        vsJet.emplace_back(1.0);
        vsJet.emplace_back(1.0);
        vsEle.emplace_back(1.0);
        vsEle.emplace_back(1.0);
        vsEle.emplace_back(1.0);
        vsMu.emplace_back(1.0);
        vsMu.emplace_back(1.0);
        vsMu.emplace_back(1.0);
        
        
    }
    
    else{
        string id;
        //std::string year = std::to_string(2017);
        // vs Jet
        id =  "DeepTau2017v2p1VSjet";
        //TString path = TString(remote_storage) + TString("data/tauSF/TauID_SF_pt_") + TString(id) + TString("_") + TString(year) + TString("ReReco") + TString(".root");
        //TFile *f = new TFile(path);
        //TFile *f =  TauID_SF_pt_DeepTau2017v2p1VSjet_2017ReReco;
        //else if (year == "2017") TFile *f =  TauID_SF_pt_DeepTau2017v2p1VSjet_2018ReReco;
        //TFile *f = TFile::Open(path);
        //TFile *f = new TFile();
        double_t pt = SelectedTau_pt;
        if (SelectedTau_genPartFlav==5){
            //TString path_down = TString(TString(vsJetwp) + TString("_down"));
            //TString path_cent = TString(TString(vsJetwp) + TString("_cent"));
            //TString path_up = TString(TString(vsJetwp) + TString("_up"));
            //TF1 * h_down = (TF1*)f->Get(path_down);
            //TF1 * h_cent =  (TF1*)f->Get(path_cent);
            //TF1 * h_up =  (TF1*)f->Get(path_up);
            //vsJet.emplace_back(h_down->Eval(pt));
            //vsJet.emplace_back(h_cent->Eval(pt));
            //vsJet.emplace_back(h_up->Eval(pt));
            vsJet.emplace_back(TauID_SF_pt_DeepTau2017v2p1VSjet_UL2017_h_down->Eval(pt));
            vsJet.emplace_back(TauID_SF_pt_DeepTau2017v2p1VSjet_UL2017_h_cent->Eval(pt));
            vsJet.emplace_back(TauID_SF_pt_DeepTau2017v2p1VSjet_UL2017_h_up->Eval(pt));
        }
        else{
            vsJet.emplace_back(1.0);
            vsJet.emplace_back(1.0);
            vsJet.emplace_back(1.0);
        }

        int bin;
        float sf, err;
        float eta = abs(SelectedTau_eta);  

        // vs ele
        id = "DeepTau2017v2p1VSe";
        //TString path_ele =  TString(remote_storage) + TString("data/tauSF/TauID_SF_eta_") + TString(id) + TString("_") + TString(year) + TString("ReReco") + TString(".root");
        //TFile *f_ele = new TFile(path_ele);
        //TFile *f_ele = TFile::Open(path_ele);

        //TFile *f_ele = TauID_SF_eta_DeepTau2017v2p1VSe_2017ReReco;
        //else if(year == "2018") TFile *f_ele = TauID_SF_eta_DeepTau2017v2p1VSe_2018ReReco;
        //TString histoname_ele = TString(vsElewp);
        //TH1F * hist = (TH1F *) f_ele->Get(histoname_ele);
        if (SelectedTau_genPartFlav == 1 || SelectedTau_genPartFlav == 3){
            //bin = hist->GetXaxis()->FindBin(eta);
            //sf  = hist->GetBinContent(bin);
            //err = hist->GetBinError(bin);
            bin = TauID_SF_eta_DeepTau2017v2p1VSe_UL2017_hist->GetXaxis()->FindBin(eta);
            sf  = TauID_SF_eta_DeepTau2017v2p1VSe_UL2017_hist->GetBinContent(bin);
            err = TauID_SF_eta_DeepTau2017v2p1VSe_UL2017_hist->GetBinError(bin);
            vsEle.emplace_back(sf-err);
            vsEle.emplace_back(sf);
            vsEle.emplace_back(sf+err);
        }
        else{
            vsEle.emplace_back(1.0);
            vsEle.emplace_back(1.0);
            vsEle.emplace_back(1.0);
        }

        //vs Mu
        id = "DeepTau2017v2p1VSmu";
        //TString path_mu =  TString(remote_storage) + TString("data/tauSF/TauID_SF_eta_") + TString(id) + TString("_") + TString(year)  + TString("ReReco")+ TString(".root");
        //TFile *f_mu = new TFile(path_mu);
        //TFile *f_mu = TFile::Open(path_mu);
        //TFile *f_mu = TauID_SF_eta_DeepTau2017v2p1VSmu_2017ReReco;
        //else if (year == "2018") TFile *f_mu = TauID_SF_eta_DeepTau2017v2p1VSmu_2018ReReco
        //TString histoname_mu = TString(vsMuwp);
        //TH1F * hist_mu = (TH1F *) f_mu->Get(histoname_mu);
        if (SelectedTau_genPartFlav == 2 || SelectedTau_genPartFlav == 4){
            //bin = hist_mu->GetXaxis()->FindBin(eta);
            //sf  = hist_mu->GetBinContent(bin);
            //err = hist_mu->GetBinError(bin);
            bin = TauID_SF_eta_DeepTau2017v2p1VSmu_UL2017_hist->GetXaxis()->FindBin(eta);
            sf  = TauID_SF_eta_DeepTau2017v2p1VSmu_UL2017_hist->GetBinContent(bin);
            err = TauID_SF_eta_DeepTau2017v2p1VSmu_UL2017_hist->GetBinError(bin);
            vsMu.emplace_back(sf-err);
            vsMu.emplace_back(sf);
            vsMu.emplace_back(sf+err);
        }
        else{
            vsMu.emplace_back(1.0);
            vsMu.emplace_back(1.0);
            vsMu.emplace_back(1.0);
        }
    }
    
    RVec<RVec<float>> result;
    result.emplace_back(vsJet);
    result.emplace_back(vsEle);
    result.emplace_back(vsMu);

    return result;
}

/*
RVec<float> getTES(float SelectedTau_pt, int SelectedTau_decayMode, int SelectedTau_genPartFlav, bool IsMC, string year){
    //cout<<"acscsascno"<<endl;
    //string year = "2017";
    string id = "DeepTau2017v2p1VSjet";
    
    float pt_low  = 34;
    float pt_high = 170;
    
    //TString path_low = TString(remote_storage) + TString("data/tauSF/TauES_dm_") + TString(id) + TString("_") + TString(year) + TString("ReReco") + TString(".root");
    //TString path_high = TString(remote_storage) + TString("data/tauSF/TauES_dm_") + TString(id) + TString("_") + TString(year) + TString("ReReco") + TString("_ptgt100.root");
    RVec<float> result(3);
    if(IsMC == false){
        //cout<<"babba"<<endl;
        result[0] = 1.;
        result[1] = 1.;
        result[2] = 1.;
    }
    else if((SelectedTau_decayMode == 0 || SelectedTau_decayMode == 1 || SelectedTau_decayMode == 10 || SelectedTau_decayMode == 11) && SelectedTau_genPartFlav == 5){ 
        //TFile *infile_low = TauES_dm_DeepTau2017v2p1VSjet_2017ReReco;
        //else if (year =="2018") TFile *infile_low = TauES_dm_DeepTau2017v2p1VSjet_2018ReReco;
        //TFile *infile_high = TauES_dm_DeepTau2017v2p1VSjet_2017ReReco_ptgt100;
        //else if (year =="2018") TFile *infile_low = TauES_dm_DeepTau2018v2p1VSjet_2018ReReco_ptgt100;        
        //TH1F * hist_low = (TH1F *) infile_low->Get("tes");
        //TH1F * hist_high = (TH1F *) infile_high->Get("tes");
        //int bin = hist_low->GetXaxis()->FindBin(SelectedTau_decayMode);
        //float tes = hist_low->GetBinContent(bin);
        //cout<<"ciarro"<<endl;
        int bin = TauES_dm_DeepTau2017v2p1VSjet_UL2017_hist_low->GetXaxis()->FindBin(SelectedTau_decayMode);
        float tes = TauES_dm_DeepTau2017v2p1VSjet_UL2017_hist_low->GetBinContent(bin);
        float err;
        if (SelectedTau_pt > pt_high){
            //cout<<"a"<<endl;
            //int bin_high = hist_high->GetXaxis()->FindBin(SelectedTau_decayMode);
            //float err = hist_high->GetBinError(bin_high);
            int bin_high = TauES_dm_DeepTau2017v2p1VSjet_UL2017_ptgt100_hist_high->GetXaxis()->FindBin(SelectedTau_decayMode);
            //float err = TauES_dm_DeepTau2017v2p1VSjet_2017ReReco_ptgt100_hist_high->GetBinError(bin_high);
            err = TauES_dm_DeepTau2017v2p1VSjet_UL2017_ptgt100_hist_high->GetBinError(bin_high);
        }
        else if (SelectedTau_pt > pt_low){
            //int bin_high = hist_high->GetXaxis()->FindBin(SelectedTau_decayMode);
            //float err_high = hist_high->GetBinError(bin_high);
            //float err_low  = hist_low->GetBinError(bin);
            //cout<<"b"<<endl;
            int bin_high = TauES_dm_DeepTau2017v2p1VSjet_UL2017_ptgt100_hist_high->GetXaxis()->FindBin(SelectedTau_decayMode);
            float err_high = TauES_dm_DeepTau2017v2p1VSjet_UL2017_ptgt100_hist_high->GetBinError(bin_high);
            float err_low  = TauES_dm_DeepTau2017v2p1VSjet_UL2017_hist_low->GetBinError(bin);
            //float err      = err_low + (err_high-err_low)/(pt_high-pt_low)*(SelectedTau_pt-pt_low);
            err      = err_low + (err_high-err_low)/(pt_high-pt_low)*(SelectedTau_pt-pt_low);
        }
        //else err = hist_low->GetBinError(bin);
        else err = TauES_dm_DeepTau2017v2p1VSjet_UL2017_hist_low->GetBinError(bin);
        
        result[0] = tes-err;
        result[1] = tes;
        result[2] = tes+err;
    }
    else{
        //cout<<"d"<<endl;
        result[0] = 1.;
        result[1] = 1.;
        result[2] = 1.;
    }
    //cout<<"oiknionono"<<endl;
    return result;
}
*/
RVec<float> getTES(rvec_f Tau_pt, rvec_i Tau_decayMode, rvec_i Tau_genPartFlav, bool IsMC, string year){
    string id = "DeepTau2017v2p1VSjet";
    
    float pt_low  = 34;
    float pt_high = 170;
    
    RVec<float> result_all;
    
    for (int i = 0; i < Tau_pt.size(); i++){
        RVec<float> result(3);
                
        int SelectedTau_decayMode = Tau_decayMode[i];
        float SelectedTau_pt = Tau_pt[i];
        int SelectedTau_genPartFlav = Tau_genPartFlav[i];
        
        if(IsMC == false){
            result[0] = 1.;
            result[1] = 1.;
            result[2] = 1.;
        }

        else if((SelectedTau_decayMode == 0 || SelectedTau_decayMode == 1 || SelectedTau_decayMode == 10 || SelectedTau_decayMode == 11) && SelectedTau_genPartFlav == 5){ 
            int bin = TauES_dm_DeepTau2017v2p1VSjet_UL2017_hist_low->GetXaxis()->FindBin(SelectedTau_decayMode);
            float tes = TauES_dm_DeepTau2017v2p1VSjet_UL2017_hist_low->GetBinContent(bin);
            float err;
            if (SelectedTau_pt > pt_high){
                int bin_high = TauES_dm_DeepTau2017v2p1VSjet_UL2017_ptgt100_hist_high->GetXaxis()->FindBin(SelectedTau_decayMode);
                err = TauES_dm_DeepTau2017v2p1VSjet_UL2017_ptgt100_hist_high->GetBinError(bin_high);
            }
            else if (SelectedTau_pt > pt_low){
                int bin_high = TauES_dm_DeepTau2017v2p1VSjet_UL2017_ptgt100_hist_high->GetXaxis()->FindBin(SelectedTau_decayMode);
                float err_high = TauES_dm_DeepTau2017v2p1VSjet_UL2017_ptgt100_hist_high->GetBinError(bin_high);
                float err_low  = TauES_dm_DeepTau2017v2p1VSjet_UL2017_hist_low->GetBinError(bin);
                err      = err_low + (err_high-err_low)/(pt_high-pt_low)*(SelectedTau_pt-pt_low);
            }
            else err = TauES_dm_DeepTau2017v2p1VSjet_UL2017_hist_low->GetBinError(bin);

            result[0] = tes-err;
            result[1] = tes;
            result[2] = tes+err;
        }
        else{
            result[0] = 1.;
            result[1] = 1.;
            result[2] = 1.;
        }
        result_all.emplace_back(result[0]);
        result_all.emplace_back(result[1]);
        result_all.emplace_back(result[2]);
    }
    return result_all;
}

/*
RVec<float> getFES(float SelectedTau_eta, int SelectedTau_decayMode, int SelectedTau_genPartFlav, bool IsMC, string year){

    //string year = "2017";
    string id = "DeepTau2017v2p1VSe";
    
    //TString path = TString(remote_storage) + TString("data/tauSF/TauFES_eta-dm_") + TString(id) + TString("_") + TString(year) + TString("ReReco") + TString(".root");
    RVec<float> result(3);
    if(IsMC == false){
        result[0] = 1.;
        result[1] = 1.;
        result[2] = 1.;
    }
    else if((SelectedTau_decayMode == 0 || SelectedTau_decayMode == 1) && (SelectedTau_genPartFlav == 1 || SelectedTau_genPartFlav == 3)){ 
        //TFile *infile = new TFile(path);
        //TFile *infile = TFile::Open(path);
        //TFile *infile = TauFES_eta_dm_DeepTau2017v2p1VSe_2017ReReco;
        //else if(year =="2018") infile = TauFES_eta-dm_DeepTau2017v2p1VSe_2018ReReco;
        //TGraphAsymmErrors * graph = (TGraphAsymmErrors *) infile->Get("fes");
        
        //float y = graph->GetY()[SelectedTau_decayMode];
        //float yup  = graph->GetErrorYhigh(SelectedTau_decayMode);
        //float ylow = graph->GetErrorYlow(SelectedTau_decayMode);
        
        //float y = TauFES_eta_dm_DeepTau2017v2p1VSe_2017ReReco_graph->GetY()[SelectedTau_decayMode];
        //float yup  = TauFES_eta_dm_DeepTau2017v2p1VSe_2017ReReco_graph->GetErrorYhigh(SelectedTau_decayMode);
        //float ylow = TauFES_eta_dm_DeepTau2017v2p1VSe_2017ReReco_graph->GetErrorYlow(SelectedTau_decayMode);
        
        int endcap_index_shift = 0;
        if(abs(SelectedTau_eta) >= 1.5) endcap_index_shift = 2;
        
        float y = TauFES_eta_dm_DeepTau2017v2p1VSe_UL2017_graph->GetY()[SelectedTau_decayMode + endcap_index_shift];
        float yup  = TauFES_eta_dm_DeepTau2017v2p1VSe_UL2017_graph->GetErrorYhigh(SelectedTau_decayMode + endcap_index_shift);
        float ylow = TauFES_eta_dm_DeepTau2017v2p1VSe_UL2017_graph->GetErrorYlow(SelectedTau_decayMode + endcap_index_shift);
        
        result[0] = y-ylow;
        result[1] = y;
        result[2] = y + yup;
    }
    else{
        result[0] = 1.;
        result[1] = 1.;
        result[2] = 1.;
    }
    return result;
}
*/

RVec<float> getFES(rvec_f Tau_eta, rvec_i Tau_decayMode, rvec_i Tau_genPartFlav, bool IsMC, string year){
    RVec<float> result_all;
    string id = "DeepTau2017v2p1VSe";
    for (int i = 0; i < Tau_eta.size(); i++){
        RVec<float> result(3);
                
        int SelectedTau_decayMode = Tau_decayMode[i];
        float SelectedTau_eta = Tau_eta[i];
        int SelectedTau_genPartFlav = Tau_genPartFlav[i];
        
        if(IsMC == false){
            result[0] = 1.;
            result[1] = 1.;
            result[2] = 1.;
        }
        
        else if((SelectedTau_decayMode == 0 || SelectedTau_decayMode == 1) && (SelectedTau_genPartFlav == 1 || SelectedTau_genPartFlav == 3)){ 

            int endcap_index_shift = 0;
            if(abs(SelectedTau_eta) >= 1.5) endcap_index_shift = 2;

            float y = TauFES_eta_dm_DeepTau2017v2p1VSe_UL2017_graph->GetY()[SelectedTau_decayMode + endcap_index_shift];
            float yup  = TauFES_eta_dm_DeepTau2017v2p1VSe_UL2017_graph->GetErrorYhigh(SelectedTau_decayMode + endcap_index_shift);
            float ylow = TauFES_eta_dm_DeepTau2017v2p1VSe_UL2017_graph->GetErrorYlow(SelectedTau_decayMode + endcap_index_shift);

            result[0] = y-ylow;
            result[1] = y;
            result[2] = y + yup;
        }
        else{
            result[0] = 1.;
            result[1] = 1.;
            result[2] = 1.;
        }
        result_all.emplace_back(result[0]);
        result_all.emplace_back(result[1]);
        result_all.emplace_back(result[2]);
    }
    return result_all;
}

float efficiency(int flv, float eta, float pt, string year){

    //TString path = TString(remote_storage) + TString("Btag_eff_") + TString(year) + TString(".root");
    //TFile *infile = new TFile(path);
    //TFile *infile = TFile::Open(path);
    TFile *infile = Btag_eff_UL2017;
    //else if(year =="2018") infile = Btag_eff_2018;
    TH2F * h;
    if(flv == 5){
        //TH2F * h = (TH2F *) infile->Get("h2_BTaggingEff_b")->CreateHistogram();
        //TEfficiency *eff = (TEfficiency *) infile->Get("h2_BTaggingEff_b");
        //h = (TH2F *) eff->CreateHistogram();
        h = Btag_eff_UL2017_h_b;
        
    }
    else if(flv == 4){
        //TEfficiency *eff = (TEfficiency *) infile->Get("h2_BTaggingEff_c");
        //h = (TH2F *) eff->CreateHistogram();
        h = Btag_eff_UL2017_h_c;
    }
    else{
        //h = (TH2F *) infile->Get("h2_BTaggingEff_udsg");
        //TEfficiency *eff = (TEfficiency *) infile->Get("h2_BTaggingEff_udsg");
        //h = (TH2F *) eff->CreateHistogram();
        h = Btag_eff_UL2017_h_udsg;
    }
    
    int binx = max(1, min(h->GetNbinsX(), h->GetXaxis()->FindBin(pt)));
    int biny = max(1, min(h->GetNbinsY(), h->GetYaxis()->FindBin(abs(eta))));
    
    return h->GetBinContent(binx,biny);
}
    
RVec<float> btagcalc(rvec_i GoodJets_idx, rvec_f Jet_pt, rvec_f Jet_eta, rvec_i Jet_partonFlavour, rvec_f Jet_btagDeepFlavB, rvec_f Jet_btagSF_deepjet_M_up, rvec_f Jet_btagSF_deepjet_M_down, rvec_f Jet_btagSF_deepjet_M, rvec_f Jet_btagDeepB, bool IsMC, string year){

    if (IsMC == false){
        RVec<float> result_dummy(5);
        result_dummy[0] = 1.;
        result_dummy[1] = 1.;
        result_dummy[2] = 1.;
        result_dummy[3] = 1.;
        result_dummy[4] = 1.;
    
    return result_dummy;
    
    }
    
    string tagger = "DeepFlv"; 
    string WP = "M";
    float threshold;
    //string year = "2017";
    
    float p_MC = 1.;
    float p_data = 1.;
    float p_data_btagUp = 1.;
    float p_data_btagDown = 1.;
    float p_data_mistagUp = 1.;
    float p_data_mistagDown = 1.;
    
    //map<string, float> WPbtagger = {
    //    {"DeepFlv_T", 0.7264}, 
    //    {"DeepFlv_M", 0.2770}, 
    //    {"DeepFlv_L", 0.0494}, 
    //    {"DeepCSV_T", 0.7527}, 
    //    {"DeepCSV_M", 0.4184}, 
    //    {"DeepCSV_L", 0.1241}
    //};


    for (size_t i = 0; i < GoodJets_idx.size(); i++) {
        int j = GoodJets_idx[i];
        
        if(tagger == "DeepFlv"){
            //threshold = WPbtagger[tagger + "_" + WP];
            threshold = 0.2770;
            if(Jet_btagDeepFlavB[j] >= threshold && Jet_pt[j] > BTAG_PT_CUT && abs(Jet_eta[j])<BTAG_ETA_CUT){
                p_MC =  p_MC * efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year);
                p_data = p_data * Jet_btagSF_deepjet_M[j] * efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year);
                if (abs(Jet_partonFlavour[j]) == 4 or abs(Jet_partonFlavour[j]) == 5){
                    p_data_btagUp = p_data_btagUp * Jet_btagSF_deepjet_M_up[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year);
                    p_data_mistagUp = p_data_mistagUp * Jet_btagSF_deepjet_M[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year);
                    p_data_btagDown = p_data_btagDown * Jet_btagSF_deepjet_M_down[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year);
                    p_data_mistagDown = p_data_mistagDown * Jet_btagSF_deepjet_M[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year);
                }
                else{
                    p_data_btagUp = p_data_btagUp *Jet_btagSF_deepjet_M[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year);
                    p_data_mistagUp = p_data_mistagUp * Jet_btagSF_deepjet_M_up[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year);
                    p_data_btagDown = p_data_btagDown * Jet_btagSF_deepjet_M[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year);
                    p_data_mistagDown = p_data_mistagDown *Jet_btagSF_deepjet_M_down[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year);
                }
            }  
            else if (Jet_btagDeepFlavB[j] < threshold && Jet_pt[j] > BTAG_PT_CUT && abs(Jet_eta[j])<BTAG_ETA_CUT){
                p_MC = p_MC *(1 - efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year));
                p_data = p_data *(1 - Jet_btagSF_deepjet_M[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year));
                if (abs(Jet_partonFlavour[j]) == 4 || abs(Jet_partonFlavour[j]) == 5){
                    p_data_btagUp =  p_data_btagUp * (1 - Jet_btagSF_deepjet_M_up[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year));
                    p_data_mistagUp = p_data_mistagUp * (1 - Jet_btagSF_deepjet_M[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year));
                    p_data_btagDown = p_data_btagDown * (1 - Jet_btagSF_deepjet_M_down[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year));
                    p_data_mistagDown = p_data_mistagDown * (1 - Jet_btagSF_deepjet_M[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year));
                }
                else{
                    p_data_btagUp = p_data_btagUp * (1 - Jet_btagSF_deepjet_M[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year));
                    p_data_mistagUp = p_data_mistagUp * (1 - Jet_btagSF_deepjet_M_up[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year));
                    p_data_btagDown = p_data_btagDown * (1 - Jet_btagSF_deepjet_M[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year));
                    p_data_mistagDown =  p_data_mistagDown * (1 - Jet_btagSF_deepjet_M_down[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year));
                }
            }
        }
        
        else if(tagger == "DeepCSV"){
            //threshold = WPbtagger[tagger + "_" + WP];
            threshold = 0.4184;
            if (Jet_btagDeepB[j] >= threshold){
                p_MC =  p_MC * efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year);
                p_data = p_data * Jet_btagSF_deepjet_M[j] * efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year);
                if (abs(Jet_partonFlavour[j]) == 4 || abs(Jet_partonFlavour[j]) == 5){
                    p_data_btagUp = p_data_btagUp * Jet_btagSF_deepjet_M_up[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year);
                    p_data_mistagUp = p_data_mistagUp * Jet_btagSF_deepjet_M[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year);
                    p_data_btagDown = p_data_btagDown * Jet_btagSF_deepjet_M_down[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year);
                    p_data_mistagDown = p_data_mistagDown * Jet_btagSF_deepjet_M[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year);
                }
                else{
                    p_data_btagUp = p_data_btagUp *Jet_btagSF_deepjet_M[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year);
                    p_data_mistagUp = p_data_mistagUp * Jet_btagSF_deepjet_M_up[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year);
                    p_data_btagDown = p_data_btagDown * Jet_btagSF_deepjet_M[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year);
                    p_data_mistagDown = p_data_mistagDown *Jet_btagSF_deepjet_M_down[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year);
                }
            }
            else if (Jet_btagDeepB[j] < threshold){
                p_MC = p_MC * (1 - efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year));
                p_data = p_data *(1 - Jet_btagSF_deepjet_M[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year));
                if (abs(Jet_partonFlavour[j]) == 4 || abs(Jet_partonFlavour[j]) == 5){
                    p_data_btagUp =  p_data_btagUp * (1 - Jet_btagSF_deepjet_M_up[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year));
                    p_data_mistagUp = p_data_mistagUp * (1 - Jet_btagSF_deepjet_M[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year));
                    p_data_btagDown = p_data_btagDown * (1 - Jet_btagSF_deepjet_M_down[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year));
                    p_data_mistagDown = p_data_mistagDown * (1 - Jet_btagSF_deepjet_M[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year));
                }
                else{
                    p_data_btagUp = p_data_btagUp * (1 - Jet_btagSF_deepjet_M[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year));
                    p_data_mistagUp = p_data_mistagUp * (1 - Jet_btagSF_deepjet_M_up[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year));
                    p_data_btagDown = p_data_btagDown * (1 - Jet_btagSF_deepjet_M[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year));
                    p_data_mistagDown =  p_data_mistagDown * (1 - Jet_btagSF_deepjet_M_down[j]*efficiency(abs(Jet_partonFlavour[j]), Jet_eta[j], Jet_pt[j], year));
                }
            }
        }
    }
    
    RVec<float> result(5);
    result[0] = p_data/p_MC;
    result[1] = p_data_btagUp/p_MC;
    result[2] = p_data_btagDown/p_MC;
    result[3] = p_data_mistagUp/p_MC;
    result[4] = p_data_mistagDown/p_MC;
    
    return result;
}

unordered_set<int> data_flags({339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, });


unordered_set<int> dataEle_flags({361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, });

unordered_set<int> dataMu_flags({339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360,});

unordered_map<int,float> xsecs({

{25,51.1}, //UL2016APV
{26,191.3},
{69,0.01538},
{18,3.697},
{19,0.5297},
{20,0.2529},
{21,0.4062},
{22,0.216},
{23,0.0758},
{28,12.178},
{50,0.2086},
{51,0.1651},
{52,0.05565},
{53,0.01398},
{54,0.2147},
{15,87.3348},
{64,27.59},
{67,6077.22},
{28,12.178},
{29,0.005724444444444444},
{30,0.005724444444444444},
{31,0.005725999999999999},
{32,0.005725999999999999},
{33,0.005725999999999999},
{34,0.005725999999999999},
{35,0.005725999999999999},
{36,0.005725999999999999},
{37,0.005725999999999999},
{38,35.85},
{40,1.0315},
{42,0.0118},
{43,2.7757},
{44,0.0896},
{45,0.237},
{46,0.212},
{47,0.952},
{1,0.564},
{2,13.74},
{11,0.003194},
{4,0.003194},
{5,0.003194},
{6,0.003194},
{8,0.003194},
{3,0.001586},
{7,0.001586},
{10,0.001586},
{71,0.002014},
{72,0.01036},
{73,0.01595},
{363,1.},
{364,1.},
{365,1.},
{366,1.},
{367,1.},
{341,1.},
{342,1.},
{343,1.},
{344,1.},
{345,1.},
    
    
{110,51.1}, //UL2016
{111,191.3},
{153,0.01538},
{103,3.697},
{104,0.5297},
{105,0.2529},
{106,0.4062},
{107,0.216},
{108,0.0758},
{113,12.178},
{134,0.2086},
{135,0.1651},
{136,0.05565},
{137,0.01398},
{138,0.2147},
{100,87.3348},
{148,27.59},
{151,6077.22},
{113,12.178},
{114,0.005724444444444444},
{115,0.005724444444444444},
{116,0.005725999999999999},
{117,0.005725999999999999},
{118,0.005725999999999999},
{119,0.005725999999999999},
{120,0.005725999999999999},
{121,0.005725999999999999},
{122,0.005725999999999999},
{123,35.85},
{125,1.0315},
{127,0.0118},
{128,2.7757},
{129,0.0896},
{130,0.237},
{131,0.212},
{132,0.952},
{86,0.9738},
{87,13.74},
{96,0.003194},
{89,0.003194},
{90,0.003194},
{91,0.003194},
{93,0.003194},
{88,0.001586},
{92,0.001586},
{95,0.001586},
{155,0.002014},
{156,0.01036},
{157,0.01595},
{369,1.},
{370,1.},
{371,1.},
{347,1.},
{348,1.},
{349,1.},
    
{194,51.1}, //UL2017
{195,191.3},
{238,0.01538},
{187,3.697},
{188,0.5297},
{189,0.2529},
{190,0.4062},
{191,0.216},
{192,0.0758},
{197,12.178},
{219,0.2086},
{220,0.1651},
{221,0.05565},
{222,0.01398},
{223,0.2147},
{184,87.3348},
{233,27.59},
{236,6077.22},
{197,12.178},
{198,0.005724444444444444},
{199,0.005724444444444444},
{200,0.005725999999999999},
{201,0.005725999999999999},
{202,0.005725999999999999},
{203,0.005725999999999999},
{204,0.005725999999999999},
{205,0.005725999999999999},
{206,0.005725999999999999},
{207,35.85},
{209,1.0315},
{211,0.0118},
{212,2.7757},
{213,0.0896},
{214,0.237},
{215,0.212},
{216,0.952},
{170,0.9738},
{171,13.74},
{180,0.003194},
{173,0.003194},
{174,0.003194},
{175,0.003194},
{177,0.003194},
{172,0.001586},
{176,0.001586},
{179,0.001586},
{240,0.002014},
{241,0.01036},
{242,0.01595},
{373,1.},
{374,1.},
{375,1.},
{376,1.},
{377,1.},
{351,1.},
{352,1.},
{353,1.},
{354,1.},
{355,1.}, 
    
{279,51.1}, //UL2018
{280,191.3},
{323,0.01538},
{272,3.697},
{273,0.5297},
{274,0.2529},
{275,0.4062},
{276,0.216},
{277,0.0758},
{282,12.178},
{304,0.2086},
{305,0.1651},
{306,0.05565},
{307,0.01398},
{308,0.2147},
{269,87.3348},
{318,27.59},
{321,6077.22},
{282,12.178},
{283,0.005724444444444444},
{284,0.005724444444444444},
{285,0.005725999999999999},
{286,0.005725999999999999},
{287,0.005725999999999999},
{288,0.005725999999999999},
{289,0.005725999999999999},
{290,0.005725999999999999},
{291,0.005725999999999999},
{292,35.85},
{294,1.0315},
{296,0.0118},
{297,2.7757},
{298,0.0896},
{299,0.237},
{300,0.212},
{301,0.952},
{255,0.9738},
{256,13.74},
{265,0.003194},
{258,0.003194},
{259,0.003194},
{260,0.003194},
{262,0.003194},
{257,0.001586},
{261,0.001586},
{264,0.001586},
{325,0.002014},
{326,0.01036},
{327,0.01595},
{379,1.},
{380,1.},
{381,1.},
{382,1.},
{357,1.},
{358,1.},
{359,1.},
{360,1.},
    });
   
unordered_map<int,float> Nevents({
{194,29890946}, //UL2017
{195,60212926},
{187,3534208},
{188,13822000},
{189,14036000},
{190,655018},
{191,2891483},
{192,9530000},
{197,7098000},
{219,9854000},
{220,178000},
{221,9898000},
{222,9524000},
{223,1736000},
{184,103204000},
{233,7889000},
{236,195529774},
{197,7098000},
{198,3780000},
{199,3876000},
{200,3998000},
{201,3962000},
{202,3978000},
{203,3985000},
{204,3840000},
{205,3984000},
{206,1910000},
{207,5649000},
{209,6828983},
{211,998000},
{212,12974000},
{213,6440000},
{214,2811630},
{215,5070989},
{216,869559},
{170,40839000},
{171,41708429},
{180,497000},
{173,500000},
{174,498000},
{175,500000},
{177,500000},
{172,997000},
{176,975090},
{179,499000},
{240,1978000},
{241,1952000},
{242,1902000},
{373,1},
{374,1},
{375,1},
{376,1},
{377,1},
{351,1},
{352,1},
{353,1},
{354,1},
{355,1},
});

bool isMC(int SampleFlag){
    bool is_in = data_flags.find(SampleFlag) != data_flags.end();
    if (is_in == true) return false;
    else return true;
}

float getLumi(string Year, bool IsMC){
    if (IsMC == false) return 1.;
    else if(Year == "UL2016" || Year == "UL2016APV") return -999;
    //else if (Year == 2017) return 41.48;
    else if (Year == "UL2017") return 41.54;
    else return 59.83;
}

float getXSec(int Sample, bool IsMC){
    if (IsMC == false) return 1.;
    else{
        return  xsecs[Sample];
    }
}

float getNevents(int Sample, bool IsMC){
    if (IsMC == false) return 1.;
    else{
        return  Nevents[Sample];
    }
}


///// old
////questa è vecchia TMVA::Experimental::RBDT<> bdt("SMxgb", "https://ttedesch.web.cern.ch/ttedesch/SMxgb.root");
//TMVA::Experimental::RBDT<> bdt("xgb_SM_v100", "https://ttedesch.web.cern.ch/ttedesch/VBS_ML_v4/xgb_SM_v100.root");
//float SMinference(float m_jj, float m_jjtaulep, float m_taulep, float mT_lep_MET, float leadjet_pt, float subleadjet_pt, float tau_mass, float MET_pt){
//    auto y1 = bdt.Compute({m_jj, m_jjtaulep, m_taulep, mT_lep_MET, leadjet_pt, subleadjet_pt, tau_mass, MET_pt});
//    return y1[0];
//}

int CountBJets(rvec_f Jet_pt, rvec_f Jet_eta, rvec_f Jet_btagDeepFlavB){
    int nb=0;
    for (int i = 0; i < Jet_pt.size(); i++) {
        if (Jet_btagDeepFlavB[i]>=BTAG_WP_VALUE && Jet_pt[i]>BTAG_PT_CUT && abs(Jet_eta[i])<BTAG_ETA_CUT) nb++;
    }
    return nb;
}

float taujet_RelPt(int selectedtau_jetIdx, float selectedtau_pt, rvec_f Jet_pt){
    if (selectedtau_jetIdx == -1) return -999.;
    else return Jet_pt[selectedtau_jetIdx]/selectedtau_pt;
}

float taujet_deltaPhi(int selectedtau_jetIdx, float selectedtau_phi, rvec_f Jet_phi){
    if (selectedtau_jetIdx == -1) return -999.;
    else return deltaPhi(selectedtau_phi, Jet_phi[selectedtau_jetIdx]);
}

//to_keep = ['leadjet_DeepFlv_b', 'subleadjet_DeepFlv_b', 'event_Zeppenfeld_over_deltaEta_jj', 'taujet_relpt', 'm_o1', 'event_RT', 'taujet_deltaPhi', 'nJets', 'mT_lep_MET', 'm_jjtau', 'm_1T', 'nBJets', 'subleadjet_pt', 'm_jjtaulep', 'm_jj', 'leadjet_pt']
//TMVA::Experimental::RBDT<> bdt("xgb_SM_v100", "https://ttedesch.web.cern.ch/ttedesch/VBS_ML_v4/optimized_model_SM_clean.root");
TMVA::Experimental::RBDT<> bdt("xgb_SM_v100", "https://vbs-pg-support.web.cern.ch/models/optimized_model_SM_clean.root");

float SMinference(float leadjet_DeepFlv_b, float subleadjet_DeepFlv_b, float event_Zeppenfeld_over_deltaEta_jj, float taujet_relpt, float m_o1, float event_RT, float taujet_deltaPhi, float nJets, float mT_lep_MET, float m_jjtau, float m_1T, float nBJets, float subleadjet_pt, float m_jjtaulep, float m_jj, float leadjet_pt){
    auto y1 = bdt.Compute({leadjet_DeepFlv_b, subleadjet_DeepFlv_b, event_Zeppenfeld_over_deltaEta_jj, taujet_relpt, m_o1, event_RT, taujet_deltaPhi, nJets, mT_lep_MET, m_jjtau, m_1T, nBJets, subleadjet_pt, m_jjtaulep, m_jj, leadjet_pt});
    return y1[0];
}

/*
R__LOAD_LIBRARY(/usr/local/lib/libtensorflow.so)
#include "cppflow/ops.h"
#include "cppflow/model.h"
cppflow::model model("./dnn_optimized_model_tagger_mjjabsdeltaetajj");

    

float DNNinference(float x1, float x2, float x3, float x4, float x5, float x6, float x7, float x8, float x9, float x10, float x11, float x12, float x13, float x14, float x15, float x16, float x17, float x18, float x19, float x20, float x21, float x22, float x23, float x24, float x25, float x26, float x27, float x28){
    auto input = cppflow::fill({1, 28}, 1.0f);
    auto output = model(input);
    auto values = output.get_data<float>();
    return values[0];
}
*/
/*
float DNNinference(float x1, float x2, float x3, float x4, float x5, float x6, float x7, float x8, float x9, float x10, float x11, float x12, float x13, float x14, float x15, float x16, float x17, float x18, float x19, float x20){
     TMVA::Reader *reader = new TMVA::Reader();
    Float_t var1, var2, var3, var4;
    reader->AddVariable( "var1", &x1 );
    reader->AddVariable( "var2", &x2 );
    reader->AddVariable( "var3", &x3 );
    reader->AddVariable( "var4", &x4 );
    reader->AddVariable( "var5", &x5 );
    reader->AddVariable( "var6", &x6 );
    reader->AddVariable( "var7", &x7 );
    reader->AddVariable( "var8", &x8 );
    reader->AddVariable( "var9", &x9 );
    reader->AddVariable( "var10", &x10 );
    reader->AddVariable( "var11", &x11 );
    reader->AddVariable( "var12", &x12 );
    reader->AddVariable( "var13", &x13 );
    reader->AddVariable( "var14", &x14 );
    reader->AddVariable( "var15", &x15 );
    reader->AddVariable( "var16", &x16 );
    reader->AddVariable( "var17", &x17 );
    reader->AddVariable( "var18", &x18 );
    reader->AddVariable( "var19", &x19 );
    reader->AddVariable( "var20", &x20 );
    
    reader->BookMVA('PyKeras', TString('dataset/weights/TMVAClassification_PyKeras.weights.xml'))
        
    auto y1 = bdt.Compute({m_jj, m_jjtaulep, m_taulep, mT_lep_MET, leadjet_pt, subleadjet_pt, tau_mass, MET_pt});
    return y1[0];
}
*/

RVec<int> GetGenMatched(rvec_i Jet_genJetIdx, rvec_i GenJet_idx){
    RVec<int> GenMatched_idx(2);
    for (int i = 0; i < Jet_genJetIdx.size(); i++) {
        if(Jet_genJetIdx[i] == GenJet_idx[0])  GenMatched_idx[0] = i;
        else if(Jet_genJetIdx[i] == GenJet_idx[1])  GenMatched_idx[1] = i;
    }
    return GenMatched_idx;
}


bool DataLeptonCheck(int SampleFlag, int GoodLeptonFamily, bool isMC){
    if(isMC == false){
        bool is_in_ele = dataEle_flags.find(SampleFlag) != dataEle_flags.end();
        if(is_in_ele == true && GoodLeptonFamily == 1) return false;
        bool is_in_mu = dataMu_flags.find(SampleFlag) != dataMu_flags.end();
        if(is_in_mu == true && GoodLeptonFamily == 0) return false;
    }
    return true;
}