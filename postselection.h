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

const float PT_CUT_MU=  35;
const float ETA_CUT_MU= 2.4;
const float ISO_CUT_MU= 0.15;

//const size_t PT_CUT_ELE=  35;
const float PT_CUT_ELE=  35;
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

const string remote_storage = "https://ttedesch.web.cern.ch/ttedesch/nanoAOD-tools/python/postprocessing/";

TFile *TauID_SF_pt_DeepTau2017v2p1VSjet_2017ReReco = TFile::Open(TString(remote_storage) + TString("data/tauSF/TauID_SF_pt_") + TString("DeepTau2017v2p1VSjet") + TString("_") + TString("2017") + TString("ReReco") + TString(".root"));
TString path_down = TString(TString(vsJetwp) + TString("_down"));
TString path_cent = TString(TString(vsJetwp) + TString("_cent"));
TString path_up = TString(TString(vsJetwp) + TString("_up"));
TF1 * TauID_SF_pt_DeepTau2017v2p1VSjet_2017ReReco_h_down = (TF1*)TauID_SF_pt_DeepTau2017v2p1VSjet_2017ReReco->Get(path_down);
TF1 * TauID_SF_pt_DeepTau2017v2p1VSjet_2017ReReco_h_cent =  (TF1*)TauID_SF_pt_DeepTau2017v2p1VSjet_2017ReReco->Get(path_cent);
TF1 * TauID_SF_pt_DeepTau2017v2p1VSjet_2017ReReco_h_up =  (TF1*)TauID_SF_pt_DeepTau2017v2p1VSjet_2017ReReco->Get(path_up);

TFile *TauID_SF_eta_DeepTau2017v2p1VSe_2017ReReco = TFile::Open(TString(remote_storage) + TString("data/tauSF/TauID_SF_eta_") + TString("DeepTau2017v2p1VSe") + TString("_") + TString("2017") + TString("ReReco") + TString(".root"));
TString histoname_ele = TString(vsElewp);
TH1F * TauID_SF_eta_DeepTau2017v2p1VSe_2017ReReco_hist = (TH1F *) TauID_SF_eta_DeepTau2017v2p1VSe_2017ReReco->Get(histoname_ele);

TFile *TauID_SF_eta_DeepTau2017v2p1VSmu_2017ReReco = TFile::Open(TString(remote_storage) + TString("data/tauSF/TauID_SF_eta_") + TString("DeepTau2017v2p1VSmu") + TString("_") + TString("2017")  + TString("ReReco")+ TString(".root"));
TString histoname_mu = TString(vsMuwp);
TH1F * TauID_SF_eta_DeepTau2017v2p1VSmu_2017ReReco_hist = (TH1F *) TauID_SF_eta_DeepTau2017v2p1VSmu_2017ReReco->Get(histoname_mu);

TFile *TauES_dm_DeepTau2017v2p1VSjet_2017ReReco = TFile::Open(TString(remote_storage) + TString("data/tauSF/TauES_dm_") + TString("DeepTau2017v2p1VSjet") + TString("_") + TString("2017") + TString("ReReco") + TString(".root"));        
TH1F * TauES_dm_DeepTau2017v2p1VSjet_2017ReReco_hist_low = (TH1F *) TauES_dm_DeepTau2017v2p1VSjet_2017ReReco->Get("tes");

TFile *TauES_dm_DeepTau2017v2p1VSjet_2017ReReco_ptgt100 = TFile::Open(TString(remote_storage) + TString("data/tauSF/TauES_dm_") + TString("DeepTau2017v2p1VSjet") + TString("_") + TString("2017") + TString("ReReco") + TString("_ptgt100.root"));
TH1F * TauES_dm_DeepTau2017v2p1VSjet_2017ReReco_ptgt100_hist_high = (TH1F *) TauES_dm_DeepTau2017v2p1VSjet_2017ReReco_ptgt100->Get("tes");

TFile *TauFES_eta_dm_DeepTau2017v2p1VSe_2017ReReco = TFile::Open(TString(remote_storage) + TString("data/tauSF/TauFES_eta-dm_") + TString("DeepTau2017v2p1VSe") + TString("_") + TString("2017") + TString("ReReco") + TString(".root"));
TGraphAsymmErrors * TauFES_eta_dm_DeepTau2017v2p1VSe_2017ReReco_graph = (TGraphAsymmErrors *) TauFES_eta_dm_DeepTau2017v2p1VSe_2017ReReco->Get("fes");

TFile *Btag_eff_2017 = TFile::Open(TString(remote_storage) + TString("Btag_eff_") + TString("2017") + TString(".root"));
TEfficiency *eff_b = (TEfficiency *) Btag_eff_2017->Get("h2_BTaggingEff_b");
TH2F *Btag_eff_2017_h_b = (TH2F *) eff_b->CreateHistogram();
TEfficiency *eff_c = (TEfficiency *) Btag_eff_2017->Get("h2_BTaggingEff_c");
TH2F *Btag_eff_2017_h_c = (TH2F *) eff_c->CreateHistogram();
TEfficiency *eff_udsg = (TEfficiency *) Btag_eff_2017->Get("h2_BTaggingEff_udsg");
TH2F *Btag_eff_2017_h_udsg = (TH2F *) eff_udsg->CreateHistogram();

//TFile *FR_vsjet2_vsmuT_ZZ = TFile::Open(TString(remote_storage) + TString("FR_vsjet2_vsmuT_ZZ.root"));
TFile *FR_vsjet2_vsmuT_ZZ = TFile::Open(TString(remote_storage) + TString("FR_vsjet2_2017.root"));
TH2F *FR_vsjet2_vsmuT_ZZ_histo_ele = (TH2F*)FR_vsjet2_vsmuT_ZZ->Get("hFRDataeledif");
TH2F *FR_vsjet2_vsmuT_ZZ_histo_mu = (TH2F*)FR_vsjet2_vsmuT_ZZ->Get("hFRDatamudif");
TH2F *FR_vsjet2_vsmuT_ZZ_histo_tau = (TH2F*)FR_vsjet2_vsmuT_ZZ->Get("hFRDatataudif");
                                        
//TFile *FR_vsjet4_vsmuT_ZZ = TFile::Open(TString(remote_storage) + TString("FR_vsjet4_vsmuT_ZZ.root"));
TFile *FR_vsjet4_vsmuT_ZZ = TFile::Open(TString(remote_storage) + TString("FR_vsjet4_2017.root"));
TH2F *FR_vsjet4_vsmuT_ZZ_histo_ele = (TH2F*)FR_vsjet4_vsmuT_ZZ->Get("hFRDataeledif");
TH2F *FR_vsjet4_vsmuT_ZZ_histo_mu = (TH2F*)FR_vsjet4_vsmuT_ZZ->Get("hFRDatamudif");
TH2F *FR_vsjet4_vsmuT_ZZ_histo_tau = (TH2F*)FR_vsjet4_vsmuT_ZZ->Get("hFRDatataudif");

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

RVec<int> SelectElectron(rvec_f lepton_pt, rvec_f lepton_eta, rvec_f lepton_phi, rvec_f lepton_jetRelIso, rvec_b lepton_mvaFall17V2Iso_WPL, rvec_f lepton_mvaFall17V2Iso_WP90, rvec_f jet_eta, rvec_f jet_phi, rvec_i VBSJets_idx){
    //setting jet-related quantities if isolation from them is needed
    //const auto jet1_idx = VBSJets_idx[0];
    //const auto jet2_idx = VBSJets_idx[1];

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


RVec<int> SelectAndVetoTaus(rvec_f Tau_pt, rvec_f Tau_eta, rvec_f Tau_phi, RVec<UChar_t> &Tau_idDeepTau2017v2p1VSjet, RVec<UChar_t> &Tau_idDeepTau2017v2p1VSe, RVec<UChar_t> &Tau_idDeepTau2017v2p1VSmu,rvec_b Tau_idDecayModeNewDMs, int GoodLeptonFamily, rvec_i Electron_idx, rvec_f Electron_eta, rvec_f Electron_phi, rvec_i Muon_idx, rvec_f Muon_eta, rvec_f Muon_phi, rvec_f Jet_eta, rvec_f Jet_phi, rvec_i VBSJet_idx)
{
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
            cutloose_vsjet = ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_ELE;
            if ((Tau_idDeepTau2017v2p1VSjet[i]>=cutloose_vsjet && Tau_idDeepTau2017v2p1VSe[i]>=ID_TAU_RECO_DEEPTAU_VSELE && Tau_idDeepTau2017v2p1VSmu[i]>=ID_TAU_RECO_DEEPTAU_VSMU && Tau_idDecayModeNewDMs[i]) && deltaR(Tau_eta[i], Tau_phi[i], Electron_eta[Electron_idx[0]], Electron_phi[Electron_idx[0]])>DR_OVERLAP_CONE_TAU && deltaR(Tau_eta[i], Tau_phi[i], jet1eta, jet1phi)>isocone && deltaR(Tau_eta[i], Tau_phi[i], jet2eta, jet2phi)>isocone && Tau_pt[i]>=PT_CUT_TAU && abs(Tau_eta[i])<=ETA_CUT_TAU){
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
            cutloose_vsjet = ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_MU;
            if ((Tau_idDeepTau2017v2p1VSjet[i]>=cutloose_vsjet && Tau_idDeepTau2017v2p1VSe[i]>=ID_TAU_RECO_DEEPTAU_VSELE && Tau_idDeepTau2017v2p1VSmu[i]>=ID_TAU_RECO_DEEPTAU_VSMU && Tau_idDecayModeNewDMs[i]) && deltaR(Tau_eta[i], Tau_phi[i], Muon_eta[Muon_idx[0]], Muon_phi[Muon_idx[0]])>DR_OVERLAP_CONE_TAU && deltaR(Tau_eta[i], Tau_phi[i], jet1eta, jet1phi)>isocone && deltaR(Tau_eta[i], Tau_phi[i], jet2eta, jet2phi)>isocone && Tau_pt[i]>=PT_CUT_TAU && abs(Tau_eta[i])<=ETA_CUT_TAU){
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

RVec<RVec<float>> getTauSF(float SelectedTau_pt, float SelectedTau_eta, int SelectedTau_genPartFlav, bool IsMC){
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
        std::string year = std::to_string(2017);
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
            vsJet.emplace_back(TauID_SF_pt_DeepTau2017v2p1VSjet_2017ReReco_h_down->Eval(pt));
            vsJet.emplace_back(TauID_SF_pt_DeepTau2017v2p1VSjet_2017ReReco_h_cent->Eval(pt));
            vsJet.emplace_back(TauID_SF_pt_DeepTau2017v2p1VSjet_2017ReReco_h_up->Eval(pt));
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
            bin = TauID_SF_eta_DeepTau2017v2p1VSe_2017ReReco_hist->GetXaxis()->FindBin(eta);
            sf  = TauID_SF_eta_DeepTau2017v2p1VSe_2017ReReco_hist->GetBinContent(bin);
            err = TauID_SF_eta_DeepTau2017v2p1VSe_2017ReReco_hist->GetBinError(bin);
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
            bin = TauID_SF_eta_DeepTau2017v2p1VSmu_2017ReReco_hist->GetXaxis()->FindBin(eta);
            sf  = TauID_SF_eta_DeepTau2017v2p1VSmu_2017ReReco_hist->GetBinContent(bin);
            err = TauID_SF_eta_DeepTau2017v2p1VSmu_2017ReReco_hist->GetBinError(bin);
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


RVec<float> getTES(float SelectedTau_pt, int SelectedTau_decayMode, int SelectedTau_genPartFlav, bool IsMC){
    //cout<<"acscsascno"<<endl;
    string year = "2017";
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
        int bin = TauES_dm_DeepTau2017v2p1VSjet_2017ReReco_hist_low->GetXaxis()->FindBin(SelectedTau_decayMode);
        float tes = TauES_dm_DeepTau2017v2p1VSjet_2017ReReco_hist_low->GetBinContent(bin);
        float err;
        if (SelectedTau_pt > pt_high){
            //cout<<"a"<<endl;
            //int bin_high = hist_high->GetXaxis()->FindBin(SelectedTau_decayMode);
            //float err = hist_high->GetBinError(bin_high);
            int bin_high = TauES_dm_DeepTau2017v2p1VSjet_2017ReReco_ptgt100_hist_high->GetXaxis()->FindBin(SelectedTau_decayMode);
            //float err = TauES_dm_DeepTau2017v2p1VSjet_2017ReReco_ptgt100_hist_high->GetBinError(bin_high);
            err = TauES_dm_DeepTau2017v2p1VSjet_2017ReReco_ptgt100_hist_high->GetBinError(bin_high);
        }
        else if (SelectedTau_pt > pt_low){
            //int bin_high = hist_high->GetXaxis()->FindBin(SelectedTau_decayMode);
            //float err_high = hist_high->GetBinError(bin_high);
            //float err_low  = hist_low->GetBinError(bin);
            //cout<<"b"<<endl;
            int bin_high = TauES_dm_DeepTau2017v2p1VSjet_2017ReReco_ptgt100_hist_high->GetXaxis()->FindBin(SelectedTau_decayMode);
            float err_high = TauES_dm_DeepTau2017v2p1VSjet_2017ReReco_ptgt100_hist_high->GetBinError(bin_high);
            float err_low  = TauES_dm_DeepTau2017v2p1VSjet_2017ReReco_hist_low->GetBinError(bin);
            //float err      = err_low + (err_high-err_low)/(pt_high-pt_low)*(SelectedTau_pt-pt_low);
            err      = err_low + (err_high-err_low)/(pt_high-pt_low)*(SelectedTau_pt-pt_low);
        }
        //else err = hist_low->GetBinError(bin);
        else err = TauES_dm_DeepTau2017v2p1VSjet_2017ReReco_hist_low->GetBinError(bin);
        
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

//RVec<float> getFES(int SelectedTau_decayMode, int SelectedTau_genPartFlav, bool IsMC){
RVec<float> getFES(float SelectedTau_eta, int SelectedTau_decayMode, int SelectedTau_genPartFlav, bool IsMC){

    string year = "2017";
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
        
        float y = TauFES_eta_dm_DeepTau2017v2p1VSe_2017ReReco_graph->GetY()[SelectedTau_decayMode + endcap_index_shift];
        float yup  = TauFES_eta_dm_DeepTau2017v2p1VSe_2017ReReco_graph->GetErrorYhigh(SelectedTau_decayMode + endcap_index_shift);
        float ylow = TauFES_eta_dm_DeepTau2017v2p1VSe_2017ReReco_graph->GetErrorYlow(SelectedTau_decayMode + endcap_index_shift);
        
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


float efficiency(int flv, float eta, float pt, string year){

    //TString path = TString(remote_storage) + TString("Btag_eff_") + TString(year) + TString(".root");
    //TFile *infile = new TFile(path);
    //TFile *infile = TFile::Open(path);
    TFile *infile = Btag_eff_2017;
    //else if(year =="2018") infile = Btag_eff_2018;
    TH2F * h;
    if(flv == 5){
        //TH2F * h = (TH2F *) infile->Get("h2_BTaggingEff_b")->CreateHistogram();
        //TEfficiency *eff = (TEfficiency *) infile->Get("h2_BTaggingEff_b");
        //h = (TH2F *) eff->CreateHistogram();
        h = Btag_eff_2017_h_b;
        
    }
    else if(flv == 4){
        //TEfficiency *eff = (TEfficiency *) infile->Get("h2_BTaggingEff_c");
        //h = (TH2F *) eff->CreateHistogram();
        h = Btag_eff_2017_h_c;
    }
    else{
        //h = (TH2F *) infile->Get("h2_BTaggingEff_udsg");
        //TEfficiency *eff = (TEfficiency *) infile->Get("h2_BTaggingEff_udsg");
        //h = (TH2F *) eff->CreateHistogram();
        h = Btag_eff_2017_h_udsg;
    }
    
    int binx = max(1, min(h->GetNbinsX(), h->GetXaxis()->FindBin(pt)));
    int biny = max(1, min(h->GetNbinsY(), h->GetYaxis()->FindBin(abs(eta))));
    
    return h->GetBinContent(binx,biny);
}
    
RVec<float> btagcalc(rvec_i GoodJets_idx, rvec_f Jet_pt, rvec_f Jet_eta, rvec_i Jet_partonFlavour, rvec_f Jet_btagDeepFlavB, rvec_f Jet_btagSF_deepjet_M_up, rvec_f Jet_btagSF_deepjet_M_down, rvec_f Jet_btagSF_deepjet_M, rvec_f Jet_btagDeepB, bool IsMC){

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
    string year = "2017";
    
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

unordered_set<int> data_flags({201, 202, 203, 204, 205, 206, 207, 212, 213, 214, 215, 216, 217, 218, 223, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 429, 430, 431, 432, 433, 438, 439, 440, 441, 442, 446, 447, 448, 449, 450});

unordered_set<int> dataEle_flags({212, 213, 214, 215, 216, 217, 218, 438, 439, 440, 441, 442});

unordered_set<int> dataMu_flags({201, 202, 203, 204, 205, 206, 207, 429, 430, 431, 432, 433});

/*
unordered_map<int,float> xsecs({
    {268,0.1191},
    {19,0.1191},
    {432,1.0},
    {203,1.0},
    {5,0.01595},
    {315,1.0},
    {66,1.0},
    {254,0.01595},
    {423,0.05565},
    {196,0.05565},
    {262,1.0},
    {13,1.0},
    {163,0.5439},
    {92,1.0},
    {122,3.185},
    {371,3.185},
    {341,1.0},
    {396,0.5439},
    {335,0.1191},
    {240,1.0},
    {447,1.0},
    {86,0.1191},
    {299,0.1191},
    {52,1.0},
    {301,1.0},
    {50,0.1191},
    {378,377.96},
    {129,377.96},
    {24,1.0},
    {338,0.1191},
    {149,0.0389136},
    {89,0.1191},
    {400,0.2432},
    {167,0.2432},
    {174,22.81},
    {219,1.0},
    {141,1.0},
    {14,1.0},
    {263,1.0},
    {63,0.1191},
    {312,0.1191},
    {148,1.60809},
    {369,3.185},
    {238,1.0},
    {325,0.1191},
    {303,0.1191},
    {54,0.1191},
    {76,0.1191},
    {120,3.185},
    {12,0.002014},
    {234,1.0},
    {448,1.0},
    {261,0.002014},
    {115,0.5644},
    {306,1.0},
    {11,0.002014},
    {278,1.0},
    {29,1.0},
    {118,14.93},
    {217,1.0},
    {36,0.1191},
    {285,0.1191},
    {21,1.0},
    {411,28.87},
    {184,28.87},
    {157,1012.0},
    {414,3.879},
    {187,3.879},
    {323,1.0},
    {74,1.0},
    {402,0.2149},
    {28,0.1191},
    {171,1.0},
    {277,0.1191},
    {272,0.1191},
    {320,0.1191},
    {419,1.0},
    {169,0.2149},
    {192,1.0},
    {158,330.4},
    {23,0.1191},
    {424,0.013989},
    {354,1.0},
    {105,1.0},
    {197,0.013989},
    {344,0.1191},
    {416,0.5269},
    {189,0.5269},
    {95,0.1191},
    {59,0.1191},
    {446,1.0},
    {258,1.0},
    {9,1.0},
    {232,1.0},
    {308,0.1191},
    {42,1.0},
    {267,0.1191},
    {150,1.0},
    {291,1.0},
    {336,1.0},
    {351,1.0},
    {87,1.0},
    {53,0.1191},
    {374,1.404},
    {49,0.1191},
    {298,0.1191},
    {302,0.1191},
    {218,1.0},
    {37,0.1191},
    {172,11.08},
    {439,1.0},
    {295,0.1191},
    {69,1.0},
    {55,0.1191},
    {71,0.1191},
    {215,1.0},
    {227,1.0},
    {442,1.0},
    {244,1.0},
    {355,347700},
    {286,0.1191},
    {3,0.002014},
    {252,0.002014},
    {106,23700000.0},
    {412,28.87},
    {130,88.287},
    {426,88.287},
    {185,28.87},
    {114,1.0},
    {133,1637.13},
    {363,1.0},
    {284,0.1191},
    {313,0.1191},
    {222,1.0},
    {421,0.2086},
    {194,0.2086},
    {445,1.0},
    {177,22.81},
    {383,59.1811},
    {35,0.1191},
    {228,1.0},
    {413,30.52},
    {273,1.0},
    {367,14.93},
    {385,6.65621},
    {138,6.65621},
    {136,59.1811},
    {72,0.1191},
    {311,0.1191},
    {17,0.1191},
    {266,0.1191},
    {62,0.1191},
    {321,0.1191},
    {347,0.1191},
    {201,1.0},
    {180,22.81},
    {358,32100},
    {109,29980.0},
    {429,1.0},
    {98,0.1191},
    {231,1.0},
    {236,1.0},
    {264,1.0},
    {450,1.0},
    {393,1.0},
    {293,0.1191},
    {410,28.87},
    {44,0.1191},
    {247,1.0},
    {193,0.1703},
    {126,0.7882},
    {324,1.0},
    {183,28.87},
    {125,1.404},
    {75,1.0},
    {375,0.7882},
    {420,0.1703},
    {206,1.0},
    {399,0.5104},
    {27,0.1191},
    {287,1.0},
    {31,0.1191},
    {38,1.0},
    {276,0.1191},
    {99,0.0002014},
    {348,0.0002014},
    {329,0.1191},
    {80,0.1191},
    {26,0.1191},
    {444,1.0},
    {221,1.0},
    {275,0.1191},
    {404,1.0},
    {372,1.404},
    {417,2.127},
    {190,2.127},
    {123,1.404},
    {111,1088.0},
    {427,47.13},
    {207,1.0},
    {83,1.0},
    {107,1547000.0},
    {269,1.0},
    {237,1.0},
    {20,1.0},
    {356,347700},
    {332,1.0},
    {199,47.13},
    {388,1.0},
    {101,1.0},
    {350,1.0},
    {146,14.5805},
    {300,1.0},
    {397,1.0},
    {250,0.01538},
    {1,0.01538},
    {164,1.0},
    {51,1.0},
    {153,101.8},
    {297,1.0},
    {48,1.0},
    {0,0.02064},
    {93,1.0},
    {145,59.1811},
    {342,1.0},
    {175,22.81},
    {191,1.0},
    {418,1.0},
    {213,1.0},
    {440,1.0},
    {211,1.0},
    {390,330.4},
    {381,1627.45},
    {134,1627.45},
    {152,330.4},
    {437,1.0},
    {25,1.0},
    {391,101.8},
    {90,0.1191},
    {339,0.1191},
    {147,6.65621},
    {274,1.0},
    {10,0.002014},
    {33,1.0},
    {282,1.0},
    {259,0.002014},
    {188,3.879},
    {415,3.879},
    {47,1.0},
    {296,1.0},
    {434,1.0},
    {166,0.5104},
    {61,1.0},
    {408,45.62},
    {256,0.002014},
    {7,0.002014},
    {179,22.81},
    {15,1.0},
    {243,1.0},
    {310,1.0},
    {208,1.0},
    {376,1.0},
    {68,0.1191},
    {317,0.1191},
    {127,1.0},
    {405,11.08},
    {326,0.1191},
    {77,0.1191},
    {230,1.0},
    {331,0.1191},
    {209,1.0},
    {435,1.0},
    {82,0.1191},
    {84,1.0},
    {30,1.0},
    {279,1.0},
    {333,1.0},
    {39,1.0},
    {202,1.0},
    {431,1.0},
    {337,1.0},
    {425,0.2147},
    {198,0.2147},
    {88,1.0},
    {253,0.01036},
    {143,1627.45},
    {173,22.81},
    {4,0.01036},
    {235,1.0},
    {283,1.0},
    {34,1.0},
    {288,1.0},
    {449,1.0},
    {384,14.5805},
    {443,1.0},
    {220,1.0},
    {137,14.5805},
    {160,54.8},
    {438,1.0},
    {156,1.0},
    {65,1.0},
    {314,1.0},
    {212,1.0},
    {58,0.1191},
    {307,0.1191},
    {103,1.0},
    {334,0.1191},
    {85,0.1191},
    {377,365.3},
    {352,1.0},
    {176,22.81},
    {139,1.60809},
    {386,1.60809},
    {407,45.62},
    {233,1.0},
    {216,1.0},
    {81,0.1191},
    {318,1.0},
    {40,0.1191},
    {392,54.8},
    {204,1.0},
    {362,25.24},
    {113,20.23},
    {433,1.0},
    {154,54.8},
    {289,0.1191},
    {292,1.0},
    {186,30.52},
    {380,1637.13},
    {43,1.0},
    {210,1.0},
    {41,0.1191},
    {290,0.1191},
    {387,0.038592},
    {280,0.1191},
    {271,0.1191},
    {22,0.1191},
    {394,1.0},
    {140,0.0389136},
    {161,1.0},
    {246,1.0},
    {389,1012.0},
    {151,1012.0},
    {46,0.1191},
    {6,1.0},
    {224,1.0},
    {16,1.0},
    {379,1.0},
    {248,1.0},
    {132,1.0},
    {265,1.0},
    {255,1.0},
    {436,1.0},
    {368,3.185},
    {67,0.1191},
    {316,0.1191},
    {119,3.185},
    {328,1.0},
    {360,1207},
    {79,1.0},
    {60,1.0},
    {309,1.0},
    {200,1.0},
    {406,45.62},
    {428,1.0},
    {294,0.1191},
    {45,0.1191},
    {178,22.81},
    {128,365.3},
    {195,0.008039},
    {340,0.1191},
    {104,1.0},
    {229,1.0},
    {353,1.0},
    {91,0.1191},
    {422,0.008039},
    {346,1.0},
    {214,1.0},
    {441,1.0},
    {97,1.0},
    {226,1.0},
    {142,1637.13},
    {18,0.1191},
    {359,6831},
    {251,1.0},
    {131,1.0},
    {110,6334.0},
    {366,0.008348},
    {144,435.237},
    {117,0.008348},
    {304,0.1191},
    {373,1.404},
    {124,1.404},
    {102,1.0},
    {239,1.0},
    {430,1.0},
    {155,81880.0},
    {57,1.0},
    {225,1.0},
    {112,99.11},
    {2,1.0},
    {116,0.0004536},
    {365,0.0004536},
    {361,119.9},
    {168,0.4316},
    {108,322600.0},
    {357,347700},
    {401,0.4316},
    {70,1.0},
    {181,22.81},
    {343,0.1191},
    {245,1.0},
    {94,0.1191},
    {223,1.0},
    {382,435.237},
    {322,0.1191},
    {364,0.5644},
    {73,0.1191},
    {135,435.237},
    {8,0.002014},
    {241,1.0},
    {330,0.1191},
    {249,0.02064},
    {205,1.0},
    {182,34.91},
    {281,0.1191},
    {170,0.07358},
    {403,0.07358},
    {32,0.1191},
    {409,34.91},
    {159,101.8},
    {395,0.1097},
    {121,14.93},
    {370,14.93},
    {162,0.1097},
    {56,1.0},
    {305,1.0},
    {260,0.002014},
    {345,1.0},
    {96,1.0},
    {242,1.0},
    {257,0.002014},
    {100,1.0},
    {349,1.0},
    {319,1.0},
    {270,1.0},
    {78,1.0},
    {327,1.0},
    {165,4.078},
    {398,4.078},
    {64,0.1191},
    });
  
unordered_map<int,float> Nevents({
    {268,1},
    {19,1},
    {432,1},
    {203,1},
    {5,1977800},
    {315,1},
    {66,1},
    {254,1991000},
    {423,250000},
    {196,250000},
    {262,1},
    {13,1},
    {163,700000},
    {92,1},
    {122,408000},
    {371,500000},
    {341,1},
    {396,678600},
    {335,1},
    {240,1},
    {447,1},
    {86,1},
    {299,1},
    {52,1},
    {301,1},
    {50,1},
    {378,128640000},
    {129,41729120},
    {24,1},
    {338,1},
    {149,21495421},
    {89,1},
    {400,13280000},
    {167,7932650},
    {174,200000},
    {219,1},
    {141,1},
    {14,1},
    {263,1},
    {63,1},
    {312,1},
    {148,20258624},
    {369,492000},
    {238,1},
    {325,1},
    {303,1},
    {54,1},
    {76,1},
    {120,55000},
    {12,490000},
    {234,1},
    {448,1},
    {261,489000},
    {115,8744768},
    {306,1},
    {11,496000},
    {278,1},
    {29,1},
    {118,500000},
    {217,1},
    {36,1},
    {285,1},
    {21,1},
    {411,984000},
    {184,992000},
    {157,42331295},
    {414,1000000},
    {187,500000},
    {323,1},
    {74,1},
    {402,4911941},
    {28,1},
    {171,1},
    {277,1},
    {272,1},
    {320,1},
    {419,1},
    {169,4994543},
    {192,1},
    {158,10037851},
    {23,1},
    {424,250000},
    {354,1},
    {105,1},
    {197,250000},
    {344,1},
    {416,7525991},
    {189,9999987},
    {95,1},
    {59,1},
    {446,1},
    {258,1},
    {9,1},
    {232,1},
    {308,1},
    {42,1},
    {267,1},
    {150,1},
    {291,1},
    {336,1},
    {351,1},
    {87,1},
    {53,1},
    {374,482400},
    {49,1},
    {298,1},
    {302,1},
    {218,1},
    {37,1},
    {172,2000000},
    {439,1},
    {295,1},
    {69,1},
    {55,1},
    {71,1},
    {215,1},
    {227,1},
    {442,1},
    {244,1},
    {355,1},
    {286,1},
    {3,1937000},
    {252,1768000},
    {106,9403812},
    {412,958000},
    {130,9000000},
    {426,64310000},
    {185,988000},
    {114,1},
    {133,22226705},
    {363,1},
    {284,1},
    {313,1},
    {222,1},
    {421,240000},
    {194,232300},
    {445,1},
    {177,200000},
    {383,5932701},
    {35,1},
    {228,1},
    {413,12575000},
    {273,1},
    {367,485000},
    {385,8402687},
    {138,20432728},
    {136,14313274},
    {72,1},
    {311,1},
    {17,989200},
    {266,1000000},
    {62,1},
    {321,1},
    {347,1},
    {201,1},
    {180,200000},
    {358,1},
    {109,1},
    {429,1},
    {98,1},
    {231,1},
    {236,1},
    {264,1},
    {450,1},
    {393,1},
    {293,1},
    {410,1000000},
    {44,1},
    {247,1},
    {193,968000},
    {126,827046},
    {324,1},
    {183,953600},
    {125,485420},
    {75,1},
    {375,1021031},
    {420,871500},
    {206,1},
    {399,8822000},
    {27,1},
    {287,1},
    {31,1},
    {38,1},
    {276,1},
    {99,100000},
    {348,100000},
    {329,1},
    {80,1},
    {26,1},
    {444,1},
    {221,1},
    {275,1},
    {404,1},
    {372,913200},
    {417,1102578},
    {190,918508},
    {123,936000},
    {111,1261997},
    {427,3885000},
    {207,1},
    {83,1},
    {107,8611681},
    {269,1},
    {237,1},
    {20,1},
    {356,1},
    {332,1},
    {199,3807850},
    {388,1},
    {101,1},
    {350,1},
    {146,21709087},
    {300,1},
    {397,1},
    {250,150000},
    {1,150000},
    {164,1},
    {51,1},
    {153,5748466},
    {297,1},
    {48,1},
    {0,145800},
    {93,1},
    {145,14313274},
    {342,1},
    {175,200000},
    {191,1},
    {418,1},
    {213,1},
    {440,1},
    {211,1},
    {390,20456037},
    {381,29425374},
    {134,35862893},
    {152,10037851},
    {437,1},
    {25,1},
    {391,5652357},
    {90,1},
    {339,1},
    {147,20432728},
    {274,1},
    {10,484000},
    {33,1},
    {282,1},
    {259,482000},
    {188,9725000},
    {415,1},
    {47,1},
    {296,1},
    {434,1},
    {166,8940000},
    {61,1},
    {408,200000},
    {256,488000},
    {7,490000},
    {179,200000},
    {15,1},
    {243,1},
    {310,1},
    {208,1},
    {376,1},
    {68,1},
    {317,1},
    {127,1},
    {405,7758900},
    {326,1},
    {77,1},
    {230,1},
    {331,1},
    {209,1},
    {435,1},
    {82,1},
    {84,1},
    {30,1},
    {279,1},
    {333,1},
    {39,1},
    {202,1},
    {431,1},
    {337,1},
    {425,1154000},
    {198,1000000},
    {88,1},
    {253,1923000},
    {143,35862893},
    {173,200000},
    {4,1976900},
    {235,1},
    {283,1},
    {34,1},
    {288,1},
    {449,1},
    {384,19771294},
    {443,1},
    {220,1},
    {137,21709087},
    {160,4328648},
    {438,1},
    {156,1},
    {65,1},
    {314,1},
    {212,1},
    {58,1},
    {307,1},
    {103,1},
    {334,1},
    {85,1},
    {377,100790000},
    {352,1},
    {176,200000},
    {139,20258624},
    {386,7633949},
    {407,200000},
    {233,1},
    {216,1},
    {81,1},
    {318,1},
    {40,1},
    {392,2809978},
    {204,1},
    {362,1},
    {113,214088},
    {433,1},
    {154,4328648},
    {289,1},
    {292,1},
    {186,2921180},
    {380,28084244},
    {43,1},
    {210,1},
    {41,1},
    {290,1},
    {387,3273980},
    {280,1},
    {271,1},
    {22,1},
    {394,1},
    {140,21495421},
    {161,1},
    {246,1},
    {389,68898175},
    {151,42331295},
    {46,1},
    {6,1},
    {224,1},
    {16,1},
    {379,1},
    {248,1},
    {132,1},
    {265,1},
    {255,1},
    {436,1},
    {368,494000},
    {67,1},
    {316,1},
    {119,1472800},
    {328,1},
    {360,1},
    {79,1},
    {60,1},
    {309,1},
    {200,1},
    {406,200000},
    {428,1},
    {294,1},
    {45,1},
    {178,200000},
    {128,43506449},
    {195,1944000},
    {340,1},
    {104,1},
    {229,1},
    {353,1},
    {91,1},
    {422,1986000},
    {346,1},
    {214,1},
    {441,1},
    {97,1},
    {226,1},
    {142,22226705},
    {18,1},
    {359,1},
    {251,1},
    {131,1},
    {110,1877745},
    {366,3774800},
    {144,21250517},
    {117,3949400},
    {304,1},
    {373,911500},
    {124,923500},
    {102,1},
    {239,1},
    {430,1},
    {155,9586029},
    {57,1},
    {225,1},
    {112,224332},
    {2,1},
    {116,3858600},
    {365,3908600},
    {361,1},
    {168,811306},
    {108,5529691},
    {357,1},
    {401,835296},
    {70,1},
    {181,200000},
    {343,1},
    {245,1},
    {94,1},
    {223,1},
    {382,25468933},
    {322,1},
    {364,8382600},
    {73,1},
    {135,21250517},
    {8,498000},
    {241,1},
    {330,1},
    {249,149400},
    {205,1},
    {182,7794186},
    {281,1},
    {170,13276146},
    {403,13736000},
    {32,1},
    {409,9598000},
    {159,5748466},
    {395,592000},
    {121,493000},
    {370,491200},
    {162,499800},
    {56,1},
    {305,1},
    {260,488000},
    {345,1},
    {96,1},
    {242,1},
    {257,476000},
    {100,1},
    {349,1},
    {319,1},
    {270,1},
    {78,1},
    {327,1},
    {165,8729288},
    {398,4691915},
    {64,1},
});
*/

unordered_map<int,float> xsecs({
{130,88.287},
{193,0.1703},
{194,0.2086},
{195,0.008039},
{196,0.05565},
{197,0.013989},
{198,0.2147},
{165,4.078},
{166,0.5104},
{167,0.2432},
{168,0.4316},
{169,0.2149},
{170,0.07358},
{115,0.5644},
{116,0.0004536},
{117,0.008348},
{118,14.93},
{119,3.185},
{120,3.185},
{121,14.93},
{122,3.185},
{123,1.404},
{124,1.404},
{125,1.404},
{172,11.08},
{173,7.09644444444},
{174,7.09644444444},
{175,7.09644444444},
{176,7.09644444444},
{177,7.09644444444},
{178,7.09644444444},
{179,7.09644444444},
{180,7.09644444444},
{181,7.09644444444},
{182,34.91},
{183,1.0315},
{185,0.0118},
{186,2.7757},
{187,0.0896},
{188,0.237},
{189,0.212},
{190,0.952},
{151,877.8},
{152,304.4},
{153,111.5},
{154,44.03},
{155,81880.0},
{3,0.002014},
{4,0.01036},
{5,0.01595},
{162,0.1097},
{163,0.5439},
{199,47.13},
{1,0.01538},   
    });
   
unordered_map<int,float> Nevents({
{193,968000},
{194,232300},
{195,1944000},
{196,250000},
{197,250000},
{198,1000000},
{165,4958713},
{166,8098000},
{167,6730495},
{168,811306},
{169,1771042},
{170,2471659},
{1,128500},
{162,499800},
{163,700000},
{115,7547567},
{116,3256700},
{117,3495100},
{118,500000},
{119,705800},
{120,55000},
{121,493000},
{122,408000},
{123,353600},
{124,923500},
{125,485420},
{199,3807850},
{3,1032500},
{4,1976900},
{5,1977800},
{151,27229712},
{152,8723495},
{153,172167},
{154,2461611},
{155,9586029},
{172,2000000},
{173,200000},
{174,200000},
{175,200000},
{176,200000},
{177,200000},
{178,200000},
{179,200000},
{180,200000},
{181,200000},
{182,7794186},
{183,953600},
{185,988000},
{186,2921180},
{187,500000},
{188,9725000},
{189,9999987},
{190,918508},
{130,1644828},
});

bool isMC(int SampleFlag){
    bool is_in = data_flags.find(SampleFlag) != data_flags.end();
    if (is_in == true) return false;
    else return true;
}

float getLumi(int Year, bool IsMC){
    if (IsMC == false) return 1.;
    else if(Year == 2016) return -999;
    //else if (Year == 2017) return 41.48;
    else if (Year == 2017) return 41.54;
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
TMVA::Experimental::RBDT<> bdt("xgb_SM_v100", "https://ttedesch.web.cern.ch/ttedesch/VBS_ML_v4/optimized_model_SM_clean.root");
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