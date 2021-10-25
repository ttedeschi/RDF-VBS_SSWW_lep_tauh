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

const size_t PT_CUT_ELE=  35;
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


bool LepVetoEle(rvec_i Electron_idx, rvec_f Electron_pt, rvec_f Electron_eta, rvec_b Iso_WPL, rvec_i jetRelIso, rvec_f Muon_pt, rvec_f Muon_eta, rvec_i Iso04_all, rvec_b Muon_looseId )
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



bool LepVetoMu(rvec_i Muon_idx, rvec_f Electron_pt, rvec_f Electron_eta, rvec_b Electron_mvaFall17V2Iso_WPL, rvec_i Electron_jetRelIso, rvec_f Muon_pt, rvec_f Muon_eta, rvec_i Muon_pfRelIso04_all, rvec_b Muon_looseId)
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
    auto best_mass = -1;
    size_t best_i1 = 0; size_t best_i2 = 0;
    for (size_t i = 0; i < idx_cmb[0].size(); i++) {
        const auto i1 = idx_cmb[0][i];
        const auto i2 = idx_cmb[1][i];
        //std::cout<<i1<<i2<<endl;
        if (abs(eta[GoodJets_idx[i1]] - eta[GoodJets_idx[i2]]) >= DELTAETA_JJ_CUT) {
            ROOT::Math::PtEtaPhiMVector p1(pt[GoodJets_idx[i1]], eta[GoodJets_idx[i1]], phi[GoodJets_idx[i1]], mass[GoodJets_idx[i1]]);
            ROOT::Math::PtEtaPhiMVector p2(pt[GoodJets_idx[i2]], eta[GoodJets_idx[i2]], phi[GoodJets_idx[i2]], mass[GoodJets_idx[i2]]);
            const auto this_mass = (p1 + p2).M();
            if (this_mass > best_mass) {
                best_mass = this_mass;
                best_i1 = GoodJets_idx[i1];
                best_i2 = GoodJets_idx[i2];
            }
        }
    } 
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

float GetLeptonTightFlag(rvec_i Electron_idx, rvec_i Muon_idx, int GoodLeptonFamily){
    if (GoodLeptonFamily == 0) return Electron_idx[1];
    else return Muon_idx[1];
}

float GetTau(rvec_f pt, rvec_i idx){
    return pt[idx[0]];
}

RVec<int> SelectElectron(rvec_f lepton_pt, rvec_f lepton_eta, rvec_f lepton_phi, rvec_f lepton_jetRelIso, rvec_b lepton_mvaFall17V2Iso_WPL, rvec_f lepton_mvaFall17V2Iso_WP90, rvec_f jet_eta, rvec_f jet_phi, rvec_i VBSJets_idx){
    //setting jet-related quantities if isolation from them is needed
    const auto jet1_idx = VBSJets_idx[0];
    const auto jet2_idx = VBSJets_idx[1];

    const auto jet1eta = jet_eta[jet1_idx];
    const auto jet2eta = jet_eta[jet2_idx];
    const auto jet1phi = jet_phi[jet1_idx];
    const auto jet2phi = jet_phi[jet2_idx];
    
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
    const auto jet1_idx = VBSJets_idx[0];
    const auto jet2_idx = VBSJets_idx[1];
    //const auto jet1eta = lepton_pt[0]
    const auto jet1eta = jet_eta[jet1_idx];
    const auto jet2eta = jet_eta[jet2_idx];
    const auto jet1phi = jet_phi[jet1_idx];
    const auto jet2phi = jet_phi[jet2_idx];
    
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
        idx[1] = 0;
    }
    else if (LnTleptons_idx.size() > 0){
        idx[0] = LnTleptons_idx[0];
        idx[1] = 1;
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
            else if(not ele_lepton_veto and mu_lepton_veto){           
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
                    idx[1] = 0;
                } 
                else{
                    idx[0] = i;
                    idx[1] = 1;
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
        if ((Jet_btagDeepFlavB[i]>=BTAG_WP_VALUE) && (Jet_pt[i]>BTAG_PT_CUT) && (abs(Jet_eta[i])<BTAG_ETA_CUT)) return true;
    }
    return false;
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
    
    TFile inFile("FR_vsjet2_vsmuT_ZZ.root"); 
    TH2F *histo;
    if (abs(pdgId) == 11) histo = (TH2F*)inFile.Get("hFRDataeledif");
    else if (abs(pdgId) == 13) histo = (TH2F*)inFile.Get("hFRDatamudif");

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
    
    TFile inFile("FR_vsjet4_vsmuT_ZZ.root"); 
    TH2F *histo;
    if (abs(pdgId) == 11) histo = (TH2F*)inFile.Get("hFRDataeledif");
    else if (abs(pdgId) == 13) histo = (TH2F*)inFile.Get("hFRDatamudif");

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
    TFile inFile("FR_vsjet2_vsmuT_ZZ.root"); 

    TH2F *histo = (TH2F*)inFile.Get("hFRDatataudif");

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
        
float SFFakeRatio_tau_calc_vsjet4(float pT, float eta){

    TFile inFile("FR_vsjet4_vsmuT_ZZ.root"); 

    TH2F *histo = (TH2F*)inFile.Get("hFRDatataudif");

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

float GetEventSFFake(float lepton_SFFake, float tau_SFFake, int lepton_LnTRegion, int tau_LnTRegion){
    if(lepton_LnTRegion==1 && tau_LnTRegion==0) return lepton_SFFake;
    else if (lepton_LnTRegion==0 && tau_LnTRegion==1) return tau_SFFake;
    else if (lepton_LnTRegion==1 && tau_LnTRegion==1) return lepton_SFFake*tau_SFFake;
    else if (lepton_LnTRegion==0 && tau_LnTRegion==0) return 0.;
}

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

//RVec<int> SelectVBSQGenJet(rvec_i GenJet_partonFlavour){
//    RVec<int> idx;
//    for (int i = 0; i < GenJet_partonFlavour.size(); i++) {
//        if(abs(GenJet_partonFlavour[i])>0 && abs(GenJet_partonFlavour[i])<10) idx.emplace_back(i);
//    }
//    return idx;
//}

bool atleast2Jets(rvec_i GenJet_idx){
    if(GenJet_idx.size() < 2) return false;
    else return true;
}

int IsGenMatched(int i, int j){
    if(i == j) return 1;
    else return 0;
}

/*
float Lepton_IDIso_SF(float Lepton_pt, float Lepton_eta, int Lepton_pdgId){
    if(abs(Lepton_pdgId)==13){
        f_muon_id = ROOT.TFile.Open("Muon_RunBCDEF_SF_ID_2017.root", "READ")
        f_muon_iso = ROOT.TFile.Open("Muon_RunBCDEF_SF_ISO_2017.root", "READ")
        h_muon_id = ROOT.TH2D(f_muon_id.Get("NUM_TightID_DEN_genTracks_pt_abseta"))
        h_muon_iso = ROOT.TH2D(f_muon_iso.Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta"))
      
        float pt = Lepton_pt;
        float abseta = abs(Lepton_eta);
        int binx_id = h_muon_id.GetXaxis().FindBin(pt);
        int biny_id = h_muon_id.GetYaxis().FindBin(abseta);
        int nxbins_id = h_muon_id.GetXaxis().GetNbins();
        int nybins_id = h_muon_id.GetYaxis().GetNbins();

        if (binx_id > nxbins_id) binx_id = copy.deepcopy(nxbins_id);
        else if (binx_id <= 0) binx_id = 1;

        if (biny_id > nybins_id) biny_id = copy.deepcopy(nybins_id);
        else if (biny_id <= 0) biny_id = 1;

        binx_iso = h_muon_iso.GetXaxis().FindBin(pt)
        biny_iso = h_muon_iso.GetYaxis().FindBin(abseta)
        nxbins_iso = h_muon_iso.GetXaxis().GetNbins()
        nybins_iso = h_muon_iso.GetYaxis().GetNbins()
        
        if (binx_iso > nxbins_iso) binx_iso = copy.deepcopy(nxbins_iso);
        else if(binx_iso <= 0) binx_iso = 1;
        if (biny_iso > nybins_iso) biny_iso = copy.deepcopy(nybins_iso);
        else if (biny_iso <= 0) biny_iso = 1;
        
        SF_muon_id = copy.deepcopy(h_muon_id.GetBinContent(binx_id, biny_id));
        SF_muon_iso = copy.deepcopy(h_muon_iso.GetBinContent(binx_iso, biny_iso));

        f_muon_id.Close();
        f_muon_iso.Close();

        return SF_muon_id*SF_muon_iso
    }
    else if (abs(Lepton_pdgId)==11){
        TFile f_electron_id("Electron_MVA90_2017.root", "READ");
        TH2F *h_electron_id = (TH2F*)f_electron_id.Get("EGamma_SF2D");

        if(Lepton_pt < 20.) TFile f_electron_reco("EGM2D_2017_passingRECO_lowEt.root", "READ");
        else if (lepton.pt >= 20.) TFile f_electron_reco("EGM2D_2017_passingRECO_highEt.root", "READ");
        TH2F *h_electron_reco = (TH2F*)f_electron_reco.Get("EGamma_SF2D");

        auto pt = Lepton_pt;
        auto eta = Lepton_eta;
        auto biny_id = h_electron_id->GetXaxis()->FindBin(pt);
        auto binx_id = h_electron_id->GetYaxis()->FindBin(eta);
        auto nxbins_id = h_electron_id->GetXaxis()->GetNbins();
        auto nybins_id = h_electron_id->GetYaxis()->GetNbins();
        if (binx_id > nxbins_id) binx_id = nxbins_id;
        else if (binx_id <= 0) binx_id = 1;
        if (biny_id > nybins_id) biny_id = nybins_id;
        else if (biny_id <= 0) biny_id = 1;

        auto biny_reco = h_electron_reco->GetXaxis()->FindBin(pt);
        auto binx_reco = h_electron_reco->GetYaxis()->FindBin(eta);
        auto nxbins_reco = h_electron_reco->GetXaxis()->GetNbins();
        auto nybins_reco = h_electron_reco->GetYaxis()->GetNbins();
        if (binx_reco > nxbins_reco) binx_reco = *nxbins_reco;
        else if (binx_reco <= 0) binx_reco = 1;
        if (biny_reco > nybins_reco) biny_reco = *nybins_reco;
        else if (biny_reco <= 0) biny_reco = 1;

        auto SF_electron_id = h_electron_id.GetBinContent(binx_id, biny_id);
        auto SF_electron_reco = h_electron_reco.GetBinContent(binx_reco, biny_reco);

        f_electron_id.Close();
        f_electron_reco.Close();

        return SF_electron_id*SF_electron_reco;
    }

    else return -1.;
}
*/

RVec<RVec<float>> getTauSF(float SelectedTau_pt, float SelectedTau_eta, int SelectedTau_genPartFlav){
    RVec<float> vsJet, vsEle, vsMu;
    string id;
    std::string year = std::to_string(2017);
    
    // vs Jet
    id =  "DeepTau2017v2p1VSjet";
    TString path = TString("tauSF/TauID_SF_pt_") + TString(id) + TString("_") + TString(year) + TString("ReReco") + TString(".root");
    TFile *f = new TFile(path);
    //TFile *f = new TFile();
    double_t pt = SelectedTau_pt;
    if (SelectedTau_genPartFlav==5){
        TString path_down = TString(TString(vsJetwp) + TString("_down"));
        TString path_cent = TString(TString(vsJetwp) + TString("_cent"));
        TString path_up = TString(TString(vsJetwp) + TString("_up"));
        TF1 * h_down = (TF1*)f->Get(path_down);
        TF1 * h_cent =  (TF1*)f->Get(path_cent);
        TF1 * h_up =  (TF1*)f->Get(path_up);
        vsJet.emplace_back(h_down->Eval(pt));
        vsJet.emplace_back(h_cent->Eval(pt));
        vsJet.emplace_back(h_up->Eval(pt));
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
    TString path_ele = TString("tauSF/TauID_SF_eta_") + TString(id) + TString("_") + TString(year) + TString("ReReco") + TString(".root");
    TFile *f_ele = new TFile(path_ele);
    TString histoname_ele = TString(vsElewp);
    TH1F * hist = (TH1F *) f_ele->Get(histoname_ele);
    if (SelectedTau_genPartFlav == 1 || SelectedTau_genPartFlav == 3){
        bin = hist->GetXaxis()->FindBin(eta);
        sf  = hist->GetBinContent(bin);
        err = hist->GetBinError(bin);
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
    TString path_mu = TString("tauSF/TauID_SF_eta_") + TString(id) + TString("_") + TString(year)  + TString("ReReco")+ TString(".root");
    TFile *f_mu = new TFile(path_mu);
    TString histoname_mu = TString(vsMuwp);
    TH1F * hist_mu = (TH1F *) f_mu->Get(histoname_mu);
    if (SelectedTau_genPartFlav == 2 || SelectedTau_genPartFlav == 4){
        bin = hist_mu->GetXaxis()->FindBin(eta);
        sf  = hist_mu->GetBinContent(bin);
        err = hist_mu->GetBinError(bin);
        vsMu.emplace_back(sf-err);
        vsMu.emplace_back(sf);
        vsMu.emplace_back(sf+err);
    }
    else{
        vsMu.emplace_back(1.0);
        vsMu.emplace_back(1.0);
        vsMu.emplace_back(1.0);
    }
    
    RVec<RVec<float>> result;
    result.emplace_back(vsJet);
    result.emplace_back(vsEle);
    result.emplace_back(vsMu);

    return result;
}

float efficiency(int flv, float eta, float pt, string year){
    TString path = TString("Btag_eff_") + TString(year) + TString(".root");
    TFile *infile = new TFile(path);
    TH2F * h;
    if(flv == 5){
        //TH2F * h = (TH2F *) infile->Get("h2_BTaggingEff_b")->CreateHistogram();
        TEfficiency *eff = (TEfficiency *) infile->Get("h2_BTaggingEff_b");
        h = (TH2F *) eff->CreateHistogram();
    }
    else if(flv == 4){
        TEfficiency *eff = (TEfficiency *) infile->Get("h2_BTaggingEff_c");
        h = (TH2F *) eff->CreateHistogram();
    }
    else{
        //h = (TH2F *) infile->Get("h2_BTaggingEff_udsg");
        TEfficiency *eff = (TEfficiency *) infile->Get("h2_BTaggingEff_udsg");
        h = (TH2F *) eff->CreateHistogram();
    }
    
    int binx = max(1, min(h->GetNbinsX(), h->GetXaxis()->FindBin(pt)));
    int biny = max(1, min(h->GetNbinsY(), h->GetYaxis()->FindBin(abs(eta))));
    
    return h->GetBinContent(binx,biny);
}
    
RVec<float> btagcalc(rvec_i GoodJets_idx, rvec_f Jet_pt, rvec_f Jet_eta, rvec_i Jet_partonFlavour, rvec_f Jet_btagDeepFlavB, rvec_f Jet_btagSF_deepjet_M_up, rvec_f Jet_btagSF_deepjet_M_down, rvec_f Jet_btagSF_deepjet_M, rvec_f Jet_btagDeepB){

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