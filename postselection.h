/// Header file with functions needed to execute the Python version of
/// postselection step of the analysis. The header is declared to the
/// ROOT C++ interpreter prior to the start of the analysis via the
/// `ROOT.gInterpreter.Declare()` function.
///
 
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "Math/Vector4D.h"
#include "TStyle.h"
#include <map>

using namespace ROOT::VecOps;
using RNode = ROOT::RDF::RNode;
using rvec_f = const RVec<float> &;
using rvec_i = const RVec<int> &;


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

/*
//RVec<int> SelectLepton(rvec_f jet_eta, rvec_f jet_phi, rvec_f lepton_pt, rvec_f lepton_eta, rvec_f lepton_phi, rvec_i VBSJets_idx){
RVec<int> SelectLepton(rvec_f lepton_pt, rvec_f lepton_eta, rvec_f lepton_phi, rvec_i lepton_pdgId, rvec_i lepton_tightId, rvec_i lepton_looseId,rvec_f lepton_jetRelIso, rvec_f Iso04_all, rvec_f lepton_mvaFall17V2Iso_WPL, rvec_f lepton_mvaFall17V2Iso_WP90, rvec_f jet_eta, rvec_f jet_phi, rvec_i VBSJets_idx){
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
        int lepflav = abs(lepton_pdgId[i]);

        if (lepflav == 11){
            IsLooseID = (bool)lepton_mvaFall17V2Iso_WPL[i];
            IsTightID = (bool)lepton_mvaFall17V2Iso_WP90[i];
            IsTightIso = lepton_jetRelIso[i]<ISO_CUT_ELE && lepton_jetRelIso[i]>=0.;
            IsLooseIso = lepton_jetRelIso[i]<1. && lepton_jetRelIso[i]>=0.;
            IsInEtaRegion = abs(lepton_eta[i])<ETA_CUT_ELE && !(abs(lepton_eta[i])>1.4442 && abs(lepton_eta[i])<1.566);
            IsInPtRegion = lepton_pt[i] > PT_CUT_ELE;
        }
        
        else if (lepflav == 13){
            IsTightID = lepton_tightId[i];
            IsLooseID = lepton_looseId[i];
            IsTightIso = Iso04_all[i]<ISO_CUT_MU && Iso04_all[i]>=0.;
            IsLooseIso = Iso04_all[i]<1. && Iso04_all[i]>=0.;
            IsInEtaRegion = abs(lepton_eta[i]) < ETA_CUT_MU;
            IsInPtRegion = lepton_pt[i] > PT_CUT_MU;
        }

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
}
*/

RVec<int> SelectElectron(rvec_f lepton_pt, rvec_f lepton_eta, rvec_f lepton_phi, rvec_f lepton_jetRelIso, rvec_f lepton_mvaFall17V2Iso_WPL, rvec_f lepton_mvaFall17V2Iso_WP90, rvec_f jet_eta, rvec_f jet_phi, rvec_i VBSJets_idx){
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
        IsLooseID = (bool)lepton_mvaFall17V2Iso_WPL[i];
        IsTightID = (bool)lepton_mvaFall17V2Iso_WP90[i];
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

RVec<int> SelectMuon(rvec_f lepton_pt, rvec_f lepton_eta, rvec_f lepton_phi, rvec_i lepton_tightId, rvec_i lepton_looseId, rvec_f lepton_jetRelIso, rvec_f Iso04_all, rvec_f jet_eta, rvec_f jet_phi, rvec_i VBSJets_idx){
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

int DetermineGoodLepton(bool HLT_IsoMu27, bool HLT_Mu50, bool HLT_Ele35_WPTight_Gsf, bool HLT_Ele32_WPTight_Gsf_L1DoubleEG, bool HLT_Photon200, bool HLT_PFHT250, bool HLT_PFHT350, rvec_i Electron_idx, rvec_f Electron_pt, rvec_f Electron_eta, rvec_f Electron_mvaFall17V2Iso_WPL, rvec_f Electron_jetRelIso, rvec_i Muon_idx, rvec_f Muon_pt, rvec_f Muon_eta, rvec_f Muon_pfRelIso04_all){
    bool passMu = false;
    bool passEle = false;
    bool passHT = false;
    
    if(HLT_IsoMu27 || HLT_Mu50) passMu = true;
    if(HLT_Ele35_WPTight_Gsf || HLT_Ele32_WPTight_Gsf_L1DoubleEG || HLT_Photon200) passEle = true;
    if(HLT_PFHT250 || HLT_PFHT350) passHT = true;
    
    if(Electron_idx[1] != -1) bool ele_lepton_veto = LepVetoEle(Electron_idx, Electron_pt, Electron_eta, Electron_mvaFall17V2Iso_WPL, jetRelIso, Muon_pt, Muon_eta, Muon_pfRelIso04_all);
    if(Muon_idx[1] != -1)bool ele_lepton_veto = LepVetoMu(Muon_idx, Electron_pt, Electron_eta, Electron_mvaFall17V2Iso_WPL, jetRelIso, Muon_pt, Muon_eta, Muon_pfRelIso04_all);

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
                if(Electron_pt[Electron_idx[0]] > Muon_pt[Muon_idx[0]])
                    GoodLeptonFamily = 0;
                    //lepton_TightRegion[0] = copy.deepcopy(ele_TightRegion)
                    SingleEle = true
                    SingleMu = false
                else:
                    GoodLeptonFamily = 1;
                    //lepton_TightRegion[0] = copy.deepcopy(mu_TightRegion)
                    SingleEle = false
                    SingleMu = true
            }
            else GoodLeptonFamily = -1;
        }
        else GoodLeptonFamily = -1;
    }
    if(!(SingleEle || SingleMu)) GoodLeptonFamily = -1;
    
    return GoodLeptonFamily;
}

RVec<int> SelectAndVetoTaus(rvec_f Tau_pt, rvec_f Tau_eta, rvec_f Tau_phi, rvec_f Tau_idDeepTau2017v2p1VSjet, rvec_f Tau_idDeepTau2017v2p1VSe, rvec_f Tau_idDeepTau2017v2p1VSmu, rvec_f Tau_idDecayModeNewDMs, rvec_f Lepton_eta, rvec_f Lepton_phi, rvec_i Lepton_pdgId, rvec_i Lepton_idx, rvec_f Jet_eta, rvec_f Jet_phi, rvec_i Jet_idx)
{
    //setting jet-related quantities if isolation from them is needed
    float jet1eta = Jet_eta[Jet_idx[0]];
    float jet2eta = Jet_eta[Jet_idx[1]];
    float jet1phi = Jet_phi[Jet_idx[0]];
    float jet2phi = Jet_phi[Jet_idx[1]];
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
        if (abs(Lepton_pdgId[Lepton_idx[0]])==11) cutloose_vsjet = ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_ELE;
        else if (abs(Lepton_pdgId[Lepton_idx[0]])==13) cutloose_vsjet = ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_MU;

        if ((Tau_idDeepTau2017v2p1VSjet[i]>=cutloose_vsjet && Tau_idDeepTau2017v2p1VSe[i]>=ID_TAU_RECO_DEEPTAU_VSELE && Tau_idDeepTau2017v2p1VSmu[i]>=ID_TAU_RECO_DEEPTAU_VSMU && Tau_idDecayModeNewDMs[i]) && deltaR(Tau_eta[i], Tau_phi[i], Lepton_eta[Lepton_idx[0]], Lepton_phi[Lepton_idx[0]])>DR_OVERLAP_CONE_TAU && deltaR(Tau_eta[i], Tau_phi[i], jet1eta, jet1phi)>isocone && deltaR(Tau_eta[i], Tau_phi[i], jet2eta, jet2phi)>isocone && Tau_pt[i]>=PT_CUT_TAU && abs(Tau_eta[i])<=ETA_CUT_TAU){
            
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
    if(nTau!=1) idx[1] = -1;                                                                                               
    
    return idx;
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

bool LepVetoEle(rvec_i Electron_idx, rvec_f Electron_pt, rvec_f Electron_eta, rvec_f Electron_mvaFall17V2Iso_WPL, rvec_f jetRelIso, rvec_f Muon_pt, rvec_f Muon_eta, rvec_f Muon_pfRelIso04_all)
{
    bool IsEleVetoPassed = true;
    bool IsMuVetoPassed = true;
    for (size_t i = 0; i < Electron_pt.size(); i++) {
        if(i != Electron_idx[0] && Electron_mvaFall17V2Iso_WPL && Electron_pt > PT_CUT_LEP_VETO_ELE && abs(Electron_eta) < ETA_CUT_LEP_VETO_ELE && !(abs(Electron_eta)>1.4442 && abs(Electron_eta)<1.566) && Electron_jetRelIso < REL_ISO_CUT_LEP_VETO_ELE) IsEleVetoPassed = false;
    }
    if(IsEleVetoPassed == true){
        for (size_t i = 0; i < Muon_pt.size(); i++) {
            if(Muon_looseId && Muon_pt > PT_CUT_LEP_VETO_MU && abs(Muon_eta) < ETA_CUT_LEP_VETO_MU && Muon_pfRelIso04_all < REL_ISO_CUT_LEP_VETO_MU) IsMuVetoPassed = false;
        }
    }
    return IsEleVetoPassed and IsMuVetoPassed;
}

bool LepVetoMu(rvec_i Muon_idx, rvec_f Electron_pt, rvec_f Electron_eta, rvec_f Electron_mvaFall17V2Iso_WPL, rvec_f jetRelIso, rvec_f Muon_pt, rvec_f Muon_eta, rvec_f Muon_pfRelIso04_all)
{
    bool IsEleVetoPassed = true;
    bool IsMuVetoPassed = true;
    for (size_t i = 0; i < Electron_pt.size(); i++) {
        if(Electron_mvaFall17V2Iso_WPL && Electron_pt > PT_CUT_LEP_VETO_ELE && abs(Electron_eta) < ETA_CUT_LEP_VETO_ELE && !(abs(Electron_eta)>1.4442 && abs(Electron_eta)<1.566) && Electron_jetRelIso < REL_ISO_CUT_LEP_VETO_ELE) IsEleVetoPassed = false;
    }
    if(IsEleVetoPassed == true){
        for (size_t i = 0; i < Muon_pt.size(); i++) {
            if(i != Muon_idx[0] && Muon_looseId && Muon_pt > PT_CUT_LEP_VETO_MU && abs(Muon_eta) < ETA_CUT_LEP_VETO_MU && Muon_pfRelIso04_all < REL_ISO_CUT_LEP_VETO_MU) IsMuVetoPassed = false;
        }
    }
    return IsEleVetoPassed and IsMuVetoPassed;
}

