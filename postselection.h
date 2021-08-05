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


const float thr = 50; 

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
const string BTAG_ALGO   =   'DeepFlv';
const string BTAG_WP     =   'M';
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


std::map<string, std::map<string, float> > WP_btagger = {
    { "CSVv2", {"L": 0.5803,"M": 0.8838,"T": 0.9693}},
    { "DeepCSV", {"L": 0.1522,"M": 0.4941,"T": 0.8001}},
    { "DeepFlv", {"L": 0.0521,"M": 0.3033,"T": 0.7489}},
};

RVec<size_t> GoodJets(rvec_i jetId, rvec_f eta, rvec_f pt, rvec_i puId)
{
   RVec<size_t> idx;
   for (size_t i = 0; i < pt.size(); i++) {
      if (jetId[i] >= 2 && abs(eta[i]) < 5. && pt[i] > thr && (pt[i] > 50. || (pt[i] <= 50. && puId[i] >= 7))) {
         idx.push_back(i);
         }
      }
   return idx;
}

RVec<size_t> SelectVBSJets_invmass(rvec_f pt, rvec_f eta, rvec_f phi, rvec_f mass, rvec_i jetId, rvec_i puId, rvec_i GoodJets_idx)
{
    RVec<size_t> idx(2);
    // Find first lepton pair with invariant mass closest to Z mass
    auto idx_cmb = Combinations(GoodJets_idx, 2);
    auto best_mass = -1;
    size_t best_i1 = 0; size_t best_i2 = 0;
    for (size_t i = 0; i < idx_cmb[0].size(); i++) {
        const auto i1 = idx_cmb[0][i];
        const auto i2 = idx_cmb[1][i];
        if (abs(eta[i1] - eta[i2]) >= DELTAETA_JJ_CUT) {
            ROOT::Math::PtEtaPhiMVector p1(pt[i1], eta[i1], phi[i1], mass[i1]);
            ROOT::Math::PtEtaPhiMVector p2(pt[i2], eta[i2], phi[i2], mass[i2]);
            const auto this_mass = (p1 + p2).M();
            if (this_mass > best_mass) {
                best_mass = this_mass;
                best_i1 = i1;
                best_i2 = i2;
            }
        }
    } 
    idx.emplace_back(best_i1);
    idx.emplace_back(best_i2);

    return idx;
}

RVec<int> SelectLepton(rvec_f lepton_pt, rvec_f lepton_eta, rvec_f lepton_phi, rvec_i lepton_pdgId, rvec_i lepton_tightId, rvec_i lepton_looseId,
                 rvec_f lepton_jetRelIso, rvec_f lepton_pfRelIso04_all, 
                 rvec_f lepton_mvaFall17V2Iso_WPL, rvec_f lepton_mvaFall17V2Iso_WP90, 
                 rvec_f jet_pt, rvec_f jet_eta, rvec_f jet_phi, rvec_i VBSJets_idx)
{

    //setting jet-related quantities if isolation from them is needed
    size_t jet1_idx = VBSJets_idx[0];
    size_t jet2_idx = VBSJets_idx[1];

    float jet1eta = jet_eta[jet1_idx];
    float jet2eta = jet_eta[jet2_idx];
    float jet1phi = jet_phi[jet1_idx];
    float jet2phi = jet_phi[jet2_idx];
    float isocone = DR_OVERLAP_CONE_OTHER;

    RVec<size_t> Tleptons_idx;
    RVec<size_t> LnTleptons_idx;

    for (size_t i = 0; i < idx_cmb[0].size(); i++) {
        //setting loose and tight, eta, and pt criteria for leptons depending on lepton flavour
        lepflav = abs(lepton_pdgId);

        if (lepflav == 11){
            auto IsLooseID = lepton_mvaFall17V2Iso_WPL[i];
            auto IsTightID = lepton_mvaFall17V2Iso_WP90[i];
            auto IsTightIso = bool(lepton_jetRelIso[i]<ISO_CUT_ELE and lepton_jetRelIso[i]>=0.);
            auto IsLooseIso = bool(lepton_jetRelIso[i]<1. and lepton_jetRelIso[i]>=0.);
            auto IsInEtaRegion = bool(abs(lepton_eta[i])<ETA_CUT_ELE && !(abs(lepton_eta)>1.4442 and abs(lepton_eta)<1.566));
            auto IsInPtRegion = lepton_pt > PT_CUT_ELE;
        }
        else if (lepflav == 13){
            auto IsTightID = lepton_tightId[i];
            auto IsLooseID = lepton_looseId[i];
            auto IsTightIso = bool(lepton_pfRelIso04_all[i]<ISO_CUT_MU and lepton_pfRelIso04_all[i]>=0.);
            auto IsLooseIso = bool(lepton_pfRelIso04_all[i]<1. and lepton_pfRelIso04_all[i]>=0.);
            auto IsInEtaRegion = abs(lepton_eta[i]) < ETA_CUT_MU;
            auto IsInPtRegion = lepton_pt[i] > PT_CUT_MU;
        }

        //find tight and loose-not-tight leptons filtering with jet-lep isolation criteria
        if (IsInEtaRegion && IsInPtRegion){
            if(IsLooseID && IsLooseIso){
                if(deltaR(lepton_eta, lepton_phi, jet1eta, jet1phi) > isocone && deltaR(lepton_eta, lepton_phi, jet2eta, jet2phi) > isocone){
                    if(IsTightID && IsTightIso) Tleptons_idx.push_back(i)
                    else LnTleptons_idx.push_back(i)
                }
            }
        }
    }
    //select leading tight/loose-not-tight lepton
    RVec<int> idx(2);
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


RVec<int> SelectAndVetoTaus(rvec_f Tau_eta, rvec_f Tau_phi,\ 
                            rvec_f Tau_idDeepTau2017v2p1VSjet, rvec_f Tau_idDeepTau2017v2p1VSe, rvec_f Tau_idDeepTau2017v2p1VSmu, rvec_f Tau_idDecayModeNewDMs,\ 
                            rvec_f Lepton_eta, rvec_f Lepton_phi, rvec_i Lepton_pdgId, rvec_i Lepton_idx, rvec_f Jet_eta, rvec_f Jet_phi, rvec_i Jet_idx)
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
    
    for (size_t i = 0; i < Tau_eta.size(); i++) {
        if (abs(Lepton_pdgId[Lepton_idx[0]])==11) float cutloose_vsjet = ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_ELE;
        else if (abs(Lepton_pdgId[Lepton_idx[0]])==13) float cutloose_vsjet = ID_TAU_RECO_DEEPTAU_VSJET_LOOSE_MU;

        if ((Tau_idDeepTau2017v2p1VSjet[i]>=cutloose_vsjet && Tau_idDeepTau2017v2p1VSe[i]>=ID_TAU_RECO_DEEPTAU_VSELE && Tau_idDeepTau2017v2p1VSmu[i]>=ID_TAU_RECO_DEEPTAU_VSMU && Tau_idDecayModeNewDMs[i])\ 
            && deltaR(Tau_eta[i], Tau_phi[i], Lepton_eta[Lepton_idx[0]], Lepton_phi[Lepton_idx[0]])>DR_OVERLAP_CONE_TAU\ 
            && deltaR(Tau_eta[i], Tau_phi[i], jet1eta, jet1phi)>isocone\ 
            && deltaR(Tau_eta[i], Tau_phi[i], jet2eta, jet2phi)>isocone && Tau_pt[i]>=PT_CUT_TAU\ 
            && abs(Tau_eta[i])<=ETA_CUT_TAU){
            
            nTau++;

            if(Tau_idDeepTau2017v2p1VSjet>=ID_TAU_RECO_DEEPTAU_VSJET){
                idx[0] = i;
                idx[1] = 0;
            } 
            else:{
                idx[0] = i;
                idx[1] = 1;
            }
        }
    
    if(nTau!=1) idx[1] = -1;                                                                                               
    
    return idx;
}

bool BVeto(rvec_f Jet_pt, rvec_f Jet_eta, rvec_i GoodJet_idx)
{
    bool veto = False;
    for (size_t i = 0; i < GoodJet_idx.size(); i++) {
        if (Jet_btagDeepFlavB[i]>=WP_btagger[BTAG_ALGO][BTAG_WP])*(Jet_pt[i]>BTAG_PT_CUT)*(abs(Jet_eta[i])<BTAG_ETA_CUT):
            return true;
    return false;
}