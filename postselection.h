/// Header file with functions needed to execute the Python version
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

using namespace ROOT::VecOps;
using RNode = ROOT::RDF::RNode;
using rvec_f = const RVec<float> &;
using rvec_i = const RVec<int> &;

const float thr = 50; 
const float DELTAETA_JJ_CUT = 0.4;
const float DR_OVERLAP_CONE_OTHER = 1;
const float ETA_CUT_MU = 1;
const float ISO_CUT_MU = 1;
const float ETA_CUT_ELE = 1;
const float ISO_CUT_ELE = 1;

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
}

def SelectLepton(rvec_f lepton_pt, rvec_f lepton_eta, rvec_f lepton_phi, rvec_i lepton_pdgId, rvec_i lepton_tightId, rvec_i lepton_looseId,
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
    if (Tleptons_idx.size() == 0){
        idx[0] = Tleptons_idx[0];
        idx[1] = -1;
    }
    else if (LnTleptons_idx.size() == 0){
        idx[0] = -1;
        idx[1] = Tleptons_idx[0];
    }
    else{
        idx[0] = -1;
        idx[1] = -1;  
    }
    
    return idx;
}
