//#ifndef CMSJMECalculators_JMESystematicsCalculators_H
//#define CMSJMECalculators_JMESystematicsCalculators_H

#include <map>
#include <ROOT/RVec.hxx>
/*
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrectorCalculator.h"
#include "CommonTools/Utils/interface/FormulaEvaluator.h"
#include "JetMETCorrections/Modules/src/JetResolution.cc"
#include "CondFormats/JetMETObjects/src/JetCorrectorParameters.cc"
#include "CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc"
#include "CondFormats/JetMETObjects/src/FactorizedJetCorrectorCalculator.cc"
#include "CommonTools/Utils/src/FormulaEvaluator.cc"
*/
#include "CondFormats/JetMETObjects/interface/Utilities.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParametersHelper.h"
#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrectorCalculator.h"
#include "CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc"
#include "CondFormats/JetMETObjects/src/JetCorrectorParameters.cc"
#include "CondFormats/JetMETObjects/src/JetCorrectorParametersHelper.cc"
#include "CondFormats/JetMETObjects/src/JetResolutionObject.cc"
#include "CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc"
#include "CondFormats/JetMETObjects/src/SimpleJetCorrector.cc"
#include "CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc"
#include "CondFormats/JetMETObjects/src/FactorizedJetCorrectorCalculator.cc"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "JetMETCorrections/Modules/src/JetResolution.cc"
#include "CommonTools/Utils/interface/FormulaEvaluator.h"
#include "CommonTools/Utils/src/FormulaEvaluator.cc"
#include "CommonTools/Utils/src/formulaBinaryOperatorEvaluator.h"
#include "CommonTools/Utils/src/formulaFunctionOneArgEvaluator.h"
#include "CommonTools/Utils/src/formulaFunctionTwoArgsEvaluator.h"
#include "CommonTools/Utils/src/formulaUnaryMinusEvaluator.h"
#include "CommonTools/Utils/src/formulaConstantEvaluator.h"
#include "CommonTools/Utils/src/formulaConstantEvaluator.cc"
#include "CommonTools/Utils/src/formulaEvaluatorBase.h"
#include "CommonTools/Utils/src/formulaEvaluatorBase.cc"
#include "CommonTools/Utils/src/formulaParameterEvaluator.h"
#include "CommonTools/Utils/src/formulaParameterEvaluator.cc"
#include "CommonTools/Utils/src/formulaVariableEvaluator.h"
#include "CommonTools/Utils/src/formulaVariableEvaluator.cc"

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

//#endif // CMSJMECalculators_JMESystematicsCalculators_H
//#include "JMESystematicsCalculators.h"
//#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrectorCalculator.h"

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
    //cout<<"ciao"<<endl;
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
