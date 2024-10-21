//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Sep 26 11:33:50 2024 by ROOT version 6.32.02
// from TTree reco/xAOD->NTuple tree
// found on file: testXX.root
//////////////////////////////////////////////////////////

#ifndef MyClass_h
#define MyClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include "Math/Vector4D.h"


// Header file for the classes stored in the TTree if any.
#include "vector"


class MyClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         weight_beamspot;
   ULong64_t       eventNumber;
   UInt_t          mcChannelNumber;
   UInt_t          runNumber;
   Bool_t          trigPassed_HLT_2e12_lhloose_L12EM10VH;
   Bool_t          trigPassed_HLT_2e17_lhvloose_nod0;
   Bool_t          trigPassed_HLT_2e17_lhvloose_nod0_L12EM15VHI;
   Bool_t          trigPassed_HLT_2e24_lhvloose_nod0;
   Bool_t          trigPassed_HLT_2mu10;
   Bool_t          trigPassed_HLT_e17_lhloose_mu14;
   Bool_t          trigPassed_HLT_e17_lhloose_nod0_mu14;
   Bool_t          trigPassed_HLT_e7_lhmedium_mu24;
   Bool_t          trigPassed_HLT_mu18_mu8noL1;
   Bool_t          trigPassed_HLT_mu22_mu8noL1;
   vector<int>     *el_IFFClass;
   vector<float>   *el_charge;
   vector<float>   *el_d0sig_loose;
   vector<float>   *el_eta;
   vector<float>   *el_phi;
   vector<float>   *el_z0sintheta_loose;
   vector<float>   *jet_eta;
   vector<int>     *jet_DL1dv01_Continuous_quantile;
   vector<char>    *jet_DL1dv01_FixedCutBEff_70_select;
   vector<char>    *jet_DL1dv01_FixedCutBEff_85_select;
   vector<float>   *jet_phi;
   vector<int>     *mu_IFFClass;
   vector<float>   *mu_charge;
   vector<float>   *mu_d0sig_loose;
   vector<float>   *mu_eta;
   vector<float>   *mu_phi;
   vector<float>   *mu_z0sintheta_loose;
   vector<int>     *tau_NNDecayMode;
   vector<float>   *tau_charge;
   vector<float>   *tau_eta;
   vector<float>   *tau_phi;
   Float_t         weight_pileup_NOSYS;
   Float_t         weight_ftag_effSF_DL1dv01_Continuous_NOSYS;
   Float_t         weight_mc_NOSYS;
   Float_t         weight_jvt_effSF_NOSYS;
   Float_t         weight_leptonSF_loose_NOSYS;
   Char_t          pass_SR2L2J_NOSYS;
   Char_t          pass_SR2L3J_NOSYS;
   Char_t          pass_SR2L4J_NOSYS;
   Char_t          pass_SR2LOS2J_NOSYS;
   Char_t          pass_SR2LOS3J_NOSYS;
   Char_t          pass_SR2LOS4J_NOSYS;
   Char_t          pass_SR2LSS2J_NOSYS;
   Char_t          pass_SR2LSS3J_NOSYS;
   Char_t          pass_SR2LSS4J_NOSYS;
   Char_t          pass_SUBcommon_NOSYS;
   vector<char>    *el_select_loose_NOSYS;
   vector<char>    *el_select_tight_NOSYS;
   vector<float>   *el_e_NOSYS;
   vector<char>    *el_select_passesOR_NOSYS;
   vector<float>   *el_pt_NOSYS;
   vector<char>    *jet_select_baselineJvt_NOSYS;
   vector<float>   *jet_e_NOSYS;
   vector<char>    *jet_select_passesOR_NOSYS;
   vector<float>   *jet_pt_NOSYS;
   vector<char>    *mu_select_loose_NOSYS;
   vector<char>    *mu_select_tight_NOSYS;
   vector<float>   *mu_e_NOSYS;
   vector<float>   *mu_TTVA_effSF_loose_NOSYS;
   vector<float>   *mu_TTVA_effSF_tight_NOSYS;
   vector<float>   *mu_isol_effSF_loose_NOSYS;
   vector<float>   *mu_isol_effSF_tight_NOSYS;
   vector<float>   *mu_reco_effSF_loose_NOSYS;
   vector<float>   *mu_reco_effSF_tight_NOSYS;
   vector<char>    *mu_select_passesOR_NOSYS;
   vector<float>   *mu_pt_NOSYS;
   vector<char>    *tau_select_tight_NOSYS;
   vector<float>   *tau_e_NOSYS;
   vector<char>    *tau_select_passesOR_NOSYS;
   vector<float>   *tau_pt_NOSYS;
   vector<float>   *tau_EvetoFakeTau_effSF_tight_NOSYS;
   vector<float>   *tau_EvetoTrueTau_effSF_tight_NOSYS;
   vector<float>   *tau_ID_effSF_tight_NOSYS;
   vector<float>   *tau_Reco_effSF_tight_NOSYS;
   Float_t         met_met_NOSYS;
   Float_t         met_phi_NOSYS;
   Float_t         met_significance_NOSYS;
   Float_t         met_sumet_NOSYS;

   // List of branches
   TBranch        *b_weight_beamspot;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_mcChannelNumber;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_trigPassed_HLT_2e12_lhloose_L12EM10VH;   //!
   TBranch        *b_trigPassed_HLT_2e17_lhvloose_nod0;   //!
   TBranch        *b_trigPassed_HLT_2e17_lhvloose_nod0_L12EM15VHI;   //!
   TBranch        *b_trigPassed_HLT_2e24_lhvloose_nod0;   //!
   TBranch        *b_trigPassed_HLT_2mu10;   //!
   TBranch        *b_trigPassed_HLT_e17_lhloose_mu14;   //!
   TBranch        *b_trigPassed_HLT_e17_lhloose_nod0_mu14;   //!
   TBranch        *b_trigPassed_HLT_e7_lhmedium_mu24;   //!
   TBranch        *b_trigPassed_HLT_mu18_mu8noL1;   //!
   TBranch        *b_trigPassed_HLT_mu22_mu8noL1;   //!
   TBranch        *b_el_IFFClass;   //!
   TBranch        *b_el_charge;   //!
   TBranch        *b_el_d0sig_loose;   //!
   TBranch        *b_el_eta;   //!
   TBranch        *b_el_phi;   //!
   TBranch        *b_el_z0sintheta_loose;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_DL1dv01_Continuous_quantile;   //!
   TBranch        *b_jet_DL1dv01_FixedCutBEff_70_select;   //!
   TBranch        *b_jet_DL1dv01_FixedCutBEff_85_select;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_mu_IFFClass;   //!
   TBranch        *b_mu_charge;   //!
   TBranch        *b_mu_d0sig_loose;   //!
   TBranch        *b_mu_eta;   //!
   TBranch        *b_mu_phi;   //!
   TBranch        *b_mu_z0sintheta_loose;   //!
   TBranch        *b_tau_NNDecayMode;   //!
   TBranch        *b_tau_charge;   //!
   TBranch        *b_tau_eta;   //!
   TBranch        *b_tau_phi;   //!
   TBranch        *b_weight_pileup_NOSYS;   //!
   TBranch        *b_weight_ftag_effSF_DL1dv01_Continuous_NOSYS;   //!
   TBranch        *b_weight_mc_NOSYS;   //!
   TBranch        *b_weight_jvt_effSF_NOSYS;   //!
   TBranch        *b_weight_leptonSF_loose_NOSYS;   //!
   TBranch        *b_pass_SR2L2J_NOSYS;   //!
   TBranch        *b_pass_SR2L3J_NOSYS;   //!
   TBranch        *b_pass_SR2L4J_NOSYS;   //!
   TBranch        *b_pass_SR2LOS2J_NOSYS;   //!
   TBranch        *b_pass_SR2LOS3J_NOSYS;   //!
   TBranch        *b_pass_SR2LOS4J_NOSYS;   //!
   TBranch        *b_pass_SR2LSS2J_NOSYS;   //!
   TBranch        *b_pass_SR2LSS3J_NOSYS;   //!
   TBranch        *b_pass_SR2LSS4J_NOSYS;   //!
   TBranch        *b_pass_SUBcommon_NOSYS;   //!
   TBranch        *b_el_select_loose_NOSYS;   //!
   TBranch        *b_el_select_tight_NOSYS;   //!
   TBranch        *b_el_e_NOSYS;   //!
   TBranch        *b_el_select_passesOR_NOSYS;   //!
   TBranch        *b_el_pt_NOSYS;   //!
   TBranch        *b_jet_select_baselineJvt_NOSYS;   //!
   TBranch        *b_jet_e_NOSYS;   //!
   TBranch        *b_jet_select_passesOR_NOSYS;   //!
   TBranch        *b_jet_pt_NOSYS;   //!
   TBranch        *b_mu_select_loose_NOSYS;   //!
   TBranch        *b_mu_select_tight_NOSYS;   //!
   TBranch        *b_mu_e_NOSYS;   //!
   TBranch        *b_mu_TTVA_effSF_loose_NOSYS;   //!
   TBranch        *b_mu_TTVA_effSF_tight_NOSYS;   //!
   TBranch        *b_mu_isol_effSF_loose_NOSYS;   //!
   TBranch        *b_mu_isol_effSF_tight_NOSYS;   //!
   TBranch        *b_mu_reco_effSF_loose_NOSYS;   //!
   TBranch        *b_mu_reco_effSF_tight_NOSYS;   //!
   TBranch        *b_mu_select_passesOR_NOSYS;   //!
   TBranch        *b_mu_pt_NOSYS;   //!
   TBranch        *b_tau_select_tight_NOSYS;   //!
   TBranch        *b_tau_e_NOSYS;   //!
   TBranch        *b_tau_select_passesOR_NOSYS;   //!
   TBranch        *b_tau_pt_NOSYS;   //!
   TBranch        *b_tau_EvetoFakeTau_effSF_tight_NOSYS;   //!
   TBranch        *b_tau_EvetoTrueTau_effSF_tight_NOSYS;   //!
   TBranch        *b_tau_ID_effSF_tight_NOSYS;   //!
   TBranch        *b_tau_Reco_effSF_tight_NOSYS;   //!
   TBranch        *b_met_met_NOSYS;   //!
   TBranch        *b_met_phi_NOSYS;   //!
   TBranch        *b_met_significance_NOSYS;   //!
   TBranch        *b_met_sumet_NOSYS;   //!

   MyClass(TTree *tree=0);
   virtual ~MyClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);

   double CalcInvMass(std::vector<ROOT::Math::PtEtaPhiEVector>& particles);
   // double CalcInvMass(const ROOT::Math::PtEtaPhiEVector& particle);

   std::vector<ROOT::Math::PtEtaPhiEVector> electrons;
   std::vector<ROOT::Math::PtEtaPhiEVector> muons;
   std::vector<ROOT::Math::PtEtaPhiEVector> jets;
   std::vector<ROOT::Math::PtEtaPhiEVector> comb1;
   std::vector<ROOT::Math::PtEtaPhiEVector> comb2;
   int e2_entries{};
   int e_entries{};
   int mu2_entries{};
   double inv_mass1{};
   double inv_mass2{};
   void DrawHistos();
   void InitFilterMap();
   bool CheckFilter(const std::string& filter = "2L4J");
   std::unordered_map<std::string, Char_t*> filter_map;

};

#endif

#ifdef MyClass_cxx

#endif // #ifdef MyClass_cxx
