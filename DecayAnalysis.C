#define DecayAnalysis_cxx
#include "DecayAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>


DecayAnalysis::DecayAnalysis(TTree *tree) : fChain(0)  {
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("sig_510664_mc20e_EL_LooseBLayerLH_Loose_VarRad_MU_Loose_Loose_VarRad.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("sig_510664_mc20e_EL_LooseBLayerLH_Loose_VarRad_MU_Loose_Loose_VarRad.root");
      }
      f->GetObject("reco",tree);

   }
   Init(tree);
}

DecayAnalysis::~DecayAnalysis() {
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t DecayAnalysis::GetEntry(Long64_t entry) {
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t DecayAnalysis::LoadTree(Long64_t entry) {
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void DecayAnalysis::Init(TTree *tree)
{  // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   el_IFFClass = 0;
   el_charge = 0;
   el_d0sig_loose = 0;
   el_eta = 0;
   el_phi = 0;
   el_z0sintheta_loose = 0;
   jet_eta = 0;
   jet_DL1dv01_Continuous_quantile = 0;
   jet_DL1dv01_FixedCutBEff_70_select = 0;
   jet_DL1dv01_FixedCutBEff_85_select = 0;
   jet_phi = 0;
   mu_IFFClass = 0;
   mu_charge = 0;
   mu_d0sig_loose = 0;
   mu_eta = 0;
   mu_phi = 0;
   mu_z0sintheta_loose = 0;
   tau_NNDecayMode = 0;
   tau_charge = 0;
   tau_eta = 0;
   tau_phi = 0;
   el_select_loose_NOSYS = 0;
   el_select_tight_NOSYS = 0;
   el_e_NOSYS = 0;
   el_select_passesOR_NOSYS = 0;
   el_pt_NOSYS = 0;
   jet_select_baselineJvt_NOSYS = 0;
   jet_e_NOSYS = 0;
   jet_select_passesOR_NOSYS = 0;
   jet_pt_NOSYS = 0;
   mu_select_loose_NOSYS = 0;
   mu_select_tight_NOSYS = 0;
   mu_e_NOSYS = 0;
   mu_TTVA_effSF_loose_NOSYS = 0;
   mu_TTVA_effSF_tight_NOSYS = 0;
   mu_isol_effSF_loose_NOSYS = 0;
   mu_isol_effSF_tight_NOSYS = 0;
   mu_reco_effSF_loose_NOSYS = 0;
   mu_reco_effSF_tight_NOSYS = 0;
   mu_select_passesOR_NOSYS = 0;
   mu_pt_NOSYS = 0;
   tau_select_tight_NOSYS = 0;
   tau_e_NOSYS = 0;
   tau_select_passesOR_NOSYS = 0;
   tau_pt_NOSYS = 0;
   tau_EvetoFakeTau_effSF_tight_NOSYS = 0;
   tau_EvetoTrueTau_effSF_tight_NOSYS = 0;
   tau_ID_effSF_tight_NOSYS = 0;
   tau_Reco_effSF_tight_NOSYS = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("weight_beamspot", &weight_beamspot, &b_weight_beamspot);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("mcChannelNumber", &mcChannelNumber, &b_mcChannelNumber);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("trigPassed_HLT_2e12_lhloose_L12EM10VH", &trigPassed_HLT_2e12_lhloose_L12EM10VH, &b_trigPassed_HLT_2e12_lhloose_L12EM10VH);
   fChain->SetBranchAddress("trigPassed_HLT_2e17_lhvloose_nod0", &trigPassed_HLT_2e17_lhvloose_nod0, &b_trigPassed_HLT_2e17_lhvloose_nod0);
   fChain->SetBranchAddress("trigPassed_HLT_2e17_lhvloose_nod0_L12EM15VHI", &trigPassed_HLT_2e17_lhvloose_nod0_L12EM15VHI, &b_trigPassed_HLT_2e17_lhvloose_nod0_L12EM15VHI);
   fChain->SetBranchAddress("trigPassed_HLT_2e24_lhvloose_nod0", &trigPassed_HLT_2e24_lhvloose_nod0, &b_trigPassed_HLT_2e24_lhvloose_nod0);
   fChain->SetBranchAddress("trigPassed_HLT_2mu10", &trigPassed_HLT_2mu10, &b_trigPassed_HLT_2mu10);
   fChain->SetBranchAddress("trigPassed_HLT_e17_lhloose_mu14", &trigPassed_HLT_e17_lhloose_mu14, &b_trigPassed_HLT_e17_lhloose_mu14);
   fChain->SetBranchAddress("trigPassed_HLT_e17_lhloose_nod0_mu14", &trigPassed_HLT_e17_lhloose_nod0_mu14, &b_trigPassed_HLT_e17_lhloose_nod0_mu14);
   fChain->SetBranchAddress("trigPassed_HLT_e7_lhmedium_mu24", &trigPassed_HLT_e7_lhmedium_mu24, &b_trigPassed_HLT_e7_lhmedium_mu24);
   fChain->SetBranchAddress("trigPassed_HLT_mu18_mu8noL1", &trigPassed_HLT_mu18_mu8noL1, &b_trigPassed_HLT_mu18_mu8noL1);
   fChain->SetBranchAddress("trigPassed_HLT_mu22_mu8noL1", &trigPassed_HLT_mu22_mu8noL1, &b_trigPassed_HLT_mu22_mu8noL1);
   fChain->SetBranchAddress("el_IFFClass", &el_IFFClass, &b_el_IFFClass);
   fChain->SetBranchAddress("el_charge", &el_charge, &b_el_charge);
   fChain->SetBranchAddress("el_d0sig_loose", &el_d0sig_loose, &b_el_d0sig_loose);
   fChain->SetBranchAddress("el_eta", &el_eta, &b_el_eta);
   fChain->SetBranchAddress("el_phi", &el_phi, &b_el_phi);
   fChain->SetBranchAddress("el_z0sintheta_loose", &el_z0sintheta_loose, &b_el_z0sintheta_loose);
   fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_DL1dv01_Continuous_quantile", &jet_DL1dv01_Continuous_quantile, &b_jet_DL1dv01_Continuous_quantile);
   fChain->SetBranchAddress("jet_DL1dv01_FixedCutBEff_70_select", &jet_DL1dv01_FixedCutBEff_70_select, &b_jet_DL1dv01_FixedCutBEff_70_select);
   fChain->SetBranchAddress("jet_DL1dv01_FixedCutBEff_85_select", &jet_DL1dv01_FixedCutBEff_85_select, &b_jet_DL1dv01_FixedCutBEff_85_select);
   fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("mu_IFFClass", &mu_IFFClass, &b_mu_IFFClass);
   fChain->SetBranchAddress("mu_charge", &mu_charge, &b_mu_charge);
   fChain->SetBranchAddress("mu_d0sig_loose", &mu_d0sig_loose, &b_mu_d0sig_loose);
   fChain->SetBranchAddress("mu_eta", &mu_eta, &b_mu_eta);
   fChain->SetBranchAddress("mu_phi", &mu_phi, &b_mu_phi);
   fChain->SetBranchAddress("mu_z0sintheta_loose", &mu_z0sintheta_loose, &b_mu_z0sintheta_loose);
   fChain->SetBranchAddress("tau_NNDecayMode", &tau_NNDecayMode, &b_tau_NNDecayMode);
   fChain->SetBranchAddress("tau_charge", &tau_charge, &b_tau_charge);
   fChain->SetBranchAddress("tau_eta", &tau_eta, &b_tau_eta);
   fChain->SetBranchAddress("tau_phi", &tau_phi, &b_tau_phi);
   fChain->SetBranchAddress("weight_pileup_NOSYS", &weight_pileup_NOSYS, &b_weight_pileup_NOSYS);
   fChain->SetBranchAddress("weight_ftag_effSF_DL1dv01_Continuous_NOSYS", &weight_ftag_effSF_DL1dv01_Continuous_NOSYS, &b_weight_ftag_effSF_DL1dv01_Continuous_NOSYS);
   fChain->SetBranchAddress("weight_mc_NOSYS", &weight_mc_NOSYS, &b_weight_mc_NOSYS);
   fChain->SetBranchAddress("weight_jvt_effSF_NOSYS", &weight_jvt_effSF_NOSYS, &b_weight_jvt_effSF_NOSYS);
   fChain->SetBranchAddress("weight_leptonSF_loose_NOSYS", &weight_leptonSF_loose_NOSYS, &b_weight_leptonSF_loose_NOSYS);
   fChain->SetBranchAddress("pass_SR2L2J_NOSYS", &pass_SR2L2J_NOSYS, &b_pass_SR2L2J_NOSYS);
   fChain->SetBranchAddress("pass_SR2L3J_NOSYS", &pass_SR2L3J_NOSYS, &b_pass_SR2L3J_NOSYS);
   fChain->SetBranchAddress("pass_SR2L4J_NOSYS", &pass_SR2L4J_NOSYS, &b_pass_SR2L4J_NOSYS);
   fChain->SetBranchAddress("pass_SR2LOS2J_NOSYS", &pass_SR2LOS2J_NOSYS, &b_pass_SR2LOS2J_NOSYS);
   fChain->SetBranchAddress("pass_SR2LOS3J_NOSYS", &pass_SR2LOS3J_NOSYS, &b_pass_SR2LOS3J_NOSYS);
   fChain->SetBranchAddress("pass_SR2LOS4J_NOSYS", &pass_SR2LOS4J_NOSYS, &b_pass_SR2LOS4J_NOSYS);
   fChain->SetBranchAddress("pass_SR2LSS2J_NOSYS", &pass_SR2LSS2J_NOSYS, &b_pass_SR2LSS2J_NOSYS);
   fChain->SetBranchAddress("pass_SR2LSS3J_NOSYS", &pass_SR2LSS3J_NOSYS, &b_pass_SR2LSS3J_NOSYS);
   fChain->SetBranchAddress("pass_SR2LSS4J_NOSYS", &pass_SR2LSS4J_NOSYS, &b_pass_SR2LSS4J_NOSYS);
   fChain->SetBranchAddress("pass_SUBcommon_NOSYS", &pass_SUBcommon_NOSYS, &b_pass_SUBcommon_NOSYS);
   fChain->SetBranchAddress("el_select_loose_NOSYS", &el_select_loose_NOSYS, &b_el_select_loose_NOSYS);
   fChain->SetBranchAddress("el_select_tight_NOSYS", &el_select_tight_NOSYS, &b_el_select_tight_NOSYS);
   fChain->SetBranchAddress("el_e_NOSYS", &el_e_NOSYS, &b_el_e_NOSYS);
   fChain->SetBranchAddress("el_select_passesOR_NOSYS", &el_select_passesOR_NOSYS, &b_el_select_passesOR_NOSYS);
   fChain->SetBranchAddress("el_pt_NOSYS", &el_pt_NOSYS, &b_el_pt_NOSYS);
   fChain->SetBranchAddress("jet_select_baselineJvt_NOSYS", &jet_select_baselineJvt_NOSYS, &b_jet_select_baselineJvt_NOSYS);
   fChain->SetBranchAddress("jet_e_NOSYS", &jet_e_NOSYS, &b_jet_e_NOSYS);
   fChain->SetBranchAddress("jet_select_passesOR_NOSYS", &jet_select_passesOR_NOSYS, &b_jet_select_passesOR_NOSYS);
   fChain->SetBranchAddress("jet_pt_NOSYS", &jet_pt_NOSYS, &b_jet_pt_NOSYS);
   fChain->SetBranchAddress("mu_select_loose_NOSYS", &mu_select_loose_NOSYS, &b_mu_select_loose_NOSYS);
   fChain->SetBranchAddress("mu_select_tight_NOSYS", &mu_select_tight_NOSYS, &b_mu_select_tight_NOSYS);
   fChain->SetBranchAddress("mu_e_NOSYS", &mu_e_NOSYS, &b_mu_e_NOSYS);
   fChain->SetBranchAddress("mu_TTVA_effSF_loose_NOSYS", &mu_TTVA_effSF_loose_NOSYS, &b_mu_TTVA_effSF_loose_NOSYS);
   fChain->SetBranchAddress("mu_TTVA_effSF_tight_NOSYS", &mu_TTVA_effSF_tight_NOSYS, &b_mu_TTVA_effSF_tight_NOSYS);
   fChain->SetBranchAddress("mu_isol_effSF_loose_NOSYS", &mu_isol_effSF_loose_NOSYS, &b_mu_isol_effSF_loose_NOSYS);
   fChain->SetBranchAddress("mu_isol_effSF_tight_NOSYS", &mu_isol_effSF_tight_NOSYS, &b_mu_isol_effSF_tight_NOSYS);
   fChain->SetBranchAddress("mu_reco_effSF_loose_NOSYS", &mu_reco_effSF_loose_NOSYS, &b_mu_reco_effSF_loose_NOSYS);
   fChain->SetBranchAddress("mu_reco_effSF_tight_NOSYS", &mu_reco_effSF_tight_NOSYS, &b_mu_reco_effSF_tight_NOSYS);
   fChain->SetBranchAddress("mu_select_passesOR_NOSYS", &mu_select_passesOR_NOSYS, &b_mu_select_passesOR_NOSYS);
   fChain->SetBranchAddress("mu_pt_NOSYS", &mu_pt_NOSYS, &b_mu_pt_NOSYS);
   fChain->SetBranchAddress("tau_select_tight_NOSYS", &tau_select_tight_NOSYS, &b_tau_select_tight_NOSYS);
   fChain->SetBranchAddress("tau_e_NOSYS", &tau_e_NOSYS, &b_tau_e_NOSYS);
   fChain->SetBranchAddress("tau_select_passesOR_NOSYS", &tau_select_passesOR_NOSYS, &b_tau_select_passesOR_NOSYS);
   fChain->SetBranchAddress("tau_pt_NOSYS", &tau_pt_NOSYS, &b_tau_pt_NOSYS);
   fChain->SetBranchAddress("tau_EvetoFakeTau_effSF_tight_NOSYS", &tau_EvetoFakeTau_effSF_tight_NOSYS, &b_tau_EvetoFakeTau_effSF_tight_NOSYS);
   fChain->SetBranchAddress("tau_EvetoTrueTau_effSF_tight_NOSYS", &tau_EvetoTrueTau_effSF_tight_NOSYS, &b_tau_EvetoTrueTau_effSF_tight_NOSYS);
   fChain->SetBranchAddress("tau_ID_effSF_tight_NOSYS", &tau_ID_effSF_tight_NOSYS, &b_tau_ID_effSF_tight_NOSYS);
   fChain->SetBranchAddress("tau_Reco_effSF_tight_NOSYS", &tau_Reco_effSF_tight_NOSYS, &b_tau_Reco_effSF_tight_NOSYS);
   fChain->SetBranchAddress("met_met_NOSYS", &met_met_NOSYS, &b_met_met_NOSYS);
   fChain->SetBranchAddress("met_phi_NOSYS", &met_phi_NOSYS, &b_met_phi_NOSYS);
   fChain->SetBranchAddress("met_significance_NOSYS", &met_significance_NOSYS, &b_met_significance_NOSYS);
   fChain->SetBranchAddress("met_sumet_NOSYS", &met_sumet_NOSYS, &b_met_sumet_NOSYS);
   Notify();
}

bool DecayAnalysis::Notify() {
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void DecayAnalysis::Show(Long64_t entry) {
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t DecayAnalysis::Cut(Long64_t entry) {
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

//########################################

double DecayAnalysis::CalcInvMass(vector<ROOT::Math::PtEtaPhiEVector>& particles){
   
   double tot_px{};
   double tot_py{};
   double tot_pz{};
   double tot_e{};

   for (int i = 0; i < particles.size(); i++){
      tot_px += particles[i].Px(); 
      tot_py += particles[i].Py(); 
      tot_pz += particles[i].Pz();
      tot_e  += particles[i].E();
   } 

   double tot_p2 = (tot_px * tot_px + tot_py * tot_py + tot_pz * tot_pz);
   double tot_e2 = tot_e * tot_e;
   //std::cout << "p2 " << tot_p2 << '\n';
   //std::cout << "E2 " << tot_e2 << '\n';

   double inv_mass = TMath::Sqrt(tot_e2 - tot_p2);
   std::cout << "inv_mass " << inv_mass/1e3 << '\n';

   return inv_mass/1e3;    
}

void DecayAnalysis::DrawHistos(){
   TCanvas *c1 = new TCanvas("c1", "Momentum of the particles", 800, 600);
   TCanvas *c2 = new TCanvas("c2", "Invariant masses", 800, 600);

   c1->Divide(2, 2);
   c2->Divide(1, 2);

   TH1F* h_e_pT = new TH1F("h_ele_pT", "Electron pT; p_{T} [GeV]; Entries", 100, 0, 500); // 100 bins from 0 to 500 GeV
   for (const auto& vec : electrons) {
      //std::cout << vec.Pt() << std::endl;
      h_e_pT->Fill(vec.Pt()/1e3); //in GeV
   }
   c1->cd(1);
   h_e_pT->Draw();

   TH1F* h_mu_pT = new TH1F("h_muon_pT", "Muon pT; p_{T} [GeV]; Entries", 100, 0, 500); // 100 bins from 0 to 500 GeV
   for (const auto& vec : muons) {
      //std::cout << vec.Pt() << std::endl;
      h_mu_pT->Fill(vec.Pt()/1e3); //in GeV
   }
   c1->cd(2);
   h_mu_pT->Draw();

   TH1F* h_jet_pT = new TH1F("h_jet_pT", "Jet pT; p_{T} [GeV]; Entries", 100, 0, 500); // 100 bins from 0 to 500 GeV
   for (const auto& vec : jets) {
      //std::cout << vec.Pt() << std::endl;
      h_jet_pT->Fill(vec.Pt()/1e3); //in GeV
   }

   c1->cd(3);
   h_jet_pT->Draw();

   TH1F* leptons = new TH1F("particles", "Leptons; Leptons combinations; Entries", 7, 0.5, 3.5);

   leptons->GetXaxis()->SetBinLabel(2, "2e");
   leptons->GetXaxis()->SetBinLabel(4, "1e1mu");
   leptons->GetXaxis()->SetBinLabel(6, "2mu");
   leptons->SetMinimum(0);

   for (int i = 0 ; i < e2_entries; i++) {leptons->Fill(1);}
   for (int i = 0 ; i < e_entries; i++) {leptons->Fill(2);}
   for (int i = 0 ; i < mu2_entries; i++) {leptons->Fill(3);}

   c1->cd(4);
   leptons->Draw();

   TH1F* all_inv_masses_hist = new TH1F("all_inv_masses", "All invariant masses; m [GeV]; Entries", 100, 0, 500); // 100 bins from 0 to 500 GeV 
   for (const auto& vec : all_inv_masses) {
      all_inv_masses_hist->Fill(vec); 
      //std::cout << vec << '\n';
   }
   TH1F* inv_masses_hist = new TH1F("inv_masses", "Selected invariant masses; m [GeV]; Entries", 100, 0, 500); // 100 bins from 0 to 500 GeV 
   for (const auto& vec : inv_masses) {
      inv_masses_hist->Fill(vec); 
      //std::cout << vec << '\n';
   }

   c2->cd(1);
   all_inv_masses_hist->Draw();

   c2->cd(2);
   inv_masses_hist->Draw();

   c1->Update();
   c2->Update();

   TFile *f = new TFile("graphs.root", "RECREATE");

   h_e_pT->Write();
   h_mu_pT->Write();
   h_jet_pT->Write();
   leptons->Write();
   all_inv_masses_hist->Write();
   inv_masses_hist->Write();

   f->Close();
}

void DecayAnalysis::InitFilterMap() {
   filter_map["2L4J"] = &pass_SR2L4J_NOSYS;
   filter_map["2L3J"] = &pass_SR2L3J_NOSYS;
   filter_map["2L2J"] = &pass_SR2L2J_NOSYS;
}

bool DecayAnalysis::CheckFilter(const std::string& filter) {
   auto it = filter_map.find(filter);
   if (it != filter_map.end()) {
      return static_cast<bool>(*(it->second));  //prendere il valore come bool
   }
   std::cout << "you set a non valid filter!" << std::endl;
   return false;  
}

void DecayAnalysis::Loop()
{
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   std::string filter = "2L4J"; //scegli tra 2L4J, 2L3J, 2L2J
   InitFilterMap(); 

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      std::vector<ROOT::Math::PtEtaPhiEVector> particles;
      std::vector<std::string> particle_types;

      bool pass = CheckFilter(filter);
      if (pass) {
         for (size_t i=0; i<el_pt_NOSYS->size(); i++){
            ROOT::Math::PtEtaPhiEVector vecTmp1(el_pt_NOSYS->at(i), el_eta->at(i), el_phi->at(i), el_e_NOSYS->at(i));
            electrons.push_back(vecTmp1);
            particles.push_back(vecTmp1);
            particle_types.push_back("electron");
         }
         if (el_pt_NOSYS->size() == 2) {
            ++e2_entries;
         } else if (el_pt_NOSYS->size() == 1) {
            ++e_entries;
         }

         for (size_t i=0; i<mu_pt_NOSYS->size(); i++){
            ROOT::Math::PtEtaPhiEVector vecTmp2(mu_pt_NOSYS->at(i), mu_eta->at(i) , mu_phi->at(i), mu_e_NOSYS->at(i));
            muons.push_back(vecTmp2);
            particles.push_back(vecTmp2);
            particle_types.push_back("muon");
         }
         if (mu_pt_NOSYS->size() == 2) {++mu2_entries;}

         for (size_t i=0; i<jet_pt_NOSYS->size(); i++){
            ROOT::Math::PtEtaPhiEVector vecTmp3(jet_pt_NOSYS->at(i), jet_eta->at(i) , jet_phi->at(i), jet_e_NOSYS->at(i));
            jets.push_back(vecTmp3);
            particles.push_back(vecTmp3);

         }
      }

      if (particles.size() < 3) {
         std::cout << "There aren't enough particles" << '\n' 
                   << "---------" << '\n';
         particles.clear();
         particle_types.clear();
         continue;
      }

      std::cout << "Event " << jentry << ": "
                  << particle_types[0] << " and " << particle_types[1] << '\n';
      //for(int i = 0; i < particle_types.size(); i++) {
      //   std::cout << particle_types[i] << '\n';
      //}
      double inv_mass1{1e6}; //perchè con 0 non si aggiornerebbero mai
      double inv_mass2{0};

      if (particles.size() > 3) {
         std::vector<double> masses_i0;
         std::vector<double> masses_i1;

         for (int i = 0; i < 2; i++) {
            ROOT::Math::PtEtaPhiEVector& p1 = particles[i];

            for (int j = 2; j < particles.size(); j++) {
               ROOT::Math::PtEtaPhiEVector& p2 = particles[j];

               for (int k = j + 1; k < particles.size(); k++){
                  ROOT::Math::PtEtaPhiEVector& p3 = particles[k];
            
                  std::vector<ROOT::Math::PtEtaPhiEVector> combination = {p1, p2, p3};
                  double tempInvMass = CalcInvMass(combination); 

                  all_inv_masses.push_back(tempInvMass);

                  //std::cout << tempInvMass << '\n';

                  if (i == 0) {masses_i0.push_back(tempInvMass);}
                  if (i == 1) {masses_i1.push_back(tempInvMass);}

                  // if (i == 0 && TMath::Abs(tempInvMass - 65 * 1e3) < TMath::Abs(inv_mass1 - 65)){
                  //    inv_mass1 = tempInvMass;
                  //    comb1 = {p1, p2, p3};
                  // }
               
                  // if (i == 1 && TMath::Abs(tempInvMass - 65 * 1e3) < TMath::Abs(inv_mass2 - 65)){
                  //    inv_mass2 = tempInvMass;
                  //    comb2  = {p1, p2, p3};
                  // }            
               }
            }   
         }
            if(filter == "2L4J"){
               for (int i = 0; i < masses_i0.size(); i++) { 
                  for (int j = 0; j < masses_i1.size(); j++) {
                     if(TMath::Abs(masses_i0[i] - masses_i1[j]) < TMath::Abs(inv_mass1 - inv_mass2)) {
                        inv_mass1 = masses_i0[i];
                        inv_mass2 = masses_i1[j];
                     }
                  }  
               }
               inv_masses.push_back(inv_mass1);
               inv_masses.push_back(inv_mass2);
         }  
      }

      particles.clear();
      particle_types.clear();
   } //END LOOP IN ENTRIES

   DrawHistos();
}