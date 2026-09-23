#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TObjArray.h"
#include "TString.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>

namespace
{
bool HasBranch(TTree& tree, const char* name)
{
  return tree.GetListOfBranches()->FindObject(name) != nullptr;
}
}

void check_h11b_reaction(const char* input = "../data/reaction_all.root",
                         const char* output_prefix = "h11b_reaction_check")
{
  TChain tr("tr");
  tr.Add(input);

  const char* required_branches[] = {
    "event",
    "resonance_id",
    "branch_id",
    "ex_8Be",
    "ex_max_8Be",
    "eaa_8Be",
    "e_alpha8Be",
    "cos_theta_primary_cm",
    "cos_chi_exit",
    "event_weight",
    "phi_lab_alpha1",
    "phi_lab_alpha2",
    "phi_lab_alpha3",
    "reaction",
  };

  bool missing_branch = false;
  for(const char* branch : required_branches){
    if(!HasBranch(tr, branch)){
      std::cerr << "Missing branch: " << branch << std::endl;
      missing_branch = true;
    }
  }
  if(missing_branch){
    std::cerr << "Stop: input file does not contain all diagnostic branches." << std::endl;
    return;
  }

  Long64_t event = 0;
  int resonance_id = 0;
  int branch_id = 0;
  double ex_8Be = 0.0;
  double ex_max_8Be = 0.0;
  double eaa_8Be = 0.0;
  double e_alpha8Be = 0.0;
  double cos_theta_primary_cm = 0.0;
  double cos_chi_exit = 0.0;
  double event_weight = 0.0;
  double phi_lab_alpha1 = 0.0;
  double phi_lab_alpha2 = 0.0;
  double phi_lab_alpha3 = 0.0;
  char reaction[128] = {};

  tr.SetBranchAddress("event", &event);
  tr.SetBranchAddress("resonance_id", &resonance_id);
  tr.SetBranchAddress("branch_id", &branch_id);
  tr.SetBranchAddress("ex_8Be", &ex_8Be);
  tr.SetBranchAddress("ex_max_8Be", &ex_max_8Be);
  tr.SetBranchAddress("eaa_8Be", &eaa_8Be);
  tr.SetBranchAddress("e_alpha8Be", &e_alpha8Be);
  tr.SetBranchAddress("cos_theta_primary_cm", &cos_theta_primary_cm);
  tr.SetBranchAddress("cos_chi_exit", &cos_chi_exit);
  tr.SetBranchAddress("event_weight", &event_weight);
  tr.SetBranchAddress("phi_lab_alpha1", &phi_lab_alpha1);
  tr.SetBranchAddress("phi_lab_alpha2", &phi_lab_alpha2);
  tr.SetBranchAddress("phi_lab_alpha3", &phi_lab_alpha3);
  tr.SetBranchAddress("reaction", reaction);

  TH1D* h_resonance = new TH1D("h_resonance", "Resonance ID;resonance_id;counts", 800, 0.0, 800.0);
  TH1D* h_branch = new TH1D("h_branch", "Branch ID;branch_id;counts", 2, -0.5, 1.5);
  TH1D* h_ex = new TH1D("h_ex", "^{8}Be excitation energy;E_{x}(^{8}Be) [MeV];counts", 500, 0.0, 10.0);
  TH1D* h_eaa = new TH1D("h_eaa", "#alpha#alpha relative energy;E_{#alpha#alpha} [MeV];counts", 500, 0.0, 10.0);
  TH1D* h_e_alpha8Be =
    new TH1D("h_e_alpha8Be", "#alpha+^{8}Be relative energy;E_{#alpha^{8}Be} [MeV];counts", 500, 0.0, 10.0);
  TH1D* h_cos_theta =
    new TH1D("h_cos_theta", "Primary alpha angle;cos#theta_{primary}^{CM};counts", 120, -1.0, 1.0);
  TH1D* h_cos_chi = new TH1D("h_cos_chi", "Exit angular correlation;cos#chi;counts", 120, -1.0, 1.0);
  TH1D* h_weight = new TH1D("h_weight", "Event weight;event_weight;counts", 120, 0.0, 1.2);
  TH1D* h_phi = new TH1D("h_phi", "Lab alpha azimuths;#phi_{lab} [rad];counts", 120, -3.2, 3.2);
  TH2D* h_branch_res =
    new TH2D("h_branch_res", "Branch vs resonance;resonance_id;branch_id", 800, 0.0, 800.0, 2, -0.5, 1.5);
  TH2D* h_ex_limit =
    new TH2D("h_ex_limit", "8Be excitation limit;E_{x,max} [MeV];E_{x}(^{8}Be) [MeV]", 250, 0.0, 10.0, 250, 0.0, 10.0);

  Long64_t n_bad_ex_limit = 0;
  Long64_t n_bad_cos = 0;
  Long64_t n_bad_branch = 0;
  Long64_t n_bad_resonance = 0;
  Long64_t n_nonfinite = 0;

  Long64_t n_165 = 0;
  Long64_t n_675 = 0;
  Long64_t n_alpha0 = 0;
  Long64_t n_alpha1 = 0;
  double min_ex = 1.0e99;
  double max_ex = -1.0e99;
  double min_ex_max = 1.0e99;
  double max_ex_max = -1.0e99;
  double min_e_alpha8Be = 1.0e99;
  double max_e_alpha8Be = -1.0e99;

  const Long64_t nentries = tr.GetEntries();
  for(Long64_t i = 0; i < nentries; ++i){
    tr.GetEntry(i);

    const bool finite = std::isfinite(ex_8Be) && std::isfinite(ex_max_8Be) && std::isfinite(eaa_8Be) &&
                        std::isfinite(e_alpha8Be) && std::isfinite(cos_theta_primary_cm) &&
                        std::isfinite(cos_chi_exit) && std::isfinite(event_weight) &&
                        std::isfinite(phi_lab_alpha1) && std::isfinite(phi_lab_alpha2) &&
                        std::isfinite(phi_lab_alpha3);
    if(!finite){
      ++n_nonfinite;
      continue;
    }

    if(ex_8Be > ex_max_8Be + 1.0e-9) ++n_bad_ex_limit;
    if(cos_theta_primary_cm < -1.000001 || cos_theta_primary_cm > 1.000001 || cos_chi_exit < -1.000001 ||
       cos_chi_exit > 1.000001)
      ++n_bad_cos;
    if(branch_id != 0 && branch_id != 1) ++n_bad_branch;
    if(resonance_id != 165 && resonance_id != 675) ++n_bad_resonance;

    if(resonance_id == 165) ++n_165;
    if(resonance_id == 675) ++n_675;
    if(branch_id == 0) ++n_alpha0;
    if(branch_id == 1) ++n_alpha1;

    min_ex = std::min(min_ex, ex_8Be);
    max_ex = std::max(max_ex, ex_8Be);
    min_ex_max = std::min(min_ex_max, ex_max_8Be);
    max_ex_max = std::max(max_ex_max, ex_max_8Be);
    min_e_alpha8Be = std::min(min_e_alpha8Be, e_alpha8Be);
    max_e_alpha8Be = std::max(max_e_alpha8Be, e_alpha8Be);

    h_resonance->Fill(resonance_id);
    h_branch->Fill(branch_id);
    h_ex->Fill(ex_8Be);
    h_eaa->Fill(eaa_8Be);
    h_e_alpha8Be->Fill(e_alpha8Be);
    h_cos_theta->Fill(cos_theta_primary_cm);
    h_cos_chi->Fill(cos_chi_exit);
    h_weight->Fill(event_weight);
    h_phi->Fill(phi_lab_alpha1);
    h_phi->Fill(phi_lab_alpha2);
    h_phi->Fill(phi_lab_alpha3);
    h_branch_res->Fill(resonance_id, branch_id);
    h_ex_limit->Fill(ex_max_8Be, ex_8Be);
  }

  std::cout << "Entries: " << nentries << std::endl;
  std::cout << "165 resonance: " << n_165 << std::endl;
  std::cout << "675 resonance: " << n_675 << std::endl;
  std::cout << "alpha0 branch: " << n_alpha0 << std::endl;
  std::cout << "alpha1 branch: " << n_alpha1 << std::endl;
  std::cout << "Bad ex_8Be > ex_max_8Be: " << n_bad_ex_limit << std::endl;
  std::cout << "Bad cos range: " << n_bad_cos << std::endl;
  std::cout << "Bad branch_id: " << n_bad_branch << std::endl;
  std::cout << "Bad resonance_id: " << n_bad_resonance << std::endl;
  std::cout << "Non-finite diagnostic values: " << n_nonfinite << std::endl;
  if(nentries > n_nonfinite){
    std::cout << "ex_8Be range: [" << min_ex << ", " << max_ex << "] MeV" << std::endl;
    std::cout << "ex_max_8Be range: [" << min_ex_max << ", " << max_ex_max << "] MeV" << std::endl;
    std::cout << "e_alpha8Be range: [" << min_e_alpha8Be << ", " << max_e_alpha8Be << "] MeV" << std::endl;
  }

  TFile output(TString::Format("%s.root", output_prefix), "RECREATE");
  h_resonance->Write();
  h_branch->Write();
  h_ex->Write();
  h_eaa->Write();
  h_e_alpha8Be->Write();
  h_cos_theta->Write();
  h_cos_chi->Write();
  h_weight->Write();
  h_phi->Write();
  h_branch_res->Write();
  h_ex_limit->Write();
  output.Close();

  TCanvas* c = new TCanvas("c_h11b_reaction_check", "H11B reaction checks", 1600, 1200);
  c->Divide(3, 3);
  c->cd(1);
  h_resonance->Draw();
  c->cd(2);
  h_branch_res->Draw("COLZ TEXT");
  c->cd(3);
  h_weight->Draw();
  c->cd(4);
  h_ex->Draw();
  c->cd(5);
  h_eaa->Draw();
  c->cd(6);
  h_e_alpha8Be->Draw();
  c->cd(7);
  h_cos_theta->Draw();
  c->cd(8);
  h_cos_chi->Draw();
  c->cd(9);
  h_ex_limit->Draw("COLZ");
  c->SaveAs(TString::Format("%s.png", output_prefix));
  c->SaveAs(TString::Format("%s.pdf", output_prefix));
}
