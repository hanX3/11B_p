#include "TFile.h"
#include "TTree.h"

#include <algorithm>
#include <cstring>
#include <iostream>

void analyze_step_limit(const char* input)
{
  auto file = TFile::Open(input);
  auto tree = static_cast<TTree*>(file->Get("tr"));

  double pre_kine_energy = 0.0;
  double post_kine_energy = 0.0;
  double length = 0.0;
  char volume[128] = {};
  char particle[128] = {};

  tree->SetBranchAddress("pre_kine_energy", &pre_kine_energy);
  tree->SetBranchAddress("post_kine_energy", &post_kine_energy);
  tree->SetBranchAddress("length", &length);
  tree->SetBranchAddress("volume", volume);
  tree->SetBranchAddress("particle", particle);

  Long64_t n_target_proton = 0;
  double max_de = 0.0;
  double min_de = 1.0e99;
  double sum_de = 0.0;
  double max_length = 0.0;
  double first_pre = 0.0;
  double last_post = 0.0;
  double max_de_150_170 = 0.0;
  Long64_t n_150_170 = 0;
  double max_de_650_700 = 0.0;
  Long64_t n_650_700 = 0;

  const Long64_t nentries = tree->GetEntries();
  for(Long64_t i = 0; i < nentries; ++i){
    tree->GetEntry(i);
    if(std::strcmp(volume, "Target") != 0 || std::strcmp(particle, "proton") != 0) continue;

    const double de_keV = 1000.0 * (pre_kine_energy - post_kine_energy);
    if(n_target_proton == 0) first_pre = 1000.0 * pre_kine_energy;
    last_post = 1000.0 * post_kine_energy;
    max_de = std::max(max_de, de_keV);
    min_de = std::min(min_de, de_keV);
    sum_de += de_keV;
    max_length = std::max(max_length, length);
    if(1000.0 * pre_kine_energy >= 150.0 && 1000.0 * pre_kine_energy <= 170.0){
      max_de_150_170 = std::max(max_de_150_170, de_keV);
      ++n_150_170;
    }
    if(1000.0 * pre_kine_energy >= 650.0 && 1000.0 * pre_kine_energy <= 700.0){
      max_de_650_700 = std::max(max_de_650_700, de_keV);
      ++n_650_700;
    }
    ++n_target_proton;
  }

  std::cout << "target proton steps: " << n_target_proton << '\n';
  std::cout << "first pre energy: " << first_pre << " keV\n";
  std::cout << "last post energy: " << last_post << " keV\n";
  std::cout << "max step length: " << max_length / 1.0e-6 << " nm\n";
  std::cout << "max dE per step: " << max_de << " keV\n";
  std::cout << "min dE per step: " << min_de << " keV\n";
  std::cout << "avg dE per step: " << (n_target_proton > 0 ? sum_de / n_target_proton : 0.0) << " keV\n";
  std::cout << "steps per 5.3 keV width: " << (max_de > 0.0 ? 5.3 / max_de : 0.0) << '\n';
  std::cout << "steps per 300 keV width: " << (max_de > 0.0 ? 300.0 / max_de : 0.0) << '\n';
  std::cout << "150-170 keV window steps: " << n_150_170 << '\n';
  std::cout << "150-170 keV max dE per step: " << max_de_150_170 << " keV\n";
  std::cout << "150-170 keV steps per 5.3 keV width: " << (max_de_150_170 > 0.0 ? 5.3 / max_de_150_170 : 0.0)
            << '\n';
  std::cout << "650-700 keV window steps: " << n_650_700 << '\n';
  std::cout << "650-700 keV max dE per step: " << max_de_650_700 << " keV\n";
}
