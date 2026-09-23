#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TString.h"

#include <algorithm>
#include <cstring>
#include <vector>

void draw_dalitz(const char* input = "../data/reaction_all.root", bool sortEnergy = false, const char* reactionCut = "")
{
  TChain tr("tr");
  tr.Add(input);

  double e1 = 0.0;
  double e2 = 0.0;
  double e3 = 0.0;
  char reaction[128];

  tr.SetBranchAddress("e_3alpha_cm_alpha1", &e1);
  tr.SetBranchAddress("e_3alpha_cm_alpha2", &e2);
  tr.SetBranchAddress("e_3alpha_cm_alpha3", &e3);
  tr.SetBranchAddress("reaction", reaction);

  TH2D* h = new TH2D("h_dalitz", "p+^{11}B #rightarrow 3#alpha Dalitz plot;X;Y", 500, -1.9, 1.9, 500, -1.2, 2.2);

  const Long64_t nentries = tr.GetEntries();

  for(Long64_t i = 0; i < nentries; ++i){
    tr.GetEntry(i);

    if(std::strlen(reactionCut) > 0){
      if(TString(reaction) != TString(reactionCut)) continue;
    }

    double E1 = e1;
    double E2 = e2;
    double E3 = e3;

    if(sortEnergy){
      std::vector<double> E = {E1, E2, E3};
      std::sort(E.begin(), E.end(), std::greater<double>());

      E1 = E[0];
      E2 = E[1];
      E3 = E[2];
    }

    const double Esum = E1 + E2 + E3;
    if(Esum <= 0.0) continue;

    const double eps1 = E1 / Esum;
    const double eps2 = E2 / Esum;
    const double eps3 = E3 / Esum;

    const double X = TMath::Sqrt(3.0) * (eps2 - eps3);
    const double Y = 2.0 * eps1 - eps2 - eps3;

    h->Fill(X, Y);
  }

  TCanvas* c = new TCanvas("c_dalitz", "Dalitz plot", 900, 800);
  h->Draw("COLZ");

  c->SaveAs("dalitz.png");
  c->SaveAs("dalitz.pdf");
}
