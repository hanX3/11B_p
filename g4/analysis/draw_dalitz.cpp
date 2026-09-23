#include "TChain.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TString.h"
#include "TStyle.h"

#include <cmath>
#include <cstring>
#include <iostream>

void draw_dalitz(const char* input = "../data/reaction_all.root",
                    const char* reactionFilter = "p11B_675_sequential1",
                    bool symmetrize = true,
                    double minEx8Be = -1.0,
                    const char* outputPrefix = "dalitz_cm")
{
    TChain tr("tr");
    tr.Add(input);

    double e1 = 0.0;
    double e2 = 0.0;
    double e3 = 0.0;
    double ex8Be = 0.0;
    char reaction[256];

    tr.SetBranchAddress("e_3alpha_cm_alpha1", &e1);
    tr.SetBranchAddress("e_3alpha_cm_alpha2", &e2);
    tr.SetBranchAddress("e_3alpha_cm_alpha3", &e3);

    const bool hasReaction = (tr.GetBranch("reaction") != nullptr);
    const bool hasEx8Be = (tr.GetBranch("ex_8Be") != nullptr);

    if (hasReaction) {
        tr.SetBranchAddress("reaction", reaction);
    }

    if (hasEx8Be) {
        tr.SetBranchAddress("ex_8Be", &ex8Be);
    }

    TH2D* h = new TH2D(
        "h_dalitz_cm",
        "p+^{11}B #rightarrow 3#alpha Dalitz plot;X;Y",
        500, -1.05, 1.05,
        500, -1.05, 1.05
    );

    auto FillDalitz = [&](double E1, double E2, double E3)
    {
        const double Esum = E1 + E2 + E3;
        if (Esum <= 0.0) return;

        const double x = std::sqrt(3.0) * (E2 - E3) / Esum;
        const double y = (2.0 * E1 - E2 - E3) / Esum;

        h->Fill(x, y);
    };

    const Long64_t nentries = tr.GetEntries();
    Long64_t nselected = 0;

    for (Long64_t i = 0; i < nentries; ++i) {
        tr.GetEntry(i);

        if (hasReaction && std::strlen(reactionFilter) > 0) {
            if (!TString(reaction).Contains(reactionFilter)) continue;
        }

        if (hasEx8Be && minEx8Be >= 0.0) {
            if (ex8Be < minEx8Be) continue;
        }

        if (e1 <= 0.0 || e2 <= 0.0 || e3 <= 0.0) continue;

        nselected++;

        if (symmetrize) {
            FillDalitz(e1, e2, e3);
            FillDalitz(e2, e3, e1);
            FillDalitz(e3, e1, e2);
        }
        else {
            FillDalitz(e1, e2, e3);
        }
    }

    gStyle->SetOptStat(0);

    TCanvas* c = new TCanvas("c_dalitz_cm", "Dalitz plot in 3-alpha CM", 800, 800);
    h->Draw("COLZ");

    TString pngName = TString::Format("%s.png", outputPrefix);
    TString pdfName = TString::Format("%s.pdf", outputPrefix);

    c->SaveAs(pngName);
    c->SaveAs(pdfName);

    std::cout << "Input file      : " << input << std::endl;
    std::cout << "Reaction filter : " << reactionFilter << std::endl;
    std::cout << "minEx8Be        : " << minEx8Be << std::endl;
    std::cout << "Total entries   : " << nentries << std::endl;
    std::cout << "Selected events : " << nselected << std::endl;
    std::cout << "Histogram fills : " << h->GetEntries() << std::endl;
    std::cout << "Output          : " << pngName << ", " << pdfName << std::endl;
}
