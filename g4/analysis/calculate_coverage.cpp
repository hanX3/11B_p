
//
void calculate_coverage(double beta)
{
  TLorentzVector lv_t(0.0, 0.0, 0.0, 0.938);
  TLorentzVector lv_p(0.0, 0.0, beta, 0.938+beta); // 700 MeV/u projectiole, for example proton
  TLorentzVector lv_add = lv_t + lv_p;
  Double_t masses[3] = {0.938, 0.938, 0.}; // proton and gamma-ray

  TGenPhaseSpace event;
  event.SetDecay(lv_add, 3, masses);

  TH1D *h1 = new TH1D("h1", TString::Format("beta=%.2f",beta).Data(), 180, 0, 180);
  for(int i=0;i<1000000;i++){
    Double_t weight = event.Generate();
    // cout << weight << endl;

    TLorentzVector *lv_t2 = event.GetDecay(0);
    TLorentzVector *lv_p2 = event.GetDecay(1);
    TLorentzVector *lv_g = event.GetDecay(2);

    TLorentzVector lv_g2 = *lv_p2 + *lv_g;
    h1->Fill(lv_g2.Theta()/TMath::Pi()*180.);
  }

  h1->Draw();
}