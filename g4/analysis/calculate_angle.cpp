//
void calculate_angle()
{
  TFile *fi = TFile::Open("../data/event_20240409_12h31m35s.root");
  if(fi->IsZombie()){
    cout << "can not open the root file" << endl;
    return;
  }

  Long64_t event;
  Int_t ring;
  Int_t sector;
  Double_t e;
  Double_t x;
  Double_t y;
  Double_t z;

  //
  TTree *tr_in = (TTree*)fi->Get("tr");
  tr_in->SetBranchAddress("event", &event);
  tr_in->SetBranchAddress("ring", &ring);
  tr_in->SetBranchAddress("sector", &sector);
  tr_in->SetBranchAddress("e", &e);
  tr_in->SetBranchAddress("x", &x);
  tr_in->SetBranchAddress("y", &y);
  tr_in->SetBranchAddress("z", &z);

  //
  TFile *fo = new TFile("angle.root", "recreate");
  TTree *tr_out = new TTree("tr", "ring_theta_phi");
  Double_t theta;
  Double_t phi; 
  tr_out->Branch("ring", &ring, "ring/I");
  tr_out->Branch("sector", &sector, "sector/I");
  tr_out->Branch("theta", &theta, "theta/D");
  tr_out->Branch("phi", &phi, "phi/D");

  TVector3 v;
  //
  Long64_t entries = tr_in->GetEntries();
  for(int i=0;i<entries;i++){
    tr_in->GetEntry(i);
    v.SetXYZ(x, y, z);
    theta = v.Theta();
    phi = v.Phi();

    // cout << ring << " " << sector << " " << x << " " << y << " " << z << endl;
    // cout << theta << " " << phi << endl;  

    tr_out->Fill();
  }

  fo->cd();
  tr_out->Write();

  fi->Close();
  fo->Close();
}