//
void calculate_gamma_angular_distribution(double beta, int type = 0)
{
  gStyle->SetOptStat(0);

  //
  TGraph *gr = new TGraph();

  double sum_dist = 0.;
  vector<double> v_sum_dist_now;
  //distribution 0: uniform; 2200: J,M 2,2-> 0,0; 2100: J,M 2,1->0,0 2000: J,M 2,0 -> 0,0
  // The zet axis is the beam axis
  for(int m2=(int)(0.5*100);m2<(int)(179.5*100+1);m2++){         
    double angle  = m2/100.;
    if(type == 2000)  
      sum_dist += 1.0*sin(angle/180.*TMath::Pi())
        *( sin(angle/180.*TMath::Pi())*sin(angle/180.*TMath::Pi())*cos(angle/180.*TMath::Pi())*cos(angle/180.*TMath::Pi()));
    else if(type == 2100)
      sum_dist += 1.0*sin(angle/180.*TMath::Pi())
        *(1 - 3*cos(angle/180.*TMath::Pi())*cos(angle/180.*TMath::Pi()) + 4*cos(angle/180.*TMath::Pi())*cos(angle/180.*TMath::Pi())*cos(angle/180.*TMath::Pi())*cos(angle/180.*TMath::Pi()));
    else if (type == 2200)
      sum_dist += 1.0*sin(angle/180.*TMath::Pi())
        *(1 - cos(angle/180.*TMath::Pi())*cos(angle/180.*TMath::Pi())*cos(angle/180.*TMath::Pi())*cos(angle/180.*TMath::Pi()));
    // Making it uniform:
    else sum_dist += 1.0*sin(angle/180.*TMath::Pi());

    v_sum_dist_now.push_back(sum_dist);

    gr->SetPoint(m2, angle, sum_dist);
  }

  auto c1 = new TCanvas();
  c1->cd();
  gr->Draw();


  //
  TH1D *h1 = new TH1D("h1", "", 180, 0., 180.);
  h1->GetXaxis()->SetTitle("#theta (degree)");
  double theta_cm, theta_lab;
  double rand;

  double x, y, z;
  for(int i=0;i<1000000;i++){
    rand = gRandom->Uniform(0, sum_dist);
    // cout << "sum_dist " << rand << endl;

    for(int j=1;j<v_sum_dist_now.size();j++){
      if(rand>=v_sum_dist_now[j-1] && rand<=v_sum_dist_now[j]){
        theta_cm = j/100.0;
        break;
      }
    }

    theta_cm = theta_cm/180.*TMath::Pi();

    theta_lab = acos((cos(theta_cm)+beta)/(1.0+beta*cos(theta_cm)));
    h1->Fill(theta_lab/TMath::Pi()*180.);
  }

  auto c2 = new TCanvas();
  c2->cd();
  h1->Draw();
}

// we also can get beta value using lise++ code
void calculate_beta(int mass, double energy)
{
  energy *= mass;
  double gamma = (931.494+energy/mass)/931.494;
  double beta = sqrt(1.0 - 1.0/gamma/gamma);

  cout << "gamma " << gamma << endl;
  cout << "beta " << beta << endl;
}