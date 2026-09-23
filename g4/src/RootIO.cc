#include "RootIO.hh"

#include <iostream>
#include <stdio.h>
#include <time.h>
#include <fstream>
#include <sstream>
#include "G4UnitsTable.hh"
#include "G4ThreeVector.hh"

//
RootIO::RootIO()
{
  // reaction
  reaction_data.Clear();

  // event
  event_data.Clear();

  // track
  track_data.Clear();

  // step
  step_data.Clear();
}

//
RootIO::~RootIO()
{
  // reaction
  if(reaction_file){
    delete reaction_file;
    reaction_file = nullptr;
  }

  // event
  if(event_file){
    delete event_file;
    event_file = nullptr;
  }

  // track
  if(track_file){
    delete track_file;
    track_file = nullptr;
  }

  // step
  if(step_file){
    delete step_file;
    step_file = nullptr;
  }
}

//
void RootIO::OpenReactionFile()
{
  time_t t;
  struct tm* tt;
  t=time(0);
  tt=localtime(&t);
  sprintf(file_name, "%d%02d%02d_%02dh%02dm%02ds", tt->tm_year+1900, tt->tm_mon+1, tt->tm_mday, tt->tm_hour, tt->tm_min, tt->tm_sec);
  G4cout << "\n----> Tree file is opened in " << file_name << G4endl;

  std::stringstream ss;
  ss.str("");
  ss << DATAPATH << "/reaction_" << file_name << ".root";

  reaction_file = new TFile(ss.str().c_str(), "RECREATE");
  //  check file_in
  if(!reaction_file){
    G4cout << " RootIO::" << " problem creating the ROOT TFile!!!" << G4endl;
    return ;
    }else{
      G4cout << " RootIO::" << " successful creating the " << ss.str().c_str() << "  !!!" << G4endl; 
    }

    //  creat tree
    reaction_tree = new TTree("tr", "reaction simulation data");

    reaction_tree->Branch("event", &reaction_data.event, "event/L");
    reaction_tree->Branch("e_alpha1", &reaction_data.e_alpha1, "e_alpha1/D");
    reaction_tree->Branch("e_alpha2", &reaction_data.e_alpha2, "e_alpha2/D");
    reaction_tree->Branch("e_alpha3", &reaction_data.e_alpha3, "e_alpha3/D");
    reaction_tree->Branch("theta_lab_alpha1", &reaction_data.theta_lab_alpha1, "theta_lab_alpha1/D");
    reaction_tree->Branch("theta_lab_alpha2", &reaction_data.theta_lab_alpha2, "theta_lab_alpha2/D");
    reaction_tree->Branch("theta_lab_alpha3", &reaction_data.theta_lab_alpha3, "theta_lab_alpha3/D");
    reaction_tree->Branch("e_3alpha_cm_alpha1", &reaction_data.e_3alpha_cm_alpha1, "e_3alpha_cm_alpha1/D");
    reaction_tree->Branch("e_3alpha_cm_alpha2", &reaction_data.e_3alpha_cm_alpha2, "e_3alpha_cm_alpha2/D");
    reaction_tree->Branch("e_3alpha_cm_alpha3", &reaction_data.e_3alpha_cm_alpha3, "e_3alpha_cm_alpha3/D");
    reaction_tree->Branch("ex_8Be", &reaction_data.ex_8Be, "ex_8Be/D");
    reaction_tree->Branch("e_8Be", &reaction_data.e_8Be, "e_8Be/D");
    reaction_tree->Branch("theta_lab_8Be", &reaction_data.theta_lab_8Be, "theta_lab_8Be/D");
    reaction_tree->Branch("phi_lab_8Be", &reaction_data.phi_lab_8Be, "phi_lab_8Be/D");
    reaction_tree->Branch("x", &reaction_data.x, "x/D");
    reaction_tree->Branch("y", &reaction_data.y, "y/D");
    reaction_tree->Branch("z", &reaction_data.z, "z/D");
    reaction_tree->Branch("reaction", reaction_data.reaction, "reaction/C");

    if(!reaction_tree){
      G4cout << "\n can't create tree" << G4endl;
      return ;
    }

    G4cout << "\n----> Tree file is opened in " << ss.str() << G4endl;
}

//
void RootIO::FillReactionTree(H11BReactionData& data)
{
  reaction_data = data; 

  reaction_tree->Fill();
}

//
void RootIO::CloseReactionFile()
{
  if(!reaction_file){
    return;
  }

  reaction_file->cd();
  reaction_tree->Write();
  reaction_file->Close();

  G4cout << "\n----> reaction tree is saved.\n" << G4endl;
}

//
void RootIO::OpenEventFile()
{
  time_t t;
  struct tm* tt;
  t=time(0);
  tt=localtime(&t);
  sprintf(file_name, "%d%02d%02d_%02dh%02dm%02ds", tt->tm_year+1900, tt->tm_mon+1, tt->tm_mday, tt->tm_hour, tt->tm_min, tt->tm_sec);
  G4cout << "\n----> Tree file is opened in " << file_name << G4endl;

  std::stringstream ss;
  ss.str("");
  ss << DATAPATH << "/event_" << file_name << ".root";

  event_file = new TFile(ss.str().c_str(), "RECREATE");
  //  check file_in
  if(!event_file){
    G4cout << " RootIO::" << " problem creating the ROOT TFile!!!" << G4endl;
    return ;
    }else{
      G4cout << " RootIO::" << " successful creating the " << ss.str().c_str() << "  !!!" << G4endl; 
    }

    //  creat tree
    event_tree = new TTree("tr", "event simulation data");

    event_tree->Branch("event", &event_data.event, "event/L");
    event_tree->Branch("ring", &event_data.ring, "ring/I");
    event_tree->Branch("sector", &event_data.sector, "sector/I");
    event_tree->Branch("e", &event_data.e, "e/D");
    event_tree->Branch("x", &event_data.x, "x/D");
    event_tree->Branch("y", &event_data.y, "y/D");
    event_tree->Branch("z", &event_data.z, "z/D");
    event_tree->Branch("detector", event_data.detector, "detector/C");

    if(!event_tree){
      G4cout << "\n can't create tree" << G4endl;
      return ;
    }

    G4cout << "\n----> Tree file is opened in " << ss.str() << G4endl;
}

//
void RootIO::FillEventTree(EventData& data)
{
  event_data = data;
  event_tree->Fill();
}

//
void RootIO::CloseEventFile()
{
  if(!event_file){
    return;
  }

  event_file->cd();
  event_tree->Write();
  event_file->Close();

  G4cout << "\n----> event tree is saved.\n" << G4endl;
}

//
void RootIO::OpenTrackFile()
{
  time_t t;
  struct tm* tt;
  t=time(0);
  tt=localtime(&t);
  sprintf(file_name, "%d%02d%02d_%02dh%02dm%02ds", tt->tm_year+1900, tt->tm_mon+1, tt->tm_mday, tt->tm_hour, tt->tm_min, tt->tm_sec);
  G4cout << "\n----> Tree file is opened in " << file_name << G4endl;
  std::stringstream ss;
  ss.str("");
  ss << DATAPATH << "/track_" << file_name << ".root";
  track_file = new TFile(ss.str().c_str(), "RECREATE");
  //  check file_in
  if(!track_file){
    G4cout << " RootIO::" << " problem creating the ROOT TFile!!!" << G4endl;
    return ;
  }else{
    G4cout << " RootIO::" << " successful creating the " << ss.str().c_str() << "  !!!" << G4endl; 
  }
  //  creat tree
  track_tree = new TTree("tr", "track simulation data");
  track_tree->Branch("event", &track_data.event, "event/L");
  track_tree->Branch("track", &track_data.track, "track/I");
  track_tree->Branch("e", &track_data.e, "e/D");
  track_tree->Branch("x", &track_data.x, "x/D");
  track_tree->Branch("y", &track_data.y, "y/D");
  track_tree->Branch("z", &track_data.z, "z/D");
  track_tree->Branch("ts", &track_data.ts, "ts/D");
  track_tree->Branch("length", &track_data.length, "length/D");
  track_tree->Branch("volume", track_data.volume, "volume/C");
  track_tree->Branch("particle", track_data.particle, "particle/C");

  if(!track_tree){
    G4cout << "\n can't create tree" << G4endl;
    return ;
  }
  G4cout << "\n----> Tree file is opened in " << ss.str() << G4endl;
}

//
void RootIO::FillTrackTree(TrackData& data)
{
  track_data = data;
  track_tree->Fill();
}

//
void RootIO::CloseTrackFile()
{
  if(!track_file){
    return;
  }

  track_file->cd();
  track_tree->Write();
  track_file->Close();

  G4cout << "\n----> track tree is saved.\n" << G4endl;
}

//
void RootIO::OpenStepFile()
{
  time_t t;
  struct tm* tt;
  t=time(0);
  tt=localtime(&t);
  sprintf(file_name, "%d%02d%02d_%02dh%02dm%02ds", tt->tm_year+1900, tt->tm_mon+1, tt->tm_mday, tt->tm_hour, tt->tm_min, tt->tm_sec);
  G4cout << "\n----> Tree file is opened in " << file_name << G4endl;
  std::stringstream ss;
  ss.str("");
  ss << DATAPATH << "/step_" << file_name << ".root";
  step_file = new TFile(ss.str().c_str(), "RECREATE");
  //  check file_in
  if(!step_file){
    G4cout << " RootIO::" << " problem creating the ROOT TFile!!!" << G4endl;
    return ;
  }else{
    G4cout << " RootIO::" << " successful creating the " << ss.str().c_str() << "  !!!" << G4endl; 
  }
  //  creat tree
  step_tree = new TTree("tr", "step simulation data");
  step_tree->Branch("event", &step_data.event, "event/L");
  step_tree->Branch("track", &step_data.track, "track/I");
  step_tree->Branch("de", &step_data.de, "de/D");
  step_tree->Branch("pre_x", &step_data.pre_x, "pre_x/D");
  step_tree->Branch("pre_y", &step_data.pre_y, "pre_y/D");
  step_tree->Branch("pre_z", &step_data.pre_z, "pre_z/D");
  step_tree->Branch("pre_total_energy", &step_data.pre_total_energy, "pre_total_energy/D");
  step_tree->Branch("pre_kine_energy", &step_data.pre_kine_energy, "pre_kine_energy/D");
  step_tree->Branch("post_x", &step_data.post_x, "post_x/D");
  step_tree->Branch("post_y", &step_data.post_y, "post_y/D");
  step_tree->Branch("post_z", &step_data.post_z, "post_z/D");
  step_tree->Branch("post_total_energy", &step_data.post_total_energy, "post_total_energy/D");
  step_tree->Branch("post_kine_energy", &step_data.post_kine_energy, "post_kine_energy/D");
  step_tree->Branch("length", &step_data.length, "length/D");
  step_tree->Branch("volume", step_data.volume, "volume/C");
  step_tree->Branch("particle", step_data.particle, "particle/C");
  step_tree->Branch("process", step_data.process, "process/C");
  
  if(!step_tree){
    G4cout << "\n can't create tree" << G4endl;
    return ;
  }
  G4cout << "\n----> Tree file is opened in " << ss.str() << G4endl;
}

//
void RootIO::FillStepTree(StepData& data)
{
  step_data = data;
  step_tree->Fill();
}

//
void RootIO::CloseStepFile()
{
  if(!step_file){
    return;
  }

  step_file->cd();
  step_tree->Write();
  step_file->Close();

  G4cout << "\n----> step tree is saved.\n" << G4endl;
}

