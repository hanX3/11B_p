#ifndef DataStructure_H
#define DataStructure_H 1

#include <globals.hh>
#include <cstring>

//
struct H11BReactionData{
  G4long event;
  G4double e_alpha1;
  G4double e_alpha2;
  G4double e_alpha3;
  G4double theta_lab_alpha1;
  G4double theta_lab_alpha2;
  G4double theta_lab_alpha3;
  G4double phi_lab_alpha1;
  G4double phi_lab_alpha2;
  G4double phi_lab_alpha3;

  G4double ex_8Be; // exciation energy of 8Be
  G4double e_8Be;
  G4double theta_lab_8Be;
  G4double phi_lab_8Be;

  G4double x;
  G4double y;
  G4double z;

  // sequential1 / sequential2 / simultaneous
  // probability: sequential1 > sequential2 > simultaneous
  char reaction[128];

  //
  H11BReactionData() : 
    event(0),
    e_alpha1(0),
    e_alpha2(0),
    e_alpha3(0),
    theta_lab_alpha1(0),
    theta_lab_alpha2(0),
    theta_lab_alpha3(0),
    phi_lab_alpha1(0),
    phi_lab_alpha2(0),
    phi_lab_alpha3(0),
    ex_8Be(0),
    e_8Be(0),
    theta_lab_8Be(0),
    phi_lab_8Be(0),
    
    x(0.), y(0.), z(0.)
  {
    std::memset(reaction, '\0', sizeof(reaction));
  }

  //
  void Clear() {
    *this = H11BReactionData();
  }
};

//
struct EventData{
  G4long event;
  G4int ring;
  G4int sector;
  G4double e;
  G4double x;
  G4double y;
  G4double z;
  char detector[128];

  //
  EventData() : 
    event(0),
    ring(0),
    sector(0),
    e(0.), 
    x(0.), y(0.), z(0.)
  {
    std::memset(detector, '\0', sizeof(detector));
  }

  //
  void Clear() {
    *this = EventData();
  }
};

//
struct TrackData{
  G4long event;
  G4int track;
  G4double e;
  G4double x;
  G4double y;
  G4double z;
  G4double ts;
  G4double length;
  char volume[128];
  char particle[128];

  //
  TrackData() :
    event(-1), 
    track(-1),
    e(0.),
    x(0.), y(0.), z(0.),
    ts(0.), 
    length(0.)
  {
    std::memset(volume, '\0', sizeof(volume));
    std::memset(particle, '\0', sizeof(particle));
  }

  //
  void Clear() {
    *this = TrackData();
  }
};

//
struct StepData{
  G4long event;
  G4int track;
  G4double de;
  G4double pre_x;
  G4double pre_y;
  G4double pre_z;
  G4double pre_total_energy;
  G4double pre_kine_energy;
  G4double post_x;
  G4double post_y;
  G4double post_z;
  G4double post_total_energy;
  G4double post_kine_energy;
  G4double length;
  char volume[128];
  char particle[128];
  char process[128];

  //
  StepData() : 
    event(-1), 
    track(-1),
    de(0.), 
    pre_x(0.), pre_y(0.), pre_z(0.),
    pre_total_energy(0.), pre_kine_energy(0.),
    post_x(0.), post_y(0.), post_z(0.), 
    post_total_energy(0.), post_kine_energy(0.),
    length(0.)
  {
    std::memset(volume, '\0', sizeof(volume));
    std::memset(particle, '\0', sizeof(particle));
    std::memset(process, '\0', sizeof(process));
  }

  //
  void Clear() {
    *this = StepData();
  }
};
#endif
