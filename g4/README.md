# HFRS_gamma

## Requirements
+ HFRS has only been tested on GEANT4 v11.1.3.
+ ROOT installation with c++17 standard is required.

## Running
+ in HFRS_gamma directory
```
mkdir build
cd build
cmake ..
make -j6
./HFRS macro/gamma.mac (for example)
```

## Geometry
+ 18 LaBr3 detectors. ring id = 1, sector id = 0,1,2,...5. ring id = 2, sector id = 0,1,2,...,11.
+ 24 GAGG detectors. ring id = 3, sector id = 0,1,2,...11. ring id = 4, sector id = 0,1,2,...,5. ring id = 5, sector id = 0,1,2,...,5. Notice: GAGG detectors are in the same theta angle, different ring id means different type of detectors.
+ 1024 CsI detectors. ring id 9,10,...,24. sector id 0,1,2,...,47.

## Functionality
+ save simulation rootfile in HFRS_gamma/data directory. A tree with branch event(id of event), ring(ring id), sector(sector id), e(energy), x,y,z(the first hit position in the sensitive detector) have been creat.
+ in analysis directory
  + calculate_angle.cpp: calculate theta and phi angle of each event.