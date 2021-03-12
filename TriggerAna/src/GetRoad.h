/*

FPGA trigger relavant functions

Detector names are hardcoded.
Please confirm with GeomSvc functions if there's change.

Top, Bottom: 0, 1 for array order; +1, -1 as a value
Pos, Neg: 0, 1 for array order; +1, -1 as a value

Minjung Kim, kminjung@umich.edu

*/

#ifndef _FPGA_TRIG_H_
#define _FPGA_TRIG_H_

#include <iostream>
#include <vector>
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <geom_svc/GeomSvc.h>

class GetRoad
{
 public:
  
  GetRoad(const char* out_name, const int nevt);
  void get_mc_roads(const int idata, const int nevt);
  void get_old_roads(const int idata);
  
  int channels_to_roadID(const std::vector<int> detectorIDs, const std::vector<int> elementIDs);
  bool roadID_to_channels(const int roadID,  std::vector<int> &detectorIDs,  std::vector<int> &elementIDs);
  void reset_var();
  void get_road_id(const int i);

 private:
  GeomSvc* gs;

  TTree* tr[3];
  TTree* tr_road[3];
  
  // dimuon truth
  int nDimuon;
  float mass;
  float xf;
  float x1;
  float x2;
  float px[2];
  int track_id[2];
  int charge[2];
  // digitized hit
  int nHit;
  int tb[2]; // 0 for not-all-4-stations hit at top/bottom
  int unique_id[2][4]; // unique ID = 1000*detectorID + elementID
  int detector_id[2][4];
  float detector_zpos[2][4];
  int element_id[2][4];
  float element_xpos[2][4];
  // road
  int road_id[2]; // 0 for out-of-accentance track
  int road_freq[2]; // roads frequency for roadset67, 1 for new roadset in this step

  // to save unique road [+/- muon]
  std::vector<int> comp_id[2];
  std::vector<int>::iterator it[2];

  int u_charge;
  int u_tb;

  int u_rid;
  int u_uid[4];
  int u_eid[4];
  float u_z[4];
  float u_x[4];
  
  int u_freq;
  float u_px_min;
  float u_px_max;
  float u_px_avg;
  float u_px_rms;

  TFile* fp_mc[2];

};


#endif /* _FPGA_TRIG_H_ */
