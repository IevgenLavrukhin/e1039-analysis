/**
 * \class AnaTrkQA
 * \Analysis module for track QA
 * \author  Abinash Pun
 *
 * Created: 07-05-2020
 */

#ifndef _H_AnaTrkQA_H_
#define _H_AnaTrkQA_H_

// ROOT
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TMatrixD.h>
#include <TVector3.h>

// Fun4All includes
#include <fun4all/SubsysReco.h>

// STL includes
#include <vector>
#include <string>
#include <iostream>
#include <list>
#include <map>

class TVector3;
class TLorentzVector;

class SQRun;
class SQSpillMap;
class SQEvent;
class SQHitMap;
class SQHitVector;
class SQHit;

class PHG4TruthInfoContainer;
class PHG4HitContainer;
class SRecEvent;
class SRecTrack;
class GeomSvc;

class TFile;
class TTree;

class AnaTrkQA: public SubsysReco {

 public:

  AnaTrkQA(const std::string &name = "AnaTrkQA.root");
  virtual ~AnaTrkQA() {
  }

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  int InitEvalTree();
  int ResetEvalVars();

  const std::string& get_hit_container_choice() const {
    return _hit_container_type;
  }

  void set_hit_container_choice(const std::string& hitContainerChoice) {
    _hit_container_type = hitContainerChoice;
  }

  const std::string& get_out_name() const {
    return _out_name;
  }

  void set_out_name(const std::string& outName) {
    _out_name = outName;
  }

 private:

  int GetNodes(PHCompositeNode *topNode);

  int TruthRecoEval(PHCompositeNode *topNode);
  
  bool FindG4HitAtStation(const int trk_id, const PHG4HitContainer* g4hc, TVector3* pos, TLorentzVector* mom);
  int  FindCommonHitIDs(std::vector<int>& hitidvec1, std::vector<int>& hitidvec2);
  SRecTrack* FindBestMomRecTrack(SRecEvent *recEvent, const float true_P); 
  bool FindG4HitAtHodo(const int trk_id, const PHG4HitContainer* g4hc);
  bool FindG4HitAtProp(const int trk_id, const PHG4HitContainer* g4hc);
 
  std::string _hit_container_type;

  size_t _event;
  SQRun* _run_header;
  SQSpillMap * _spill_map;

  SQEvent * _event_header;
  SQHitMap *_hit_map;
  SQHitVector *_hit_vector;

  PHG4TruthInfoContainer* _truth;
  SRecEvent* _recEvent;

  PHG4HitContainer *g4hc_d1x;
  PHG4HitContainer *g4hc_d2xp;
  PHG4HitContainer *g4hc_d3px;
  PHG4HitContainer *g4hc_d3mx;

  PHG4HitContainer *g4hc_h1t;
  PHG4HitContainer *g4hc_h1b;
  PHG4HitContainer *g4hc_h2t;
  PHG4HitContainer *g4hc_h2b;
  PHG4HitContainer *g4hc_h3t;
  PHG4HitContainer *g4hc_h3b;
  PHG4HitContainer *g4hc_h4t;
  PHG4HitContainer *g4hc_h4b;

  PHG4HitContainer *g4hc_p1y1;
  PHG4HitContainer *g4hc_p1y2;
  PHG4HitContainer *g4hc_p1x1;
  PHG4HitContainer *g4hc_p1x2;
  PHG4HitContainer *g4hc_p2x1;
  PHG4HitContainer *g4hc_p2x2;
  PHG4HitContainer *g4hc_p2y1;
  PHG4HitContainer *g4hc_p2y2;

  std::string _out_name;
	
  TTree* _tout_reco;
  TTree* _qa_tree;
  TFile *file;

  int run_id;
  int spill_id;
  float target_pos;
  int event_id;
  int krecstat;
  int kalman_stat;
  unsigned short emu_trigger;

  int n_hits;
  int nhits_st1[100];
  int nhits_st2[100];
  int nhits_st3[100];
  int hit_id[100];
  int detector_id[100];
  int element_id[100];
  int hodo_mask[100];
  float drift_distance[100];
  float pos[100];
  float detector_z[100];

  float truth_x[100];
  float truth_y[100];
  float truth_z[100];
  float truth_pos[100];

  int n_tracks;
  int n_recTracks;
  int rec_id[100];
  int par_id[100];
  int pid[100];

  float gvx[100];
  float gvy[100];
  float gvz[100];
  float gpx[100];
  float gpy[100];
  float gpz[100];
  float gx_st1[100];
  float gy_st1[100];
  float gz_st1[100];
        
  float ac_gpx[100];
  float ac_gpy[100];
  float ac_gpz[100];  
   
  float sq_pos_st1[100];
  float sq_drift_st1[100];
  float sq_decID[100];        
  float sq_x_st1[100];
  float sq_y_st1[100];
  float sq_z_st1[100];
  float sq_px_st1[100];
  float sq_py_st1[100];
  float sq_pz_st1[100];


  float sq_pos_st2[100];
  float sq_drift_st2[100];
  float sq_x_st2[100];
  float sq_y_st2[100];
  float sq_z_st2[100];
  float sq_px_st2[100];
  float sq_py_st2[100];
  float sq_pz_st2[100];


  float sq_pos_st3[100];
  float sq_drift_st3[100];
  float sq_x_st3[100];
  float sq_y_st3[100];
  float sq_z_st3[100];
  float sq_px_st3[100];
  float sq_py_st3[100];
  float sq_pz_st3[100];


  float chisq_st1[100];
  float prob_st1[100];
  float quality[100];
  float fit_time[100];
  
  float pull_q2p_st1[100];  
  float pull_q2p_st2[100];
  float pull_q2p_st3[100];


  float rec_drift_st1[100];
  float rec_px_st1[100];
  float rec_py_st1[100];
  float rec_pz_st1[100];
  float rec_p_st1[100];
  float rec_x_st1[100];
  float rec_y_st1[100];  	
  float rec_z_st1[100];

  float rec_drift_st2[100];
  float rec_px_st2[100];
  float rec_py_st2[100];
  float rec_pz_st2[100];
  float rec_p_st2[100];
  float rec_x_st2[100];
  float rec_y_st2[100];
  float rec_z_st2[100];


  float rec_drift_st3[100];
  float rec_px_st3[100];
  float rec_py_st3[100];
  float rec_pz_st3[100];
  float rec_p_st3[100];
  float rec_x_st3[100];
  float rec_y_st3[100];
  float rec_z_st3[100];
     
  float gpx_st1[100];
  float gpy_st1[100];
  float gpz_st1[100];
  float gpt[100];
  float geta[100];
  float gphi[100];
  int gnhits[100];
  int gndc[100];
  int gnhodo[100];
  int gnprop[100];
  int gndp[100];
  int ntruhits[100];
  int nhits[100];
  int charge[100];
  float rec_vx[100];
  float rec_vy[100];
  float rec_vz[100];
  float rec_px[100];
  float rec_py[100];
  float rec_pz[100];
  float rec_pt[100];
  float rec_eta[100];
  float rec_phi[100];
  float pull_state00[100];

  int gelmid[1000][128];



  GeomSvc *p_geomSvc;
};


#endif /* _H_AnaTrkQA_H_ */
