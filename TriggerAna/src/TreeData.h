#ifndef _MJ_TREE_DATA_H__ 
#define _MJ_TREE_DATA_H__
#include <TObject.h>
#include <TLorentzVector.h>


class MCEventData : public TObject
{
 public:
  /// MC Event Info
  int process_id; // 223
  int particle_id[4]; // 2->2
  TLorentzVector particle_mom[4];

  MCEventData();
  virtual ~MCEventData() {;}

  ClassDef(MCEventData, 1);
};

class EventData : public TObject
{
 public:
  int event_id;
  int trig_bits; // trigger bits

  EventData();
  virtual ~EventData() {;}

  ClassDef(EventData, 1);
};

class DimuonData : public TObject
{
 public:
  int dim_id;
  int pdg_id;
  double mass;
  double x1;
  double x2;
  double xf;
  int track_id_pos;
  int track_id_neg;
  TLorentzVector mom;
  TLorentzVector mom_pos;
  TLorentzVector mom_neg;

  DimuonData();
  virtual ~DimuonData() {;}

  ClassDef(DimuonData, 1);
};

class TrackData : public TObject
{
 public:
  int track_id;
  int charge;
  float zvtx;

  TrackData();
  virtual ~TrackData() {;}
  
  ClassDef(TrackData, 1);
};

class HitData : public TObject
{
 public:
  int hit_id;
  int track_id;
  int detector_id;
  int element_id;
  float tdc_time;
  bool in_time;

  HitData();
  virtual ~HitData() {;}

  ClassDef(HitData, 1);
};

class TrackletData : public TObject
{
 public:
  int n_hits;
  float chisq;
  int charge;
  int detector_id[24];
  float detector_zpos[24];
  float x_exp[24]; // track x
  float y_exp[24]; // track y
  bool in_plane[24]; // track at detector acceptance
  int element_id_exp[24]; // expected hit element
  int hit_id[24];
  int track_id[24];
  int element_id_pos[24]; // expected hit element position
  int element_id_closest[24]; // exp, nearest, or -1
	     
  TrackletData();
  virtual ~TrackletData() {;}
  
  ClassDef(TrackletData, 1);
};

typedef std::vector<DimuonData>   DimuonList;
typedef std::vector<TrackData>    TrackList;
typedef std::vector<HitData>      HitList;
typedef std::vector<TrackletData> TrackletList;
#endif /* _MJ_TREE_DATA_H__ */
