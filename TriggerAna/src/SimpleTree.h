#ifndef _SIMPLE_TREE_H_
#define _SIMPLE_TREE_H_
#include <fun4all/SubsysReco.h>
#include <ktracker/FastTracklet.h>
#include <geom_svc/GeomSvc.h>
#include <interface_main/SQHit_v1.h>
#include "TreeData.h"

class TFile;
class TTree;
class SQEvent;
class SQMCEvent;
class SQTrackVector;
class SQDimuonVector;
class SQHitVector;
class GeomSvc;
class Tracklet;
class SQHit;

class SimpleTree: public SubsysReco
{
  
 private:
  int GetNodes(PHCompositeNode *topNode);
  void MakeTree();
  SQHit* findHit(const int detID, const int eleID);

  /// Detector Geometry
  GeomSvc *gs;

  /// Input
  SQMCEvent     * sqMCEvent;
  SQEvent       * sqEvent;
  SQTrackVector * sqTrackVector;
  TrackletVector* trackletVector;
  SQDimuonVector* sqDimuonVector;
  SQHitVector   * sqHitVector;

  /// Output
  TFile* file;
  TTree* tree;
  MCEventData dMCEvent;
  EventData dEvent;
  DimuonList lDimuon;
  TrackList lTrack;
  HitList lHit;
  TrackletList lTracklet;
 
  int nDimuons;
  int nHits;
  int nTracks;
  int nTracklets;

 public:
  SimpleTree(const std::string& name = "SimpleTree");
  virtual ~SimpleTree() {;}
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
  void SetOutput(const char* out_name);
};

#endif /* _SIMPLE_TREE_H_ */
