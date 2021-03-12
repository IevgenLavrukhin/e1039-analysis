#include <iomanip>
#include <TFile.h>
#include <TTree.h>
#include <interface_main/SQEvent.h>
#include <interface_main/SQMCEvent.h>
#include <interface_main/SQTrackVector.h>
#include <interface_main/SQDimuonVector.h>
#include <interface_main/SQHitVector.h>
#include <ktracker/FastTracklet.h>
#include <geom_svc/GeomSvc.h>
#include <ktracker/SRecEvent.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include "SimpleTree.h"

using namespace std;

// XT, XB, YL, YR, DP1, DP2
const int dID[24] =
  {32, 38, 40, 46, 31, 37, 39, 45,
   33, 35, 41, 43, 34, 36, 42, 44,
   55, 56, 57, 58, 59, 60, 61, 62};

SimpleTree::SimpleTree(const std::string& name)
  : SubsysReco(name),
    gs(GeomSvc::instance())
{}


int SimpleTree::Init(PHCompositeNode* topNode)
{
  cout<<"SimpleTree::Init()"<<endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int SimpleTree::InitRun(PHCompositeNode* topNode)
{
  cout<<"SimpleTree::InitRun() - GetNodes()"<<endl;
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;
    
  cout<<"SimpleTree::InitRun() - MakeTree()"<<endl;
  MakeTree();
  cout<<"SimpleTree::IninRun() Done"<<endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int SimpleTree::process_event(PHCompositeNode* topNode)
{
  static unsigned int n_evt = 0;
  if    (++n_evt % 100000 == 0) cout << n_evt << endl;
  else if (n_evt %  10000 == 0) cout << " . " << flush;


  /// Get MC Event Info 
  //cout<<"MC"<<endl;
  dMCEvent.process_id = sqMCEvent->get_process_id();
  for(int ii=0; ii<4; ii++)
    {
      dMCEvent.particle_id[ii] = sqMCEvent->get_particle_id(ii);
      dMCEvent.particle_mom[ii] = sqMCEvent->get_particle_momentum(ii);
    }


  /// Get DAQ Event Info
  //cout<<"Event"<<endl;
  dEvent.event_id   = sqEvent->get_event_id();
  dEvent.trig_bits  = sqEvent->get_trigger();


  /// Dimuon info
  //cout<<"Dimuon"<<endl;
  lDimuon.clear();
  nDimuons = sqDimuonVector->size();

  for (unsigned int ii = 0; ii < nDimuons; ii++) 
    {
      SQDimuon* dim = sqDimuonVector->at(ii);
      DimuonData dd;

      dd.dim_id       = dim->get_dimuon_id();
      dd.pdg_id       = dim->get_pdg_id();
      dd.mass         = dim->get_mass();
      dd.x1           = dim->get_x1();
      dd.x2           = dim->get_x2();
      dd.xf           = dim->get_xf();
      dd.track_id_pos = dim->get_track_id_pos();
      dd.track_id_neg = dim->get_track_id_neg();
      dd.mom          = dim->get_mom();
      dd.mom_pos      = dim->get_mom_pos();
      dd.mom_neg      = dim->get_mom_neg();
      
      lDimuon.push_back(dd);
    }

  
  /// Track info
  //cout<<"Track"<<endl;
  lTrack.clear();
  nTracks = sqTrackVector->size();
  
  for (unsigned int ii = 0; ii < nTracks; ii++) 
    {
      SQTrack* trk = sqTrackVector->at(ii); 
      TrackData dd;

      dd.track_id     = trk->get_track_id();
      dd.charge = trk->get_charge();
      dd.zvtx   = trk->get_pos_vtx().Z();

      lTrack.push_back(dd);
     }


  /// Hit info
  //cout<<"Hit"<<endl;
  lHit.clear();
  nHits = sqHitVector->size();
  
  for (unsigned int ii = 0; ii < nHits; ii++) 
    {
      SQHit* hit = sqHitVector->at(ii);
      HitData dd;

      dd.hit_id      = hit->get_hit_id();
      dd.track_id    = hit->get_track_id();
      dd.detector_id = hit->get_detector_id();
      dd.element_id  = hit->get_element_id();
      dd.tdc_time    = hit->get_tdc_time();
      dd.in_time     = hit->is_in_time();

      lHit.push_back(dd);
    }

  /// Tracklet info (+ Sorted Hit info)
  //cout<<"Tracklet"<<endl;
  lTracklet.clear();
  nTracklets = trackletVector->size();
  //cout<<"nTracklets = "<<nTracklets<<endl;

  for (unsigned int ii = 0; ii < nTracklets; ii++)
    {
      Tracklet* tracklet = trackletVector->at(ii);
      TrackletData dd;
      //cout<<"get data ii = "<<ii<<endl;
      
      dd.n_hits = tracklet->getNHits();
      dd.chisq = tracklet->getChisq();
      dd.charge = tracklet->getCharge();

      for(int jj=0; jj<24; jj++)
	{
	  const int d_id = dID[jj];
	  dd.detector_id[jj] = d_id;
	 	  
	  const float z = gs->getPlanePosition(d_id);
	  dd.detector_zpos[jj] = z;

	  const float x = tracklet->getExpPositionX(z);
	  const float y = tracklet->getExpPositionY(z);

	  dd.x_exp[jj] = x;
          dd.y_exp[jj] = y;
	  
	  if( !gs->isInPlane(d_id, x, y) ) continue;
	    
	  dd.in_plane[jj] = true;
	  
	  const int e_id = gs->getExpElementID(d_id, tracklet->getExpPositionW(d_id));
	  dd.element_id_exp[jj] = e_id;

	  if(e_id < 1 || e_id > gs->getPlaneNElements(d_id) ) continue;
	  

	  SQHit* hit = findHit(d_id, e_id);
	  if( hit == nullptr ) continue;

	  dd.hit_id[jj] = hit->get_hit_id();
	  dd.track_id[jj] = hit->get_track_id();
	  dd.element_id_pos[jj] = hit->get_pos();
	  dd.element_id_closest[jj] = hit->get_element_id();
	}

      lTracklet.push_back(dd);
    } 

  //cout<<"Fill"<<endl;
  tree->Fill();
  //cout<<"process_event done"<<endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int SimpleTree::End(PHCompositeNode* topNode)
{
  cout<<"SimpleTree::End()"<<endl;
  file->cd();
  file->Write();
  file->Close();

  cout<<"SimpleTree::End() - Done"<<endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int SimpleTree::GetNodes(PHCompositeNode *topNode)
{
  sqMCEvent      = findNode::getClass<SQMCEvent     >(topNode, "SQMCEvent");
  sqEvent        = findNode::getClass<SQEvent       >(topNode, "SQEvent");
  sqTrackVector  = findNode::getClass<SQTrackVector >(topNode, "SQTruthTrackVector");
  sqDimuonVector = findNode::getClass<SQDimuonVector>(topNode, "SQTruthDimuonVector");
  sqHitVector    = findNode::getClass<SQHitVector   >(topNode, "SQHitVector");
  trackletVector = findNode::getClass<TrackletVector>(topNode, "TrackletVector");

  if (!sqMCEvent || !sqEvent || !sqTrackVector || !sqDimuonVector || !sqHitVector || !trackletVector) {
  return Fun4AllReturnCodes::ABORTEVENT;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void SimpleTree::SetOutput(const char* out_name)
{
  cout<<"SimpleTree::SetOutput()"<<endl;
  file= new TFile(out_name, "RECREATE");
  tree = new TTree("T", "Created by SimpleTree");
  cout<<out_name<<" -- Created"<<endl;
}

void SimpleTree::MakeTree()
{
  tree->Branch("mcevent", &dMCEvent);
  tree->Branch("event", &dEvent);
  tree->Branch("dimuon", &lDimuon);
  tree->Branch("track", &lTrack);
  tree->Branch("hit", &lHit);
  tree->Branch("tracklet", &lTracklet);
  tree->Branch("nTracks", &nTracks);
  tree->Branch("nHits", &nHits);
  tree->Branch("nTracklets", &nTracklets);
  //tree->Branch("", &);
}

SQHit* SimpleTree::findHit(const int detID, const int eleID)
{
  int delta = 999;
  SQHit* hit = nullptr;
  for(SQHitVector::Iter it = sqHitVector->begin(); it != sqHitVector->end(); ++it)
    {
      if((*it)->get_detector_id() != detID) continue;
      int delta_temp = abs(eleID - (*it)->get_element_id());
      if(delta > delta_temp)
	{
	  delta = delta_temp;
	  hit = (*it);
	}
    }
  return hit;
}
