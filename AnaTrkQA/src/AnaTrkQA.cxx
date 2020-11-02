/**
 * \class AnaTrkQA
 * \ module for track quality assurance
 * \author Abinash Pun
 *
 * 
 */


#include "AnaTrkQA.h"

#include <interface_main/SQHit.h>
#include <interface_main/SQHit_v1.h>
#include <interface_main/SQMCHit_v1.h>
#include <interface_main/SQHitMap_v1.h>
#include <interface_main/SQHitVector_v1.h>
#include <interface_main/SQEvent_v1.h>
#include <interface_main/SQRun_v1.h>
#include <interface_main/SQSpill_v1.h>
#include <interface_main/SQSpillMap_v1.h>

#include <ktracker/SRecEvent.h>
#include <geom_svc/GeomSvc.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4HitDefs.h>
#include <g4main/PHG4VtxPoint.h>

#include <TFile.h>
#include <TTree.h>

#include <cstring>
#include <cmath>
#include <cfloat>
#include <stdexcept>
#include <limits>
#include <tuple>

#include <boost/lexical_cast.hpp>

#define NDET 62
#define LogDebug(exp)		std::cout<<"DEBUG: "  <<__FILE__<<": "<<__LINE__<<": "<< exp << std::endl
#define LogError(exp)		std::cout<<"ERROR: "  <<__FILE__<<": "<<__LINE__<<": "<< exp << std::endl
#define LogWarning(exp)	    std::cout<<"WARNING: "<<__FILE__<<": "<<__LINE__<<": "<< exp << std::endl

using namespace std;

AnaTrkQA::AnaTrkQA(const std::string& name) :
  SubsysReco(name),
  _hit_container_type("Vector"),
  _event(0),
  _run_header(nullptr),
  _spill_map(nullptr),
  _event_header(nullptr),
  _hit_map(nullptr),
  _hit_vector(nullptr),
  _out_name("eval.root")
{
 
}

int AnaTrkQA::Init(PHCompositeNode* topNode) {
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaTrkQA::InitRun(PHCompositeNode* topNode) {

  ResetEvalVars();
  InitEvalTree();

  p_geomSvc = GeomSvc::instance();

  int ret = GetNodes(topNode);
  if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaTrkQA::process_event(PHCompositeNode* topNode) {
  
  int ret = Fun4AllReturnCodes::ABORTRUN;

  if(_recEvent) {    
    //ret = RecoEval(topNode);
    ret = TruthRecoEval(topNode);
    if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;
  }

  ++_event;

  return ret;
}


//play ground for Abi==============================================
int AnaTrkQA::TruthRecoEval(PHCompositeNode* topNode)
{
  ResetEvalVars();
	
  if(_truth) {
    for(auto iter=_truth->GetPrimaryParticleRange().first;
	iter!=_truth->GetPrimaryParticleRange().second;
	++iter) {
      PHG4Particle * par = iter->second; 
     
      int trk_id = par->get_track_id();

      pid[n_tracks] = par->get_pid();
      int vtx_id =  par->get_vtx_id();
		
      PHG4VtxPoint* vtx = _truth->GetVtx(vtx_id);
      gvx[n_tracks] = vtx->get_x();
      gvy[n_tracks] = vtx->get_y();
      gvz[n_tracks] = vtx->get_z();

      TVector3 mom_truth(par->get_px(), par->get_py(), par->get_pz());
      gpx[n_tracks] = par->get_px();
      gpy[n_tracks] = par->get_py();
      gpz[n_tracks] = par->get_pz();
      gpt[n_tracks] = mom_truth.Pt();
      geta[n_tracks] = mom_truth.Eta();
      gphi[n_tracks] = mom_truth.Phi();

       
 /// G4Hits at different stations----  (G4Hit information are not stored in tree for now. Can be extended to use later on)   
      TVector3 g_pos_st1;
      TLorentzVector g_mom_st1;
      bool st1hit = FindG4HitAtStation(trk_id, g4hc_d1x, &g_pos_st1, &g_mom_st1);
      if(st1hit){
	
	gx_st1[n_tracks]  = g_pos_st1.X();
	gy_st1[n_tracks]  = g_pos_st1.Y();
	gz_st1[n_tracks]  = g_pos_st1.Z();
     		
	gpx_st1[n_tracks] = g_mom_st1.Px();
	gpy_st1[n_tracks] = g_mom_st1.Py();
	gpz_st1[n_tracks] = g_mom_st1.Pz();
      }

     TVector3 g_pos_st2;
     TLorentzVector g_mom_st2;
     bool st2hit =  FindG4HitAtStation(trk_id, g4hc_d2xp, &g_pos_st2, &g_mom_st2);
    
     TVector3 g_pos_st3;
     TLorentzVector g_mom_st3;

      bool st3hit = FindG4HitAtStation(trk_id, g4hc_d3px, &g_pos_st3, &g_mom_st3)|| FindG4HitAtStation(trk_id, g4hc_d3mx, &g_pos_st3, &g_mom_st3);
		        
     bool prophit = (FindG4HitAtProp(trk_id,g4hc_p1x1)||FindG4HitAtProp(trk_id,g4hc_p1x2)||FindG4HitAtProp(trk_id,g4hc_p1y1)||FindG4HitAtProp(trk_id,g4hc_p1y2)) && (FindG4HitAtProp(trk_id,g4hc_p2x1)||FindG4HitAtProp(trk_id,g4hc_p2x2)||FindG4HitAtProp(trk_id,g4hc_p2y1)||FindG4HitAtProp(trk_id,g4hc_p2y2));//truth having prop hits at station 4

/// Detector acceptance: Truth particle passing through all the 4 stations------------------
    if(st1hit && st2hit && st3hit && prophit){
    	ac_gpx[n_tracks] = gpx[n_tracks];
    	ac_gpy[n_tracks] = gpy[n_tracks];
    	ac_gpz[n_tracks] = gpz[n_tracks];
    }


///==========Implementing functions to catch best track===
      SRecTrack* Best_recTrack = NULL;
      n_recTracks = _recEvent->getNTracks();
	
     // if(_recEvent->getNTracks()>0) Best_recTrack = FindBestMomRecTrack(_recEvent, mom_truth.Mag());      
  
///@Hit matching condition for choosing reconstructed track ---------------   
       vector<int> sqhit_idvec;
       map<int, vector<int> > rtrkid_hitidvec;

	///fill sqhit hit_id vector    
	 for(int ihit=0; ihit<_hit_vector->size(); ++ihit) {
	 	SQHit *sqhit = _hit_vector->at(ihit);
	   	int sq_detid = sqhit->get_detector_id();
           	if(sq_detid > nChamberPlanes || (sq_detid >= 7 && sq_detid <= 12)) continue;
	   	if(trk_id != sqhit->get_track_id()) continue;
		sqhit_idvec.push_back(sqhit->get_hit_id());
				
	  }
	
         sort(sqhit_idvec.begin(), sqhit_idvec.end());

	///fill map of reco track id and hitindex
	 if(n_recTracks>0){
	 	for(int i = 0; i < n_recTracks; ++i) {
      		 	SRecTrack* recTrack = &_recEvent->getTrack(i);                      
      		 	rtrkid_hitidvec[i] = vector<int>();
      
       			int n_rhits = recTrack->getNHits();
      			for(int j = 0; j < n_rhits; ++j) { 
        		rtrkid_hitidvec[i].push_back(fabs(recTrack->getHitIndex(j)));
			
      			}

      		sort(rtrkid_hitidvec[i].begin(), rtrkid_hitidvec[i].end());
    		}
         }
	
///Now try to find the matching reco id
	double m_matching_threshold = 0.75;
	int rtrkid = -1;
     	unsigned int n_match = 0;
	
	for(auto it = rtrkid_hitidvec.begin(); it != rtrkid_hitidvec.end(); ++it) {
      		int n_match_new = FindCommonHitIDs(sqhit_idvec, it->second);
      		if(n_match < n_match_new) {
        		n_match = n_match_new;
        		rtrkid = it->first;
      		}
    	}

    	if(rtrkid >= 0 && double(n_match)/double(sqhit_idvec.size()) > m_matching_threshold) {
      	Best_recTrack = &_recEvent->getTrack(rtrkid);
    	}
///@Hitmatching...ends...........................................


      if(Best_recTrack){

	TVector3 rec_mom = Best_recTrack->getTargetMom();

	nhits[n_tracks] = Best_recTrack->getNHits();	
	charge[n_tracks] = Best_recTrack->getCharge();

	TVector3 rec_vtx = Best_recTrack->getTargetPos();
 
	rec_vx[n_tracks]  = rec_vtx.X();
	rec_vy[n_tracks]  = rec_vtx.Y();
	rec_vz[n_tracks]  = rec_vtx.Z();
    		
	rec_px[n_tracks]  = rec_mom.Px();
	rec_py[n_tracks]  = rec_mom.Py();
	rec_pz[n_tracks]  = rec_mom.Pz();
	rec_pt[n_tracks]  = rec_mom.Pt();
	rec_eta[n_tracks] = rec_mom.Eta();
	rec_phi[n_tracks] = rec_mom.Phi();

	nhits_st1[n_tracks] = Best_recTrack->getNHitsInStation(1);
	nhits_st2[n_tracks] = Best_recTrack->getNHitsInStation(2);
        nhits_st3[n_tracks] = Best_recTrack->getNHitsInStation(3);
	    
	chisq_st1[n_tracks] = Best_recTrack->getChisq();
	prob_st1[n_tracks] = Best_recTrack->getProb();                
	quality[n_tracks] = Best_recTrack->getQuality(); 


		
	//Old way of getting st. 1 reco values   
	/*double tx, ty, tz;
	Best_recTrack->getMomentumSt1(tx, ty, tz);
	px_st1[n_tracks] = tx;
	py_st1[n_tracks] = ty;
	pz_st1[n_tracks] = tz;
		
	double x, y;
	Best_recTrack->getPositionSt1(x, y);		
        x_st1[n_tracks] = x;		
	y_st1[n_tracks] = y;
        */
      }//if best reco track


      

      ///=================Digitized hit information at different stations and corresponding reco values (if any)

      if(_hit_vector) {
	for(int ihit=0; ihit<_hit_vector->size(); ++ihit) {
	  SQHit *hit = _hit_vector->at(ihit);
           
           int detid = hit->get_detector_id();      
           if(detid > nChamberPlanes || (detid >= 7 && detid <= 12)) continue; //only for the chamber hits
	  
	  /// D0X is considered as st1
	  if(hit->get_track_id() == trk_id && hit->get_detector_id() ==3 ) {
            		
	    sq_px_st1[n_tracks] = hit->get_truth_px();
	    sq_py_st1[n_tracks] = hit->get_truth_py();
	    sq_pz_st1[n_tracks] = hit->get_truth_pz();
		
	    sq_x_st1[n_tracks] = hit->get_truth_x();
	    sq_y_st1[n_tracks] = hit->get_truth_y();
	    sq_z_st1[n_tracks] = hit->get_truth_z();
	    sq_pos_st1[n_tracks]=hit->get_pos();
	    sq_drift_st1[n_tracks] = hit->get_drift_distance();
 
	    ///if the best reco track available
	    if (Best_recTrack){						
	      double sq_z_st1 = hit->get_truth_z();
	      int rec_index_st1 = Best_recTrack->getNearestNode(sq_z_st1);
			
	      double z_rec_st1 = Best_recTrack->getZ(rec_index_st1);
						
	      if(fabs(sq_z_st1- z_rec_st1>1.)) continue;//to avoid mismatch of nodes
			
	      double p_rec_st1 = fabs(1./Best_recTrack->getStateVector(rec_index_st1)[0][0]);
	      double tx_rec_st1 = Best_recTrack->getStateVector(rec_index_st1)[1][0];
	      double ty_rec_st1 = Best_recTrack->getStateVector(rec_index_st1)[2][0];
	      double x_rec_st1 = Best_recTrack->getStateVector(rec_index_st1)[3][0];
	      double y_rec_st1 = Best_recTrack->getStateVector(rec_index_st1)[4][0];

	      double x0_st1 = x_rec_st1 - tx_rec_st1 *z_rec_st1;
	      double y0_st1 = y_rec_st1 - ty_rec_st1 *z_rec_st1;

	      rec_p_st1[n_tracks] =  p_rec_st1;
	      rec_pz_st1[n_tracks] = p_rec_st1/sqrt(1.+tx_rec_st1*tx_rec_st1+ty_rec_st1*ty_rec_st1);
	      rec_px_st1[n_tracks] = rec_pz_st1[n_tracks]* tx_rec_st1;
	      rec_py_st1[n_tracks] = rec_pz_st1[n_tracks]* ty_rec_st1;

	      rec_x_st1[n_tracks] = x_rec_st1;
	      rec_y_st1[n_tracks] = y_rec_st1;
              rec_z_st1[n_tracks] = z_rec_st1;  
             							
	      rec_drift_st1[n_tracks] = p_geomSvc->getDCA(hit->get_detector_id(), hit->get_element_id(),tx_rec_st1, ty_rec_st1, x0_st1,y0_st1);			
	      // Pull distribution work 
	      double cov00_st1 = Best_recTrack->getCovariance(rec_index_st1)[0][0];
	      double sq_mom_st1 = sqrt(sq_px_st1[n_tracks]*sq_px_st1[n_tracks]+sq_py_st1[n_tracks]*sq_py_st1[n_tracks]+sq_pz_st1[n_tracks]*sq_pz_st1[n_tracks]);
	      pull_q2p_st1[n_tracks] = (fabs(Best_recTrack->getStateVector(rec_index_st1)[0][0]) - 1./sq_mom_st1)/sqrt(cov00_st1);		        
		
	  
	    }///best reco condition
	  }///st.1 work ends

	  ///==========

	  ///st. 2 is now D2Xp
	  if(hit->get_track_id() == trk_id && hit->get_detector_id() ==15 ) {
			           		
	    sq_px_st2[n_tracks] = hit->get_truth_px();
	    sq_py_st2[n_tracks] = hit->get_truth_py();
	    sq_pz_st2[n_tracks] = hit->get_truth_pz();
		
	    sq_x_st2[n_tracks] = hit->get_truth_x();
	    sq_y_st2[n_tracks] = hit->get_truth_y();
	    sq_z_st2[n_tracks] = hit->get_truth_z();

	    sq_pos_st2[n_tracks]=hit->get_pos();
	    sq_drift_st2[n_tracks] = hit->get_drift_distance();
 
	    ///if the best reco track available
	    if (Best_recTrack){						
	      double sq_z = hit->get_truth_z();
	      int rec_index = Best_recTrack->getNearestNode(sq_z);
			
	      double rec_z = Best_recTrack->getZ(rec_index);
						
	      if(fabs(sq_z- rec_z>1.)) continue;///to avid mismatch of node
			
	      double p_rec = fabs(1./Best_recTrack->getStateVector(rec_index)[0][0]);
	      double tx_rec = Best_recTrack->getStateVector(rec_index)[1][0];
	      double ty_rec = Best_recTrack->getStateVector(rec_index)[2][0];
	      double x_rec = Best_recTrack->getStateVector(rec_index)[3][0];
	      double y_rec = Best_recTrack->getStateVector(rec_index)[4][0];

	      double x0 = x_rec - tx_rec *rec_z;
	      double y0 = y_rec - ty_rec *rec_z;

	      rec_p_st2[n_tracks] =  p_rec;
	      rec_pz_st2[n_tracks] = p_rec/sqrt(1.+tx_rec*tx_rec+ty_rec*ty_rec);
	      rec_px_st2[n_tracks] = rec_pz_st2[n_tracks]* tx_rec;
	      rec_py_st2 [n_tracks]= rec_pz_st2[n_tracks]* ty_rec;


	      rec_x_st2[n_tracks] = x_rec;
	      rec_y_st2[n_tracks] = y_rec;
	      rec_z_st2[n_tracks] = rec_z;
						
  	      rec_drift_st2[n_tracks] = p_geomSvc->getDCA(hit->get_detector_id(), hit->get_element_id(),tx_rec, ty_rec, x0,y0);			

	      ///Pull distribution work 
	      double cov00_st2 = Best_recTrack->getCovariance(rec_index)[0][0];
	      double sq_mom_st2 = sqrt(sq_px_st2[n_tracks]*sq_px_st2[n_tracks]+sq_py_st2[n_tracks]*sq_py_st2[n_tracks]+sq_pz_st2[n_tracks]*sq_pz_st2[n_tracks]);
	      pull_q2p_st2[n_tracks] = (fabs(Best_recTrack->getStateVector(rec_index)[0][0]) - 1./sq_mom_st2)/sqrt(cov00_st2);		        		
	  
	    }///if best reco track	
	
	  }//st2. work done

	  //=======================
	  			 	

	  //st. 3 is now D3mXp(id = 27) or D3pXp (id=21)
	  if(hit->get_track_id() == trk_id && (hit->get_detector_id() == 27 ||hit->get_detector_id() == 21)) {          		
	  sq_px_st3[n_tracks] = hit->get_truth_px();
	  sq_py_st3[n_tracks] = hit->get_truth_py();
	  sq_pz_st3[n_tracks] = hit->get_truth_pz();
		
	  sq_x_st3[n_tracks] = hit->get_truth_x();
	  sq_y_st3[n_tracks] = hit->get_truth_y();
	  sq_z_st3[n_tracks] = hit->get_truth_z();
	  sq_pos_st3[n_tracks]=hit->get_pos();
	  sq_drift_st3[n_tracks] = hit->get_drift_distance();
 
	  ///if the best reco track available
	  if (Best_recTrack){						
	  double sq_z = hit->get_truth_z();
	  int rec_index = Best_recTrack->getNearestNode(sq_z);
			
	  double rec_z = Best_recTrack->getZ(rec_index);
						
	  if(fabs(sq_z- rec_z>1.)) continue;///to avoid mismatch of nodes
			
	  double p_rec = fabs(1./Best_recTrack->getStateVector(rec_index)[0][0]);
	  double tx_rec = Best_recTrack->getStateVector(rec_index)[1][0];
	  double ty_rec = Best_recTrack->getStateVector(rec_index)[2][0];
	  double x_rec = Best_recTrack->getStateVector(rec_index)[3][0];
	  double y_rec = Best_recTrack->getStateVector(rec_index)[4][0];

	  double x0 = x_rec - tx_rec *rec_z;
	  double y0 = y_rec - ty_rec *rec_z;

	  rec_p_st3[n_tracks] =  p_rec;
	  rec_pz_st3[n_tracks] = p_rec/sqrt(1.+tx_rec*tx_rec+ty_rec*ty_rec);
	  rec_px_st3[n_tracks] = rec_pz_st3[n_tracks]* tx_rec;
	  rec_py_st3[n_tracks] = rec_pz_st3[n_tracks]* ty_rec;

	  rec_x_st3[n_tracks] = x_rec;
	  rec_y_st3[n_tracks] = y_rec;
          rec_z_st3[n_tracks] = rec_z;							
	  rec_drift_st3[n_tracks] = p_geomSvc->getDCA(hit->get_detector_id(), hit->get_element_id(),tx_rec, ty_rec, x0,y0);			

	  // Pull distribution work 
	  double cov00_st3 = Best_recTrack->getCovariance(rec_index)[0][0];
	  double sq_mom_st3 = sqrt(sq_px_st3[n_tracks]*sq_px_st3[n_tracks]+sq_py_st3[n_tracks]*sq_py_st3[n_tracks]+sq_pz_st3[n_tracks]*sq_pz_st3[n_tracks]);
	  pull_q2p_st3[n_tracks] = (fabs(Best_recTrack->getStateVector(rec_index)[0][0]) - 1./sq_mom_st3)/sqrt(cov00_st3);		        		
	  
	  
	  }//if best reco track	
	
	  }//st3. work done
	  
	  //=======================		
			
	}//sqhit vector loop

      }//if hit vector

 	      
      
      ++n_tracks;
      if(n_tracks>=100) break;
             
    }//truth loop
	

  }//truth condition
 
  _qa_tree->Fill();


  return Fun4AllReturnCodes::EVENT_OK;
}

///===========================
int AnaTrkQA::End(PHCompositeNode* topNode) {
  if(Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT)
    std::cout << "AnaTrkQA::End" << std::endl;

  PHTFileServer::get().cd(_out_name.c_str());
  _qa_tree->Write();
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaTrkQA::InitEvalTree() {
 
  PHTFileServer::get().open(_out_name.c_str(), "RECREATE");

  ///For the QA tree======================
  _qa_tree = new TTree("QA_ana", "QA analsysis of reconstruction and simulation");

  _qa_tree->Branch("n_tracks",      &n_tracks,           "n_tracks/I");
  _qa_tree->Branch("n_recTracks",   &n_recTracks,        "n_recTracks/I");

  ///Generated Truth info
  _qa_tree->Branch("pid",           pid,                 "pid[n_tracks]/I");
  _qa_tree->Branch("gvx",           gvx,                 "gvx[n_tracks]/F");
  _qa_tree->Branch("gvy",           gvy,                 "gvy[n_tracks]/F");
  _qa_tree->Branch("gvz",           gvz,                 "gvz[n_tracks]/F");
  _qa_tree->Branch("gpx",           gpx,                 "gpx[n_tracks]/F");
  _qa_tree->Branch("gpy",           gpy,                 "gpy[n_tracks]/F");
  _qa_tree->Branch("gpz",           gpz,                 "gpz[n_tracks]/F");
  _qa_tree->Branch("gpt",           gpt,                 "gpt[n_tracks]/F");
  _qa_tree->Branch("geta",          geta,                "geta[n_tracks]/F");
  _qa_tree->Branch("gphi",          gphi,                "gphi[n_tracks]/F");
/*
  _qa_tree->Branch("gx_st1",        gx_st1,              "gx_st1[n_tracks]/F");
  _qa_tree->Branch("gy_st1",        gy_st1,              "gy_st1[n_tracks]/F");
  _qa_tree->Branch("gz_st1",        gz_st1,              "gz_st1[n_tracks]/F");
  _qa_tree->Branch("gpx_st1",       gpx_st1,             "gpx_st1[n_tracks]/F");
  _qa_tree->Branch("gpy_st1",       gpy_st1,             "gpy_st1[n_tracks]/F");
  _qa_tree->Branch("gpz_st1",       gpz_st1,             "gpz_st1[n_tracks]/F");
*/
  ///Accepted truth info
  _qa_tree->Branch("ac_gpx",        ac_gpx,              "ac_gpx[n_tracks]/F");
  _qa_tree->Branch("ac_gpy",        ac_gpy,              "ac_gpy[n_tracks]/F");
  _qa_tree->Branch("ac_gpz",        ac_gpz,              "ac_gpz[n_tracks]/F");


  ///Reco info in vertex
  _qa_tree->Branch("rec_vx",            rec_vx,                  "rec_vx[n_tracks]/F");
  _qa_tree->Branch("rec_vy",            rec_vy,                  "rec_vy[n_tracks]/F");
  _qa_tree->Branch("rec_vz",            rec_vz,                  "rec_vz[n_tracks]/F");
  _qa_tree->Branch("rec_px",            rec_px,                  "rec_px[n_tracks]/F");
  _qa_tree->Branch("rec_py",            rec_py,                  "rec_py[n_tracks]/F");
  _qa_tree->Branch("rec_pz",            rec_pz,                  "rec_pz[n_tracks]/F");
  _qa_tree->Branch("rec_pt",            rec_pt,                  "rec_pt[n_tracks]/F");
  _qa_tree->Branch("rec_eta",           rec_eta,                 "rec_eta[n_tracks]/F");
  _qa_tree->Branch("rec_phi",           rec_phi,                 "rec_phi[n_tracks]/F");

 
  ///station 1 truth and reco info----- 
  _qa_tree->Branch("sq_pos_st1",      sq_pos_st1,         "sq_pos_st1[n_tracks]/F");
  _qa_tree->Branch("sq_drift_st1",    sq_drift_st1,       "sq_drift_st1[n_tracks]/F");
  _qa_tree->Branch("sq_x_st1",        sq_x_st1,              "sq_x_st1[n_tracks]/F");
  _qa_tree->Branch("sq_y_st1",        sq_y_st1,              "sq_y_st1[n_tracks]/F");
  _qa_tree->Branch("sq_z_st1",        sq_z_st1,              "sq_z_st1[n_tracks]/F");
  _qa_tree->Branch("sq_px_st1",       sq_px_st1,              "sq_px_st1[n_tracks]/F");
  _qa_tree->Branch("sq_py_st1",       sq_py_st1,              "sq_py_st1[n_tracks]/F");
  _qa_tree->Branch("sq_pz_st1",       sq_pz_st1,              "sq_pz_st1[n_tracks]/F");

  _qa_tree->Branch("rec_drift_st1",   rec_drift_st1,      "rec_drift_st1[n_tracks]/F");
  _qa_tree->Branch("rec_px_st1",      rec_px_st1,      "rec_px_st1[n_tracks]/F");
  _qa_tree->Branch("rec_py_st1",      rec_py_st1,      "rec_py_st1[n_tracks]/F");
  _qa_tree->Branch("rec_pz_st1",      rec_pz_st1,      "rec_pz_st1[n_tracks]/F");
  _qa_tree->Branch("rec_x_st1",       rec_x_st1,      "rec_x_st1[n_tracks]/F");
  _qa_tree->Branch("rec_y_st1",       rec_y_st1,      "rec_y_st1[n_tracks]/F");
  _qa_tree->Branch("rec_z_st1",      rec_z_st1,      "rec_z_st1[n_tracks]/F");

 
 /// Station 2 truth and reco info
  _qa_tree->Branch("sq_pos_st2",      sq_pos_st2,         "sq_pos_st2[n_tracks]/F");
  _qa_tree->Branch("sq_drift_st2",    sq_drift_st2,       "sq_drift_st2[n_tracks]/F");
  _qa_tree->Branch("sq_x_st2",        sq_x_st2,              "sq_x_st2[n_tracks]/F");
  _qa_tree->Branch("sq_y_st2",        sq_y_st2,              "sq_y_st2[n_tracks]/F");
  _qa_tree->Branch("sq_z_st2",       sq_z_st2,              "sq_z_st2/[n_tracks]F");
  _qa_tree->Branch("sq_px_st2",        sq_px_st2,              "sq_px_st2[n_tracks]/F");
  _qa_tree->Branch("sq_py_st2",        sq_py_st2,              "sq_py_st2[n_tracks]/F");
  _qa_tree->Branch("sq_pz_st2",        sq_pz_st2,              "sq_pz_st2[n_tracks]/F");

  _qa_tree->Branch("rec_drift_st2",   rec_drift_st2,      "rec_drift_st2[n_tracks]/F");
  _qa_tree->Branch("rec_px_st2",      rec_px_st2,      "rec_px_st2[n_tracks]/F");
  _qa_tree->Branch("rec_py_st2",      rec_py_st2,      "rec_py_st2[n_tracks]/F");
  _qa_tree->Branch("rec_pz_st2",      rec_pz_st2,      "rec_pz_st2[n_tracks]/F");
  _qa_tree->Branch("rec_x_st2",      rec_x_st2,      "rec_x_st2[n_tracks]/F");
  _qa_tree->Branch("rec_y_st2",      rec_y_st2,      "rec_y_st2[n_tracks]/F");
  _qa_tree->Branch("rec_z_st2",      rec_z_st2,      "rec_z_st2[n_tracks]/F");


/// Station 3 truth and reco info
  _qa_tree->Branch("sq_pos_st3",      sq_pos_st3,         "sq_pos_st3[n_tracks]/F");
  _qa_tree->Branch("sq_drift_st3",    sq_drift_st3,       "sq_drift_st3[n_tracks]/F");
  _qa_tree->Branch("sq_x_st3",        sq_x_st3,              "sq_x_st3[n_tracks]/F");
  _qa_tree->Branch("sq_y_st3",        sq_y_st3,              "sq_y_st3[n_tracks]/F");
  _qa_tree->Branch("sq_z_st3",        sq_z_st3,              "sq_z_st3[n_tracks]/F");
  _qa_tree->Branch("sq_px_st3",        sq_px_st3,              "sq_px_st3[n_tracks]/F");
  _qa_tree->Branch("sq_py_st3",        sq_py_st3,              "sq_py_st3[n_tracks]/F");
  _qa_tree->Branch("sq_pz_st3",        sq_pz_st3,              "sq_pz_st3[n_tracks]/F");

  _qa_tree->Branch("rec_drift_st3",   rec_drift_st3,      "rec_drift_st3[n_tracks]/F");
  _qa_tree->Branch("rec_px_st3",      rec_px_st3,      "rec_px_st3[n_tracks]/F");
  _qa_tree->Branch("rec_py_st3",      rec_py_st3,      "rec_py_st3[n_tracks]/F");
  _qa_tree->Branch("rec_pz_st3",      rec_pz_st3,      "rec_pz_st3[n_tracks]/F");
  _qa_tree->Branch("rec_x_st3",      rec_x_st3,      "rec_x_st3[n_tracks]/F");
  _qa_tree->Branch("rec_y_st3",      rec_y_st3,      "rec_y_st3[n_tracks]/F");
  _qa_tree->Branch("rec_z_st3",      rec_z_st3,      "rec_z_st3[n_tracks]/F");


///quality info of reconstructed tracks
  _qa_tree->Branch("pull_q2p_st1",      pull_q2p_st1,      "pull_q2p_st1[n_tracks]/F");
  _qa_tree->Branch("pull_q2p_st2",      pull_q2p_st2,      "pull_q2p_st2[n_tracks]/F");
  _qa_tree->Branch("pull_q2p_st3",      pull_q2p_st3,      "pull_q2p_st3[n_tracks]/F");

  _qa_tree->Branch("chisq_st1",      chisq_st1,            "chisq_st1[n_tracks]/F");
  _qa_tree->Branch("prob_st1",       prob_st1,             "prob_st1[n_tracks]/F");
  _qa_tree->Branch("quality",        quality,             "quality[n_tracks]/F");

  _qa_tree->Branch("nhits",         nhits,               "nhits[n_tracks]/I");
  _qa_tree->Branch("nhits_st1",         nhits_st1,               "nhits_st1[n_tracks]/I");
  _qa_tree->Branch("nhits_st2",         nhits_st2,               "nhits_st2[n_tracks]/I");
  _qa_tree->Branch("nhits_st3",         nhits_st3,               "nhits_st3[n_tracks]/I");
  _qa_tree->Branch("charge",        charge,              "charge[n_tracks]/I");

 

  

  return 0;
}

int AnaTrkQA::ResetEvalVars() {
  run_id = std::numeric_limits<int>::max();
  spill_id = std::numeric_limits<int>::max();
  target_pos = std::numeric_limits<float>::max();
  event_id = std::numeric_limits<int>::max();
  emu_trigger = 0;
  krecstat = std::numeric_limits<int>::max();

  n_hits = 0;
  for(int i=0; i<100; ++i) {
    detector_id[i]    = std::numeric_limits<short>::max();
    element_id[i]     = std::numeric_limits<short>::max();
    hodo_mask[i]      = std::numeric_limits<short>::max();
    drift_distance[i] = std::numeric_limits<float>::max();
    pos[i]            = std::numeric_limits<float>::max();
    detector_z[i]     = std::numeric_limits<float>::max();

    truth_x[i]       = std::numeric_limits<float>::max();
    truth_y[i]       = std::numeric_limits<float>::max();
    truth_z[i]       = std::numeric_limits<float>::max();
    truth_pos[i]     = std::numeric_limits<float>::max();
  }

    n_tracks = 0;
   
  for(int i=0; i<100; ++i) {
    rec_id[i]    = std::numeric_limits<int>::max();
    par_id[i]     = std::numeric_limits<int>::max();
    pid[i]       = std::numeric_limits<int>::max();
    gvx[i]        = std::numeric_limits<float>::max();
    gvy[i]        = std::numeric_limits<float>::max();
    gvz[i]        = std::numeric_limits<float>::max();
    gpx[i]        = std::numeric_limits<float>::max();
    gpy[i]        = std::numeric_limits<float>::max();
    gpz[i]        = std::numeric_limits<float>::max();
    gpt[i]        = std::numeric_limits<float>::max();
    geta[i]       = std::numeric_limits<float>::max();
    gphi[i]       = std::numeric_limits<float>::max();
    gnhits[i]     = std::numeric_limits<int>::max();
    gx_st1[i]     = std::numeric_limits<float>::max();
    gy_st1[i]     = std::numeric_limits<float>::max();
    gz_st1[i]     = std::numeric_limits<float>::max();
    gpx_st1[i]    = std::numeric_limits<float>::max();
    gpy_st1[i]    = std::numeric_limits<float>::max();
    gpz_st1[i]    = std::numeric_limits<float>::max();
    gndc[i]       = std::numeric_limits<int>::max();
    gnhodo[i]     = std::numeric_limits<int>::max();
    gnprop[i]     = std::numeric_limits<int>::max();
    gndp[i]       = std::numeric_limits<int>::max();


    ac_gpx[i]        = std::numeric_limits<float>::max();
    ac_gpy[i]        = std::numeric_limits<float>::max();
    ac_gpz[i]        = std::numeric_limits<float>::max();

    

    /*for(int j=0; j<NDET+1; ++j) {
      gelmid[i][j] = std::numeric_limits<int>::max();
    }
*/
 
    nhits[i]      = std::numeric_limits<int>::max();
    charge[i]     = std::numeric_limits<int>::max();
    rec_vx[i]         = std::numeric_limits<float>::max();
    rec_vy[i]         = std::numeric_limits<float>::max();
    rec_vz[i]         = std::numeric_limits<float>::max();
    rec_px[i]         = std::numeric_limits<float>::max();
    rec_py[i]         = std::numeric_limits<float>::max();
    rec_pz[i]         = std::numeric_limits<float>::max();
    rec_pt[i]         = std::numeric_limits<float>::max();
    rec_eta[i]        = std::numeric_limits<float>::max();
    rec_phi[i]        = std::numeric_limits<float>::max();
  /*  x_st1[i]     = std::numeric_limits<float>::max();
    y_st1[i]     = std::numeric_limits<float>::max();
    px_st1[i]     = std::numeric_limits<float>::max();
    py_st1[i]     = std::numeric_limits<float>::max();
    pz_st1[i]     = std::numeric_limits<float>::max();
*/
    sq_x_st1[i]     = std::numeric_limits<float>::max();
    sq_y_st1[i]     = std::numeric_limits<float>::max();
    sq_z_st1[i]     = std::numeric_limits<float>::max();
    sq_px_st1[i]     = std::numeric_limits<float>::max();
    sq_py_st1[i]     = std::numeric_limits<float>::max();
    sq_pz_st1[i]     = std::numeric_limits<float>::max();
    sq_pos_st1[i]     = std::numeric_limits<float>::max();
    sq_drift_st1[i]     = std::numeric_limits<float>::max(); 
 
    rec_x_st1[i]     = std::numeric_limits<float>::max();
    rec_y_st1[i]     = std::numeric_limits<float>::max();
    rec_z_st1[i]     = std::numeric_limits<float>::max();
    rec_px_st1[i]     = std::numeric_limits<float>::max();
    rec_py_st1[i]     = std::numeric_limits<float>::max();
    rec_pz_st1[i]     = std::numeric_limits<float>::max();
    rec_drift_st1[i]     = std::numeric_limits<float>::max();
 

    sq_x_st2[i]     = std::numeric_limits<float>::max();
    sq_y_st2[i]     = std::numeric_limits<float>::max();
    sq_z_st2[i]     = std::numeric_limits<float>::max();
    sq_px_st2[i]     = std::numeric_limits<float>::max();
    sq_py_st2[i]     = std::numeric_limits<float>::max();
    sq_pz_st2[i]     = std::numeric_limits<float>::max();
    sq_pos_st2[i]     = std::numeric_limits<float>::max();
    sq_drift_st2[i]     = std::numeric_limits<float>::max();

    rec_x_st2[i]     = std::numeric_limits<float>::max();
    rec_y_st2[i]     = std::numeric_limits<float>::max();
    rec_z_st2[i]     = std::numeric_limits<float>::max();
    rec_px_st2[i]     = std::numeric_limits<float>::max();    
    rec_py_st2[i]     = std::numeric_limits<float>::max();
    rec_pz_st2[i]     = std::numeric_limits<float>::max();
    rec_drift_st2[i]     = std::numeric_limits<float>::max();
 

    sq_x_st3[i]     = std::numeric_limits<float>::max();
    sq_y_st3[i]     = std::numeric_limits<float>::max();
    sq_z_st3[i]     = std::numeric_limits<float>::max();
    sq_px_st3[i]     = std::numeric_limits<float>::max();
    sq_py_st3[i]     = std::numeric_limits<float>::max();
    sq_pz_st3[i]     = std::numeric_limits<float>::max();
    sq_pos_st3[i]     = std::numeric_limits<float>::max();
    sq_drift_st3[i]     = std::numeric_limits<float>::max();

    rec_x_st3[i]     = std::numeric_limits<float>::max();
    rec_y_st3[i]     = std::numeric_limits<float>::max();
    rec_z_st3[i]     = std::numeric_limits<float>::max();
    rec_px_st3[i]     = std::numeric_limits<float>::max();
    rec_py_st3[i]     = std::numeric_limits<float>::max();
    rec_pz_st3[i]     = std::numeric_limits<float>::max();
    rec_drift_st3[i]     = std::numeric_limits<float>::max();

    nhits_st1[i] = std::numeric_limits<float>::max();
    nhits_st2[i] = std::numeric_limits<float>::max();
    nhits_st3[i] = std::numeric_limits<float>::max();

  }

  return 0;
}

int AnaTrkQA::GetNodes(PHCompositeNode* topNode) {

  _run_header = findNode::getClass<SQRun>(topNode, "SQRun");
  if (!_run_header) {
    LogError("!_run_header");
    //return Fun4AllReturnCodes::ABORTEVENT;
  }

  _spill_map = findNode::getClass<SQSpillMap>(topNode, "SQSpillMap");
  if (!_spill_map) {
    LogError("!_spill_map");
    //return Fun4AllReturnCodes::ABORTEVENT;
  }

  _event_header = findNode::getClass<SQEvent>(topNode, "SQEvent");
  if (!_event_header) {
    LogError("!_event_header");
    //return Fun4AllReturnCodes::ABORTEVENT;
  }

  if(_hit_container_type.find("Map") != std::string::npos) {
    _hit_map = findNode::getClass<SQHitMap>(topNode, "SQHitMap");
    if (!_hit_map) {
      LogError("!_hit_map");
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  if(_hit_container_type.find("Vector") != std::string::npos) {
    _hit_vector = findNode::getClass<SQHitVector>(topNode, "SQHitVector");
    if (!_hit_vector) {
      LogError("!_hit_vector");
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  _truth = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!_truth) {
    LogError("!_truth");
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _recEvent = findNode::getClass<SRecEvent>(topNode, "SRecEvent");
  if (!_recEvent) {
    LogError("!_recEvent");
    //return Fun4AllReturnCodes::ABORTEVENT;
  }

 
  g4hc_d1x  = findNode::getClass<PHG4HitContainer      >(topNode, "G4HIT_D1X");
  g4hc_d2xp  = findNode::getClass<PHG4HitContainer      >(topNode, "G4HIT_D2Xp");
  g4hc_d3px = findNode::getClass<PHG4HitContainer      >(topNode, "G4HIT_D3pXp");
  g4hc_d3mx = findNode::getClass<PHG4HitContainer      >(topNode, "G4HIT_D3mXp");
  if (! g4hc_d1x) g4hc_d1x = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_D0X");

  if ( !g4hc_d1x || !g4hc_d3px || !g4hc_d3mx) {
    cout << "Failed at getting nodes: "<< g4hc_d1x << " " << g4hc_d3px << " " << g4hc_d3mx << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

//hodoscope for the acceptance study 
  g4hc_h1t  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H1T");
  g4hc_h1b  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H1B");
  g4hc_h2t  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H2T");
  g4hc_h2b  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H2B");
  g4hc_h3t  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H3T");
  g4hc_h3b  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H3B");
  g4hc_h4t  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H4T");
  g4hc_h4b  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H4B");

  if (!g4hc_h1t || !g4hc_h1b || !g4hc_h2t || !g4hc_h2b ||
      !g4hc_h3t || !g4hc_h3b || !g4hc_h4t || !g4hc_h4b   ) {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

///Prop tubes hits
  g4hc_p1y1  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P1Y1");
  g4hc_p1y2  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P1Y2");
  g4hc_p1x1  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P1X1");
  g4hc_p1x2  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P1X2");
  g4hc_p2x1  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P2X1");
  g4hc_p2x2  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P2X2");
  g4hc_p2y1  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P2Y1");
  g4hc_p2y2  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P2Y2");

  if (!g4hc_p1y1 || !g4hc_p1y2 || !g4hc_p1x1 || !g4hc_p1x2 ||
      !g4hc_p2x1 || !g4hc_p2x2 || !g4hc_p2y1 || !g4hc_p2y2   ) {
    return Fun4AllReturnCodes::ABORTEVENT;
  }



  return Fun4AllReturnCodes::EVENT_OK;
}


//For finding g4hit information in stations (following Kenichi's truth node maker)
bool AnaTrkQA::FindG4HitAtStation(const int trk_id, const PHG4HitContainer* g4hc, TVector3* pos, TLorentzVector* mom)
{
  //const double M_MU = 0.1056583745;
  PHG4HitContainer::ConstRange range = g4hc->getHits();
  for (PHG4HitContainer::ConstIterator it = range.first; it != range.second; it++) {
    PHG4Hit* hit = it->second;
    if (hit->get_trkid() == trk_id) {
      pos->SetXYZ (hit->get_x(0)     , hit->get_y(0)     , hit->get_z(0)           );
      mom->SetXYZM(hit->get_px(0),     hit->get_py(0),     hit->get_pz(0), M_MU);
      return true;
    }
  }
  return false;
}


//Function for finding best reco track
SRecTrack* AnaTrkQA::FindBestMomRecTrack(SRecEvent *recEvent,  const float true_TargetP)
{
  double dP = 100.;
  double hold_dP = 99999.;
  
  SRecTrack* Best_recTrack =  NULL;
  for(int itrack=0; itrack<recEvent->getNTracks(); ++itrack){   
    if (hold_dP>dP) hold_dP = dP;
    SRecTrack *recTrack = &recEvent->getTrack(itrack);
    dP = fabs(true_TargetP - recTrack->getTargetMom().Mag());
   
    //Finding out best match track in terms of energy
    if(dP-hold_dP<0.) Best_recTrack = recTrack;  
  }
  return Best_recTrack;
  
}


//Function to find common hit ids for reco and truth tracks
int AnaTrkQA::FindCommonHitIDs(vector<int>& hitidvec1, vector<int>& hitidvec2)
{
  //This function assumes the input vectors have been sorted
  auto iter = hitidvec1.begin();
  auto jter = hitidvec2.begin();

  int nCommon = 0;
  while(iter != hitidvec1.end() && jter != hitidvec2.end()) {
    if(*iter < *jter) {
      ++iter;
    } else {
      if(!(*jter < *iter)) {
        ++nCommon;
        ++iter;
      }
      ++jter;
    }
  }

  return nCommon;
}


//functions for the acceptance
bool AnaTrkQA::FindG4HitAtHodo(const int trk_id, const PHG4HitContainer* g4hc)
{
  //const double M_MU = 0.1056583745; 
  PHG4HitContainer::ConstRange range = g4hc->getHits();
  for (PHG4HitContainer::ConstIterator it = range.first; it != range.second; it++) {
    PHG4Hit* hit = it->second;
    if (hit->get_trkid() == trk_id) {
      return true;
    }
  }
  return false;
}


bool AnaTrkQA::FindG4HitAtProp(const int trk_id, const PHG4HitContainer* g4hc)
{
  PHG4HitContainer::ConstRange range = g4hc->getHits();
  for (PHG4HitContainer::ConstIterator it = range.first; it != range.second; it++) {
    PHG4Hit* hit = it->second;
    if (hit->get_trkid() == trk_id) {
      return true;
    }
  }
  return false;
}

