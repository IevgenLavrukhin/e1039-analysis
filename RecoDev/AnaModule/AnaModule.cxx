#include <iomanip>
#include <TFile.h>
#include <TTree.h>

#include <interface_main/SQHitVector.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

#include "AnaModule.h"

int AnaModule::Init(PHCompositeNode* topNode)
{
  pos1 = new TVector3();
  pos2 = new TVector3();
  pos3 = new TVector3();
  mom1 = new TVector3();
  mom2 = new TVector3();
  mom3 = new TVector3();
  mom  = new TVector3();
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaModule::InitRun(PHCompositeNode* topNode)
{
  int ret = GetNodes(topNode);
  if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  eventID = 0;
  MakeTree();
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaModule::process_event(PHCompositeNode* topNode)
{
  int st1 = -1;
  int st2 = -1;
  int st3 = -1;
  for(size_t i = 0; i < hitVector->size(); ++i)
  {
    SQHit* hit = hitVector->at(i);

    int detectorID = hit->get_detector_id();
    if(detectorID == 3)
    {
      pos1->SetXYZ(hit->get_truth_x(), hit->get_truth_y(), hit->get_truth_z());
      mom1->SetXYZ(hit->get_truth_px(), hit->get_truth_py(), hit->get_truth_pz());
      st1 = 1;
    }
    else if(detectorID == 15)
    {
      pos2->SetXYZ(hit->get_truth_x(), hit->get_truth_y(), hit->get_truth_z());
      mom2->SetXYZ(hit->get_truth_px(), hit->get_truth_py(), hit->get_truth_pz());
      st2 = 1;
    }
    else if(detectorID == 21 || detectorID == 27)
    {
      pos3->SetXYZ(hit->get_truth_x(), hit->get_truth_y(), hit->get_truth_z());
      mom3->SetXYZ(hit->get_truth_px(), hit->get_truth_py(), hit->get_truth_pz());
      st3 = 1;
    }
    if(st1 > 0 && st2 > 0 && st3 > 0) break;
  }

  int rec = -1;
  if(recEvent->getNTracks() == 1)
  {
    SRecTrack track = recEvent->getTrack(0);
    *mom = track.getMomentumVecSt1();
    rec = 1;
  }

  if(st1 > 0 && st2 > 0 && st3 > 0 && rec > 0) saveTree->Fill();

  ++eventID;
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaModule::End(PHCompositeNode* topNode)
{
  saveFile->cd();
  saveFile->Write();
  saveFile->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}

void AnaModule::set_output_filename(const TString& n)
{
  saveName = n;
}

int AnaModule::GetNodes(PHCompositeNode* topNode)
{
  hitVector = findNode::getClass<SQHitVector>(topNode, "SQHitVector");
  if(!hitVector)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  recEvent = findNode::getClass<SRecEvent>(topNode, "SRecEvent");
  if(!recEvent)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void AnaModule::MakeTree()
{
  saveFile = new TFile(saveName, "RECREATE");
  saveTree = new TTree("save", "Created by AnaModule");

  saveTree->Branch("pos1", &pos1, 256000, 99);
  saveTree->Branch("pos2", &pos2, 256000, 99);
  saveTree->Branch("pos3", &pos3, 256000, 99);
  saveTree->Branch("mom1", &mom1, 256000, 99);
  saveTree->Branch("mom2", &mom2, 256000, 99);
  saveTree->Branch("mom3", &mom3, 256000, 99);
  saveTree->Branch("mom", &mom, 256000, 99);
}