#include <iomanip>
#include <TFile.h>
#include <TTree.h>

#include <interface_main/SQHitVector.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

#include "AnaModule.h"

int AnaModule::Init(PHCompositeNode* topNode)
{
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
  nHits = 0;
  for(size_t i = 0; i < hitVector->size(); ++i)
  {
    SQHit* hit = hitVector->at(i);

    int elementID = hit->get_element_id();
    edep[nHits] = hit->get_edep();
    vID[nHits] = elementID % (3*6*2);
    hID[nHits] = (elementID - hID[nHits])/(3*6*2);

    ++nHits;
  } 
  saveTree->Fill();

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

int AnaModule::GetNodes(PHCompositeNode *topNode)
{
  hitVector = findNode::getClass<SQHitVector>(topNode, "SQHitVector");
  if(!hitVector)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void AnaModule::MakeTree()
{
  saveFile = new TFile(saveName, "RECREATE");
  saveTree = new TTree("save", "Created by AnaModule");

  saveTree->Branch("eventID", &eventID, "eventID/I");
  saveTree->Branch("nHits", &nHits, "nHits/I");
  saveTree->Branch("edep", edep, "edep[nHits]/D");
  saveTree->Branch("hID", hID, "hID[nHits]/I");
  saveTree->Branch("vID", vID, "vID[nHits]/I");
}