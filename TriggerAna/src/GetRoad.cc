#include <math.h> 
#include <iomanip>
#include <fstream>
#include <TTree.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <geom_svc/GeomSvc.h>

#include "TreeData.h"
//#include "MatrixData.h" // end up not using at all

#include "GetRoad.h"

using namespace std;

//typedef vector<int, Road> SortRoad;
//typedef vector<int, RoadPair> SortRoadPair;

const vector<int> topIDs = {32, 38, 40, 46};
const vector<int> bottomIDs = {31, 37, 39, 45};

GetRoad::GetRoad(const char* out_name, const int nevt=0) : gs(GeomSvc::instance())
{
  /// Save Road data  
  TFile* fp_out = new TFile(out_name,"RECREATE");
  
  /// Event-by-event road trees
  tr[0] = new TTree("tr_new","new mc event tree for e1039 setup");
  tr[1] = new TTree("tr_tgt","new mc event tree for target at e906 location");
  tr[2] = new TTree("tr_old","old mc event from e906 roadset");
  
  for(int i=0; i<3; i++)
    {
      tr[i] -> Branch("nDimuon",&nDimuon);
      tr[i] -> Branch("mass",&mass);
      tr[i] -> Branch("xf",&xf);
      tr[i] -> Branch("x1",&x1);
      tr[i] -> Branch("x2",&x2);
      tr[i] -> Branch("px",px,"px[2]/F");
      tr[i] -> Branch("track_id",track_id,"track_id[2]/I");  
      tr[i] -> Branch("charge",charge,"charge[2]/I");

      tr[i] -> Branch("nHit",&nHit);
      tr[i] -> Branch("tb",tb,"tb[2]/I");
      tr[i] -> Branch("unique_id",unique_id,"unique_id[2][4]/I");
      tr[i] -> Branch("detector_id",detector_id,"detector_id[2][4]/I");
      tr[i] -> Branch("detector_zpos",detector_zpos,"detector_zpos[2][4]/F");
      tr[i] -> Branch("element_id",element_id,"element_id[2][4]/I");
      tr[i] -> Branch("element_xpos",element_xpos,"element_xpos[2][4]/F");  

      tr[i] -> Branch("road_id",road_id,"road_id[2]/I");
      tr[i] -> Branch("road_freq",road_freq,"road_freq[2]/I");
    }

  /// unique road tree
  tr_road[0] = new TTree("rd_new","new mc road tree for e1039 setup");
  tr_road[1] = new TTree("rd_tgt","new mc road tree for target at e906 location");
  tr_road[2] = new TTree("rd_old","old mc road tree from e906 roadset");

  for(int i=0; i<3; i++)
    {
      tr_road[i] -> Branch("charge",&u_charge);
      tr_road[i] -> Branch("tb",&u_tb);

      tr_road[i] -> Branch("road_id",&u_rid);
      tr_road[i] -> Branch("u_id",u_uid,"u_id[4]/I");
      tr_road[i] -> Branch("e_id",u_eid,"e_id[4]/I");
      tr_road[i] -> Branch("z",u_z,"z[4]/F");
      tr_road[i] -> Branch("x",u_x,"x[4]/F");

      tr_road[i] -> Branch("freq",&u_freq);
      tr_road[i] -> Branch("px_min",&u_px_min);
      tr_road[i] -> Branch("px_max",&u_px_max);
      tr_road[i] -> Branch("px_avg",&u_px_avg);
      tr_road[i] -> Branch("px_rms",&u_px_rms);
    }

  /// MC data
  fp_mc[0] = new TFile("sim_e1039.root"); 
  fp_mc[1] = new TFile("sim_e906.root");

  cout<<"get mc new"<<endl;
  get_mc_roads(0, nevt);

  cout<<"get mc tgt"<<endl;
  get_mc_roads(1, nevt);

  cout<<"get mc old"<<endl;
  get_old_roads(2);
  
  cout<<"Writing objects"<<endl;
  fp_out->cd();
  for(int i=0; i<3; i++)
    tr[i]->Write();

  cout<<"Close file"<<endl;
  fp_out->Close();

  delete fp_out;
  delete fp_mc[0];
  delete fp_mc[1];
  
  cout<<"Done"<<endl;
}

void GetRoad::get_mc_roads(const int idata, const int nevt=0) 
{
  // Read SimpleTree
  DimuonList *lDimuon = new DimuonList();
  HitList * lHit = new HitList();
  
  TTree* T = (TTree*)fp_mc[idata]->Get("T");
  T->SetBranchAddress("dimuon", &lDimuon);
  T->SetBranchAddress("hit", &lHit);
  
  const int nEntries = nevt==0 ? T->GetEntries() : nevt;

  for(int ievt=0; ievt<nEntries; ievt++)
    {
      T->GetEntry(ievt);
      if( ievt % 100000 == 0 ){ cout<<"ievt = "<<ievt<<endl; }

      reset_var();
      
      nDimuon = lDimuon->size();
      nHit = lHit->size();

      if( nDimuon!=1 ){ continue; }

      DimuonData dd = lDimuon->at(0); // either 0 or 1 dimuons
      mass = dd.mass;
      xf = dd.xf;
      x1 = dd.x1;
      x2 = dd.x2;

      charge[0] = 1;
      charge[1] = -1;
 
      px[0] = dd.mom_pos.X();
      px[1] = dd.mom_neg.X();

      track_id[0] = dd.track_id_pos;
      track_id[1] = dd.track_id_neg;

      for(int icharge=0; icharge<2; icharge++)
	{
	  int temp_id[2][4] = {{0,0,0,0}, {0,0,0,0}}; // element_id[top,bottom][station]

	  for(HitList::iterator it = lHit->begin(); it != lHit->end(); it++)
	    {
	      HitData* hh = &(*it);
	      
	      // find a hit from this muon track
	      if( hh->track_id != track_id[icharge] ){ continue; }
	      
	      // find a hodoscope hits
	      if( hh->detector_id == 31 ){ temp_id[1][0] = hh->element_id; }
	      if( hh->detector_id == 37 ){ temp_id[1][1] = hh->element_id; }
	      if( hh->detector_id == 39 ){ temp_id[1][2] = hh->element_id; }
	      if( hh->detector_id == 45 ){ temp_id[1][3] = hh->element_id; }

	      if( hh->detector_id == 32 ){ temp_id[0][0] = hh->element_id; }
	      if( hh->detector_id == 38 ){ temp_id[0][1] = hh->element_id; }
	      if( hh->detector_id == 40 ){ temp_id[0][2] = hh->element_id; }
	      if( hh->detector_id == 46 ){ temp_id[0][3] = hh->element_id; }
	    } // SQHit list
	  
	  // is this track top or bottom?
	  const bool is_top = (temp_id[0][0]*temp_id[0][1]*temp_id[0][2]*temp_id[0][3]);
	  const bool is_bot = (temp_id[1][0]*temp_id[1][1]*temp_id[1][2]*temp_id[1][3]);
	  
	  if( !is_top && !is_bot ){ continue; } // road not in acceptance for this muon track 
	  
	  if( is_top & is_bot )
	    {
	      cerr<<is_top<<" "<<is_bot<<endl;
	      continue;
	    }

	  const int itb = ( is_top ? 0 : 1 );

	  tb[icharge] = 1-2*itb;
	      
	  for(int i=0; i<4; i++)
	    {
	      detector_id[icharge][i] = ( itb==0 ? topIDs[i] : bottomIDs[i] );
	      element_id [icharge][i] = temp_id[itb][i];
	      unique_id  [icharge][i] = 1000*detector_id[icharge][i] + element_id[icharge][i];
	      element_xpos[icharge][i] = gs->getMeasurement( detector_id[icharge][i], element_id[icharge][i] );
	      detector_zpos[icharge][i] = gs->getPlanePosition( detector_id[icharge][i] );
	    }
	  
	  // non zero only if this single muon track is in hodoscope acceptance
	  get_road_id( icharge );
	  road_freq[icharge] = 1; // dummy
	 
	} // each muon track

      tr[idata]->Fill();
      
    }// event loop
  
  delete lDimuon;
  delete lHit;

}


void GetRoad::get_old_roads(const int idata=2)
{

  //cout<<"read previous roads file"<<endl;
  ifstream ifs("/e906/app/users/kimmj/e1039/e1039-analysis/TriggerAna/work/matrix/trigger_67.txt");
  // charge, roadID, uIDs, pXmin, sig, bg

  // save only either pos or neg muon per tree entry
  //cout<<"ifstream loop"<<endl;
  while(!ifs.eof())
    {
      reset_var(); 
      
      char line[255];
      ifs.getline(line,255);
      
      int chg;
      int rID;
      int uID[4];
      float pXmin;
      float freq;

      sscanf(line, "%d %d %d %d %d %d %f %f", &chg, &rID, &uID[0], &uID[1], &uID[2], &uID[3], &pXmin, &freq);
      //cout<<charge<<" "<<roadID<<" "<<px<<" "<<freq<<endl;
      
      const int cc = ( (1-chg)/2 );
      
      road_id[cc] = rID;
      road_freq[cc] = freq;
      px[cc] = pXmin;
      charge[cc] = chg;

      tb[cc] = ( rID>0 ? 1 : -1 );

      for(int i=0; i<4; i++)
	{
	  unique_id[cc][i] = uID[i];
	  detector_id[cc][i] = (uID[i]/1000);
	  element_id[cc][i] = (uID[i]%100);
	  
	  element_xpos[cc][i] = gs->getMeasurement( detector_id[cc][i], element_id[cc][i] );
	  detector_zpos[cc][i] = gs->getPlanePosition( detector_id[cc][i] );
	}
      
      tr[idata]->Fill();      
    }
  //cout<<"close ifstream"<<endl;
}

/*
void GetRoad::get_unique_road()
{
  // for mc
  for(int itr=0; itr<2; itr++)
    {
      const int nevt = tr[itr]->GetEntries();

      for(int ievt=0; ievt<nevt; ievt++)
        {
          tr[itr]->GetEntry(ievt);
	  
	  if( nDimuon!=1 ) continue;
	  
	  if( (road_id[0]*road_id[1])==0 ) continue;
      
	  for(int icharge=0; icharge<2; icharge++)
            {
	      u_charge = 1-2*icharge;
	      u_tb = tb[icharge];
	      u_rid = road_id[icharge];

              for(int i=0; i<4; i++)
		{ 
		  u_uid[i] = unique_id[icharge][i]; 
		  u_eid[i] = element_id[icharge][i];
		  u_z[i] = detector_zpos[icharge][i];
		  u_x[i] = element_xpos[icharge][i];
		}



              it[icharge] = find(comp_id[icharge].begin(), comp_id[icharge].end(), rid)\
		;

              if( (it[icharge] == comp_id.end()) || (comp_id[icharge].size()==0) )
                {
                  comp_id[icharge].push_back(rid);
                }

	    }
*/
int GetRoad::channels_to_roadID(const vector<int> detectorIDs, const vector<int> elementIDs)
/// argument should be passed after sorted by detectorID number from station 1 to 4
{
  int tb=0;

  if( detectorIDs == topIDs )
    {
      tb = 1;
    }
  else if( detectorIDs == bottomIDs )
    {
      tb = -1;
    }
  else
    {
      return 0;
    }

  return tb*( (elementIDs[0]-1)*pow(16,3) + (elementIDs[1]-1)*pow(16,2) + (elementIDs[2]-1)*16 + elementIDs[3] );
}


bool GetRoad::roadID_to_channels(const int roadID, vector<int> &detectorIDs, vector<int> &elementIDs)
{
  if( roadID==0 ) return false;

  detectorIDs.clear();
  elementIDs.clear();

  detectorIDs = (roadID > 0) ? topIDs : bottomIDs;

  const int _roadID = fabs(roadID);

  elementIDs.push_back( _roadID/(int)pow(16,3) + 1 );
  if( elementIDs[0] > 24 ) return false;

  elementIDs.push_back( _roadID%(int)pow(16,3)/(int)pow(16,2) + 1 );
  if( elementIDs[1] > 16 ) return false;

  elementIDs.push_back( _roadID%(int)pow(16,2)/16 + 1 );
  if( elementIDs[2] > 16 ) return false;

  elementIDs.push_back( _roadID%16 );
  if( elementIDs[3] > 16 ) return false;

  return true;
};

void GetRoad::reset_var()
{
  mass = -9999;
  xf = -9999;
  x1 = -9999;
  x2 = -9999;
  nDimuon = 0;
  nHit = 0;

  for(int i=0; i<2; i++)
    {
      road_id[i] = 0;
      road_freq[i] = 0;
      tb[i] = 0;
      px[i] = 0;
      charge[i] = 0;
      track_id[i] = -9999;

      for(int j=0; j<4; j++)
	{
	  unique_id[i][j] = 0;
	  detector_id[i][j] = 0;
	  element_id[i][j] = 0;
	  detector_zpos[i][j] = -9999;
	  element_xpos[i][j] = -9999;
	}
    }
}

void GetRoad::get_road_id(const int i)
{
  road_id[i] = tb[i]*( (element_id[i][0]-1)*pow(16,3) + (element_id[i][1]-1)*pow(16,2) + (element_id[i][2]-1)*16 + element_id[i][3] );
}
