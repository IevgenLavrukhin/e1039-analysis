void gen_roadset()
{
  TFile* fp = new TFile("roads.root");
  TTree* tr[3];
  tr[0] = (TTree*)fp->Get("tr_new");
  tr[1] = (TTree*)fp->Get("tr_tgt");
  tr[2] = (TTree*)fp->Get("tr_old");

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

  for(int itr=0; itr<3; itr++)
    {
      char name[255];
      sprintf(name,"roads_%d.txt",itr);

      FILE* fo = fopen(name, "w");

      tr[itr] -> SetBranchAddress("nDimuon",&nDimuon);
      tr[itr] -> SetBranchAddress("mass",&mass);
      tr[itr] -> SetBranchAddress("xf",&xf);
      tr[itr] -> SetBranchAddress("x1",&x1);
      tr[itr] -> SetBranchAddress("x2",&x2);
      tr[itr] -> SetBranchAddress("px",px);
      tr[itr] -> SetBranchAddress("track_id",track_id);
      tr[itr] -> SetBranchAddress("charge",charge);
      tr[itr] -> SetBranchAddress("nHit",&nHit);
      tr[itr] -> SetBranchAddress("tb",tb);
      tr[itr] -> SetBranchAddress("unique_id",unique_id);
      tr[itr] -> SetBranchAddress("detector_id",detector_id);
      tr[itr] -> SetBranchAddress("detector_zpos",detector_zpos);
      tr[itr] -> SetBranchAddress("element_id",element_id);
      tr[itr] -> SetBranchAddress("element_xpos",element_xpos);
      tr[itr] -> SetBranchAddress("road_id",road_id);
      tr[itr] -> SetBranchAddress("road_freq",road_freq);
      
      
      const int nevt = tr[itr]->GetEntries();
      for(int ievt=0; ievt<nevt; ievt++)
	{
	  tr[itr]->GetEntry(ievt);
	  
	  if( itr<2 )
	    {
	      if( nDimuon!=1 ) continue;
	  
	      if( (road_id[0]*road_id[1])==0 ) continue;
	    }
	  
	  for(int icharge=0; icharge<2; icharge++)
	    {
	      int rid = road_id[icharge];
	      int uid[4];
	      for(int i=0; i<4; i++){ uid[i] = unique_id[icharge][i]; }
	      
	      it[icharge] = find(comp_id[icharge].begin(), comp_id[icharge].end(), rid);
	      
	      if( (it[icharge] == comp_id.end()) || (comp_id[icharge].size()==0) )
		{
		  comp_id[icharge].push_back(rid);
		}
	      


	      fprintf(fo, "%d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t \n",charge, rid, uid[0], uid[1], uid[2], uid[3], 0, 0, 0);
	    }
	}
      fclose(fo);
      
      cout<<"Done for loop "<<itr<<endl;
    }
}

	  
	  
