/// Fun4SimDst.C:  Fun4all macro to analyze the E1039 simulated DST files.
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
R__LOAD_LIBRARY(libinterface_main)
R__LOAD_LIBRARY(libfun4all)
R__LOAD_LIBRARY(libg4detectors)
R__LOAD_LIBRARY(libg4eval)
R__LOAD_LIBRARY(libktracker)
R__LOAD_LIBRARY(libtrigger_ana)
#endif

int Fun4SimDst(const int   n_dst_ana=0,
               const char* fn_list_dst="list_dst.txt",
               const char* fn_udst="uDST.root")
{
  //gSystem->Load("libtrigger_ana.so");
 
  GeomSvc::UseDbSvc(true);
  GeomSvc* geom_svc = GeomSvc::instance();

  Fun4AllServer* se = Fun4AllServer::instance();
  //se->Verbosity(1);
  Fun4AllInputManager *in = new Fun4AllDstInputManager("SimDst");
  se->registerInputManager(in);


  Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", fn_udst);
  se->registerOutputManager(out);
  out->AddNode("SQEvent");
  out->AddNode("SQMCEvent");
  //out->AddNode("SRecEvent");
  out->AddNode("SQTruthTrackVector");
  out->AddNode("SQTruthDimuonVector");
  out->AddNode("SQHitVector");
  out->AddNode("TrackletVector");


  vector<string> list_dst;
  string fn_dst;
  ifstream ifs(fn_list_dst);
  while (ifs >> fn_dst) list_dst.push_back(fn_dst);
  ifs.close();


  int n_dst = list_dst.size();
  cout << "N of DSTs: all = " << n_dst;
  if (n_dst_ana > 0 && n_dst > n_dst_ana) n_dst = n_dst_ana;
  cout << ", to be analyzed = " << n_dst << endl;
  for (int i_dst = 0; i_dst < n_dst; i_dst++) {
    string fn_dst = list_dst[i_dst];
    cout << "DST: " << i_dst << "/" << n_dst << ": " << fn_dst << endl;

    in->fileopen(fn_dst);

    se->run();

  }

  
  se->End();
  delete se;
  return 0;
}
