/// Fun4SimMicroDst.C:  Fun4all macro to analyze you uDST.
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
R__LOAD_LIBRARY(libinterface_main)
R__LOAD_LIBRARY(libfun4all)
R__LOAD_LIBRARY(libg4detectors)
R__LOAD_LIBRARY(libg4eval)
R__LOAD_LIBRARY(libktracker)
R__LOAD_LIBRARY(libtrigger_ana)
#endif

int Fun4SimMicroDst(const char* fn_udst="uDST.root", const char* out_name = "sim.root")
{
  //gSystem->Load("libtrigger_ana.so");

  GeomSvc::UseDbSvc(true);
  GeomSvc* geom_svc = GeomSvc::instance();

  Fun4AllServer* se = Fun4AllServer::instance();
  //se->Verbosity(1);
  Fun4AllInputManager *in = new Fun4AllDstInputManager("SimMicroDst");
  in->fileopen(fn_udst);
  se->registerInputManager(in);

  SimpleTree *st = new SimpleTree();
  st->SetOutput(out_name);
  se->registerSubsystem(st);
  
  se->run();
  se->End();

  delete se;
  return 0;
}
