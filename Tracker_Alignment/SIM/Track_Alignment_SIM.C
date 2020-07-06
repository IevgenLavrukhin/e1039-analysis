#include "G4_InsensitiveVolumes.C"
#include "G4_SensitiveDetectors.C"
#include "G4_Beamline.C"
#include "G4_Target.C"
#include <TSystem.h>
#include "Alignment/Alignment.h"

R__LOAD_LIBRARY(libfun4all)
R__LOAD_LIBRARY(libPHPythia8)
R__LOAD_LIBRARY(libg4detectors)
R__LOAD_LIBRARY(libg4testbench)
R__LOAD_LIBRARY(libg4eval)
R__LOAD_LIBRARY(libg4dst)
R__LOAD_LIBRARY(libdptrigger)
R__LOAD_LIBRARY(libembedding)
R__LOAD_LIBRARY(libevt_filter)
R__LOAD_LIBRARY(libktracker)
R__LOAD_LIBRARY(libAlignment)



using namespace std;


int Track_Alignment_SIM(const int nevent = 10)
{
  const bool dimuon = false;
  const bool single = !dimuon;

  const bool do_collimator = true;
  const bool do_target     = true;
  const bool do_shielding  = true;
  const bool do_fmag       = true;
  const bool do_kmag       = true;
  const bool do_absorber   = true;

  const double collimator_pos_z = -602.36;
  const double target_coil_pos_z = -300.;
  const double target_l = 7.9; //cm
  const double target_z = (7.9-target_l)/2.; //cm

  const double FMAGSTR = -1.054;
  const double KMAGSTR = -0.951;

  recoConsts *rc = recoConsts::instance();
  rc->set_DoubleFlag("FMAGSTR", FMAGSTR);
  rc->set_DoubleFlag("KMAGSTR", KMAGSTR);
  rc->set_IntFlag("RANDOMSEED", 12345);
	rc->Print();


  JobOptsSvc* jobopt_svc = JobOptsSvc::instance();
  jobopt_svc->init("config_files/e1039_sim.opts");


  GeomSvc::UseDbSvc(true);
  GeomSvc* geom_svc = GeomSvc::instance();

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

	//use pythia for dimuon generaion. Also should be changes in the configuration file!!!!
	if(dimuon){
		printf("Dimuon is not implemented here!");
		return -1;
	}
  if(single)   //change the hard-coded numbers to change the initial vertex/momentum distribution
  {
    PHG4SimpleEventGenerator* genp = new PHG4SimpleEventGenerator("MUP");
    genp->set_seed(123);
    genp->add_particles("mu+", 1);  // mu+,e+,proton,pi+,Upsilon
    genp->set_vertex_distribution_function(PHG4SimpleEventGenerator::Uniform, PHG4SimpleEventGenerator::Uniform, PHG4SimpleEventGenerator::Uniform);
    genp->set_vertex_distribution_mean(0.0, 0.0, target_coil_pos_z);
    genp->set_vertex_distribution_width(0.0, 0.0, 0.0);
    genp->set_vertex_size_function(PHG4SimpleEventGenerator::Uniform);
    genp->set_vertex_size_parameters(0.0, 0.0);
    genp->set_pxpypz_range(-6., 6., -3. ,3., 25., 100.);
    se->registerSubsystem(genp);
  }

  //------------------------------------- Fun4All G4 module ---------------------------------------------------//
	printf("Configuring G4 Module for Fun4All .......\n");
	
  PHG4Reco* g4Reco = new PHG4Reco();
	g4Reco->set_field_map(
      jobopt_svc->m_fMagFile+" "+
      jobopt_svc->m_kMagFile+" "+
      Form("%f",FMAGSTR) + " " +
      Form("%f",KMAGSTR) + " " +
      "5.0",
      PHFieldConfig::RegionalConst);
 	// size of the world - every detector has to fit in here
  g4Reco->SetWorldSizeX(1000);
  g4Reco->SetWorldSizeY(1000);
	g4Reco->SetWorldSizeZ(5000);
  // shape of our world - it is a box
  g4Reco->SetWorldShape("G4BOX");
	// this is what our world is filled with
  g4Reco->SetWorldMaterial("G4_AIR"); //G4_Galactic, G4_AIR
  // Geant4 Physics list to use
  g4Reco->SetPhysicsList("FTFP_BERT");

  SetupBeamline(g4Reco, do_collimator, collimator_pos_z);
  SetupTarget(g4Reco, target_coil_pos_z, target_l, target_z, 1, 0);
  SetupInsensitiveVolumes(g4Reco, do_shielding, do_fmag, do_kmag, do_absorber);
  SetupSensitiveDetectors(g4Reco);
  se->registerSubsystem(g4Reco);

  printf("\t\t Done!\n");

  //--------------------------------  Save truth info to the Node Tree -----------------------------------------//
	printf("Saving TRUTH info into Node Tree .......\n");

  PHG4TruthSubsystem* truth = new PHG4TruthSubsystem();
  g4Reco->registerSubsystem(truth);

  printf("\t\t Done!\n");


  //-------------------------------------- Apply in-acceptance cut ---------------------------------------------//
	printf("Apply in-acceptance cut .......\n");

  RequireParticlesInAcc* inacc = new RequireParticlesInAcc();
  if(dimuon)
  {
    inacc->SetNumParticlesPerEvent(2);
  }
  else if(single)
  {
    inacc->SetNumParticlesPerEvent(1);
  }
  se->registerSubsystem(inacc);

  printf("\t\t Done!\n");

  //-------------------------------------- Digitizer 	----------------------------------------------------------//
  printf("Register Digitizer .......\n");

	DPDigitizer* digitizer = new DPDigitizer("DPDigitizer", 0);
  se->registerSubsystem(digitizer);

  printf("\t\t Done!\n");


  //-------------------------------------- Tracking module  ---------------------------------------------------//
  printf("Setup Tracking module .......\n");

 	SQReco* reco = new SQReco();
  reco->Verbosity(0);
  reco->set_enable_KF(true);           //Kalman filter not needed for the track finding, disabling KF saves a lot of initialization time
  reco->setInputTy(SQReco::E1039);     //options are SQReco::E906 and SQReco::E1039
  reco->setFitterTy(SQReco::KFREF);    //not relavant for the track finding
  reco->set_evt_reducer_opt("none");   //if not provided, event reducer will be using JobOptsSvc to intialize; to turn off, set it to "none", for normal tracking, set to something like "aoc"
  reco->set_enable_eval(true);          //set to true to generate evaluation file which includes final track candidates 
  reco->set_eval_file_name("eval.root");
  reco->set_enable_eval_dst(false);     //set to true to include final track cnadidates in the DST tree
  se->registerSubsystem(reco);


  printf("\t\t Done!\n");



	//-------------------------------------- Alignment ----------------------------------------------------------//
  printf("Setup Alignment module .......\n");

	Alignment* al = new Alignment();
	se-> registerSubsystem(al);

  printf("\t\t Done!\n");

  //--------------------------------------- OUTPUT: -----------------------------------------------------------//
  Fun4AllDstOutputManager* out = new Fun4AllDstOutputManager("DSTOUT", "DST.root");
  se->registerOutputManager(out);

  se->run(nevent);

  PHGeomUtility::ExportGeomtry(se->topNode(), "geom.root");

	se->End();
  se->PrintTimer();
	printf("ALL DONE! :) \n");

  delete se;
  gSystem->Exit(0);



	return 0; 
}



