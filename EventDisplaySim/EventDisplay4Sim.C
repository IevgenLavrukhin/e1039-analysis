#include "G4_InsensitiveVolumes.C"
#include "G4_SensitiveDetectors.C"
#include "G4_Beamline.C"
#include "G4_Target.C"
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
R__LOAD_LIBRARY(libSQPrimaryGen)
R__LOAD_LIBRARY(libpheve_display)
R__LOAD_LIBRARY(libpheve_modules)
R__LOAD_LIBRARY(libpheve_interface)
using namespace std;

Fun4AllServer* se;
EvtFilter* evt_filter;
TGLViewer* view;
TGNumberEntry* ne_evt_id;
TGNumberEntry* ne_trig;

class EvNavHandler
{
  public:
    void NextEvent()
    {
      printf("se->run(1)\n");
      se->run(1, true);
    }
    void ReqEvtID() {
      printf("ReqEvtID: %ld\n",ne_evt_id->GetNumberEntry()->GetIntNumber());
      evt_filter->set_event_id_req((int)(ne_evt_id->GetNumberEntry()->GetIntNumber()));
    }
    void ReqTrig() {
      printf("ReqTrig: %ld\n",ne_trig->GetNumberEntry()->GetIntNumber());
      evt_filter->set_trigger_req((int)(ne_trig->GetNumberEntry()->GetIntNumber()));
    }
    void TopView()
    {
      printf("Top View\n");
      view->ResetCurrentCamera();
      view->CurrentCamera().RotateRad(-3.14/2.0, 0);
      view->CurrentCamera().Zoom(400, 0, 0);
      view->CurrentCamera().Truck(2800,0);
      view->DoDraw();
    }
    void SideView()
    {
      printf("Side View\n");
      view->ResetCurrentCamera();
      view->CurrentCamera().Zoom(400, 0, 0);
      view->CurrentCamera().Truck(2800,0);
      view->DoDraw();
    }
    void View3D()
    {
      printf("3D View\n");
      view->ResetCurrentCamera();
      view->CurrentCamera().RotateRad(-3.14/4., -3.14/4.);
      view->CurrentCamera().Zoom(350, 0, 0);
      view->CurrentCamera().Truck(2000,-1500);
      view->DoDraw();
    }

};

////////////////////////////////////////////////////////////////
//
// Main Function
//
int EventDisplay4Sim(const int nevent = 10)
{
  const double target_coil_pos_z = -300;
  const int nmu = 1;

  const bool do_collimator = true;
  const bool do_target = true;
  const bool do_e1039_shielding = true;
  const bool do_fmag = true;
  const bool do_kmag = true;
  const bool do_absorber = true;

  const double target_l = 7.9; //cm
  const double target_z = (7.9-target_l)/2.; //cm
  const int use_g4steps = 1;

  const double FMAGSTR = -1.054;
  const double KMAGSTR = -0.951;

  const bool gen_pythia8  = true;
  const bool gen_gun      = false;
  const bool gen_particle = false;
  const bool read_hepmc   = false;
  const bool gen_e906legacy = false;

  recoConsts *rc = recoConsts::instance();
  rc->set_DoubleFlag("FMAGSTR", FMAGSTR);
  rc->set_DoubleFlag("KMAGSTR", KMAGSTR);
  rc->Print();

  JobOptsSvc *jobopt_svc = JobOptsSvc::instance();
  jobopt_svc->init("run7_sim.opts");

  GeomSvc::UseDbSvc(true);
  GeomSvc *geom_svc = GeomSvc::instance();

  ///////////////////////////////////////////
  // Make the Server
  //////////////////////////////////////////
  se = Fun4AllServer::instance();
  se->Verbosity(0);

  // pythia8
  if(gen_pythia8) {
    PHPythia8 *pythia8 = new PHPythia8();
    //pythia8->Verbosity(99);
    pythia8->set_config_file("phpythia8_DY.cfg");
    pythia8->set_vertex_distribution_mean(0, 0, target_coil_pos_z, 0);
    pythia8->set_embedding_id(1);
    se->registerSubsystem(pythia8);

    pythia8->set_trigger_AND();

    PHPy8ParticleTrigger* trigger_mup = new PHPy8ParticleTrigger();
    trigger_mup->AddParticles("-13");
    trigger_mup->SetPzHighLow(120, 10);
    pythia8->register_trigger(trigger_mup);

    PHPy8ParticleTrigger* trigger_mum = new PHPy8ParticleTrigger();
    trigger_mum->AddParticles("13");
    trigger_mum->SetPzHighLow(120, 10);
    pythia8->register_trigger(trigger_mum);
  }
  
  if(gen_pythia8 || read_hepmc) {
    HepMCNodeReader *hr = new HepMCNodeReader();
    hr->set_particle_filter_on(true);
    hr->insert_particle_filter_pid(13);
    hr->insert_particle_filter_pid(-13);
    se->registerSubsystem(hr);
  }

  // single gun
  if(gen_gun) {
    PHG4ParticleGun *gun = new PHG4ParticleGun("GUN");
    gun->set_name("mu-");
    gun->set_vtx(30, 10, 590);
    gun->set_mom(-0.3, 2, 50);
    se->registerSubsystem(gun);
  }

  // multi particle gun
  if(gen_particle) {
    PHG4SimpleEventGenerator *genp = new PHG4SimpleEventGenerator("MUP");
    //genp->set_seed(123);
    genp->add_particles("mu+", nmu);  // mu+,e+,proton,pi+,Upsilon
    genp->set_vertex_distribution_function(PHG4SimpleEventGenerator::Uniform,
        PHG4SimpleEventGenerator::Uniform,
        PHG4SimpleEventGenerator::Uniform);
    genp->set_vertex_distribution_mean(0.0, 0.0, target_coil_pos_z);
    genp->set_vertex_distribution_width(0.0, 0.0, 0.0);
    genp->set_vertex_size_function(PHG4SimpleEventGenerator::Uniform);
    genp->set_vertex_size_parameters(0.0, 0.0);

    if(FMAGSTR>0)
      //genp->set_pxpypz_range(0,6, -6,6, 10,100);
      genp->set_pxpypz_range(-3,6, -3,3, 10,100);
    else
      //genp->set_pxpypz_range(-6,0, -6,6, 10,100);
      genp->set_pxpypz_range(-6,3, -3,3, 10,100);


    genp->Verbosity(0);
    se->registerSubsystem(genp);
  }

  if(gen_particle) {
    PHG4SimpleEventGenerator *genm = new PHG4SimpleEventGenerator("MUP");
    //genm->set_seed(123);
    genm->add_particles("mu-", nmu);  // mu+,e+,proton,pi+,Upsilon
    genm->set_vertex_distribution_function(PHG4SimpleEventGenerator::Uniform,
        PHG4SimpleEventGenerator::Uniform,
        PHG4SimpleEventGenerator::Uniform);
    genm->set_vertex_distribution_mean(0.0, 0.0, target_coil_pos_z);
    genm->set_vertex_distribution_width(0.0, 0.0, 0.0);
    genm->set_vertex_size_function(PHG4SimpleEventGenerator::Uniform);
    genm->set_vertex_size_parameters(0.0, 0.0);

    if(FMAGSTR>0)
      //genm->set_pxpypz_range(-6,0, -6,6, 10,100);
      genm->set_pxpypz_range(-6,3, -3,3, 10,100);
    else
      //genm->set_pxpypz_range(0,6, -6,6, 10,100);
      genm->set_pxpypz_range(-3,6, -3,3, 10,100);

    genm->Verbosity(0);
    se->registerSubsystem(genm);
  }

  // E906LegacyGen
  //@
  if(gen_e906legacy){
    SQPrimaryParticleGen *e906legacy = new  SQPrimaryParticleGen();
    const bool pythia_gen = false;
    const bool drellyan_gen = true;
    const bool JPsi_gen = false;
    const bool Psip_gen = false;  
    if(drellyan_gen){
      e906legacy->set_xfRange(0.1, 0.5); //[-1.,1.]
      e906legacy->set_massRange(0.23, 10.0);// 0.22 and above     
      e906legacy->enableDrellYanGen();
    }
    if(Psip_gen){ 
      e906legacy->set_xfRange(0.1, 0.5); //[-1.,1.]
      e906legacy->enablePsipGen();
    }
    if(JPsi_gen){
      e906legacy->set_xfRange(0.1, 0.5); //[-1.,1.]
      e906legacy->enableJPsiGen();
    }
    if(pythia_gen) e906legacy->enablePythia();
    se->registerSubsystem(e906legacy);
  }
  //@

  // Fun4All G4 module
  PHG4Reco *g4Reco = new PHG4Reco();
  //PHG4Reco::G4Seed(123);
  //g4Reco->set_field(5.);
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
  // shape of our world - it is a tube
  g4Reco->SetWorldShape("G4BOX");
  // this is what our world is filled with
  g4Reco->SetWorldMaterial("G4_AIR"); //G4_Galactic, G4_AIR
  // Geant4 Physics list to use
  g4Reco->SetPhysicsList("FTFP_BERT");

  // insensitive elements of the spectrometer
  SetupInsensitiveVolumes(g4Reco, do_e1039_shielding, do_fmag, do_kmag, do_absorber);

  // collimator, targer and shielding between target and FMag
  SetupBeamline(g4Reco, do_collimator, target_coil_pos_z - 302.36); // Is the position correct??

  if (do_target) {
    SetupTarget(g4Reco, target_coil_pos_z, target_l, target_z, use_g4steps);
  }

  // sensitive elements of the spectrometer
  SetupSensitiveDetectors(g4Reco, 0);

  se->registerSubsystem(g4Reco);

  // save truth info to the Node Tree
  PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
  g4Reco->registerSubsystem(truth);

  // Make SQ nodes for truth info
  se->registerSubsystem(new TruthNodeMaker());

  // digitizer
  DPDigitizer *digitizer = new DPDigitizer("DPDigitizer", 0);
  //digitizer->Verbosity(99);
  se->registerSubsystem(digitizer);

  // Trigger Emulator
  DPTriggerAnalyzer* dptrigger = new DPTriggerAnalyzer();
  dptrigger->set_road_set_file_name("$E1039_RESOURCE/trigger/trigger_67.txt");
  //dptrigger->Verbosity(99);
  se->registerSubsystem(dptrigger);

  // Event Filter
  evt_filter = new EvtFilter();
  //evt_filter->Verbosity(10);
  //evt_filter->set_trigger_req(1<<5);
  se->registerSubsystem(evt_filter);

  // trakcing module
  KalmanFastTrackingWrapper *ktracker = new KalmanFastTrackingWrapper();
  //ktracker->Verbosity(99);
  ktracker->set_enable_event_reducer(true);
  ktracker->set_DS_level(0);
  ktracker->set_pattern_db_name(gSystem->ExpandPathName("$E1039_RESOURCE/dsearch/v1/pattern.root"));
  //ktracker->set_sim_db_name(gSystem->ExpandPathName("$E1039_RESOURCE/dsearch/v1/sim.root"));
  //PatternDBUtil::ResScaleDC3(3);
  //PatternDBUtil::LooseMode(false);
  se->registerSubsystem(ktracker);

  VertexFit* vertexing = new VertexFit();
  se->registerSubsystem(vertexing);

  // Event display
  // (width, height, use_fieldmap, use_geofile, field-map name, geo-file name)
  PHEventDisplay* disp = new PHEventDisplay(1920, 1080, false, false, "", "geom.root");
  disp->set_verbosity(3);
  se->registerSubsystem(disp);

  // input - we need a dummy to drive the event loop
  if(read_hepmc) {
    Fun4AllHepMCInputManager *in = new Fun4AllHepMCInputManager("HEPMCIN");
    in->Verbosity(10);
    in->set_vertex_distribution_mean(0,0,target_coil_pos_z,0);
    se->registerInputManager(in);
    in->fileopen("hepmcout.txt");
  } else {
    Fun4AllInputManager *in = new Fun4AllDummyInputManager("DUMMY");
    se->registerInputManager(in);
  }

  //
  // Make GUI
  //
  g4Reco->InitRun(se->topNode());
  disp  ->InitRun(se->topNode());
  view = gEve->GetDefaultGLViewer();

  TEveBrowser* browser = gEve->GetBrowser();
  browser->StartEmbedding(TRootBrowser::kLeft);

  TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
  frmMain->SetWindowName("Event Display");
  frmMain->SetCleanup(kDeepCleanup);

  TGVerticalFrame* frmVert = new TGVerticalFrame(frmMain);
  {

    TGTextButton* b = 0;
    EvNavHandler* handler = new EvNavHandler;
    TGHorizontalFrame* frm1 = 0;

    TGLabel* lab = 0;

    frm1 = new TGHorizontalFrame(frmVert);
    lab = new TGLabel(frm1, "Event ID");
    frm1->AddFrame(lab, new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 5, 5, 3, 4));
    ne_evt_id = new TGNumberEntry(frm1, -1, 9, 999, TGNumberFormat::kNESInteger,
        TGNumberFormat::kNEAAnyNumber,
        TGNumberFormat::kNELLimitMinMax,
        -999999, 999999);
    ne_evt_id->Connect("ValueSet(Long_t)", "EvNavHandler", handler, "ReqEvtID()");
    frm1->AddFrame(ne_evt_id, new TGLayoutHints(kLHintsCenterY | kLHintsRight, 5, 5, 5, 5));
    frmVert->AddFrame(frm1);


    frm1 = new TGHorizontalFrame(frmVert);
    lab = new TGLabel(frm1, "Trigger");
    frm1->AddFrame(lab, new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 5, 5, 3, 4));
    ne_trig = new TGNumberEntry(frm1, -1, 9, 999, TGNumberFormat::kNESInteger,
        TGNumberFormat::kNEAAnyNumber,
        TGNumberFormat::kNELLimitMinMax,
        -999, 999);
    ne_trig->Connect("ValueSet(Long_t)", "EvNavHandler", handler, "ReqTrig()");
    frm1->AddFrame(ne_trig, new TGLayoutHints(kLHintsCenterY | kLHintsRight, 5, 5, 5, 5));
    frmVert->AddFrame(frm1);


    b = new TGTextButton(frmVert, "Next Event");
    frmVert->AddFrame(b, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
    b->Connect("Clicked()", "EvNavHandler", handler, "NextEvent()");

    b = new TGTextButton(frmVert, "Top View");
    frmVert->AddFrame(b, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
    b->Connect("Clicked()", "EvNavHandler", handler, "TopView()");

    b = new TGTextButton(frmVert, "Side View");
    frmVert->AddFrame(b, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
    b->Connect("Clicked()", "EvNavHandler", handler, "SideView()");

    b = new TGTextButton(frmVert, "3D View");
    frmVert->AddFrame(b, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
    b->Connect("Clicked()", "EvNavHandler", handler, "View3D()");

  }
  frmMain->AddFrame(frmVert);

  frmMain->MapSubwindows();
  frmMain->Resize();
  frmMain->MapWindow();

  browser->StopEmbedding();
  browser->SetTabTitle("Event Control", 0);

  return 0;
}
