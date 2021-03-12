#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
R__LOAD_LIBRARY(libtrigger_ana)
#endif

int run_GetRoad(const char* fname = "roads.root")
{
  GeomSvc::UseDbSvc(true);
  GeomSvc* geom_svc = GeomSvc::instance();

  GetRoad *gr = new GetRoad(fname, 0);

  delete gr;
  return 0;
}
