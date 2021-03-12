#include "TreeData.h"
using namespace std;

MCEventData::MCEventData()
  : process_id(-9999)
{
  for(int ii=0; ii<4; ii++)
    {
      particle_id[ii] = -9999;
      particle_mom[ii].SetXYZT(0, 0, 0, 0);
    }
}

EventData::EventData()
  : event_id(-9999), 
    trig_bits(-9999)
{
  ;
}

DimuonData::DimuonData()
  : dim_id(-9999),
    pdg_id(-9999),
    mass(-9999),
    x1(-9999),
    x2(-9999),
    xf(-9999),
    track_id_pos(-9999),
    track_id_neg(-9999)
{
  mom.SetXYZT(0, 0, 0, 0);
  mom_pos.SetXYZT(0, 0, 0, 0);
  mom_neg.SetXYZT(0, 0, 0, 0);
}

TrackData::TrackData()
  : track_id(-9999),
    charge(-9999),
    zvtx(-9999)
{
  ;
}

HitData::HitData()
  : hit_id(-9999),
    track_id(-9999),
    detector_id(-9999),
    element_id(-9999),
    tdc_time(-9999),
    in_time(false)
{
  ;
}

TrackletData::TrackletData()
  : n_hits(-9999),
    chisq(-9999),
    charge(-9999)
{
  for(int ii=0; ii<24; ii++)
    {
      detector_id[ii] = -9999;
      detector_zpos[ii] = -9999;
      x_exp[ii] = -9999;
      y_exp[ii] = -9999;
      in_plane[ii] = false;
      element_id_exp[ii] = -9999; 
      hit_id[ii] = -9999;
      track_id[ii] = -9999;
      element_id_pos[ii] = -9999; 
      element_id_closest[ii] = -9999; 
    }
}
