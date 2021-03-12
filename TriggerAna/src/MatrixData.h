#ifndef _MJ_MATRIX_DATA_H__ 
#define _MJ_MATRIX_DATA_H__
#include <TObject.h>

class Road : public TObject
{
 public:
  
  int id;
  int charge; //+1 for +, -1 for -
  int tb; // +1 for top, -1 for bottom
  float pX; 
  float pY;
  float pZ;
  float chisq;
  /// followings are [x,y,dp][st]
  int detectorIDs[3][4]; 
  int elementIDs[3][4];
  float el_pos[3][4];
  float zpos[3][4];
  /// at X hodo stations only
  float exp_xpos[4];
  float exp_ypos[4];

  Road();
  virtual ~Road() {;}

  void clear()
  {
    id = 0; charge =0; pX=0;
  }

  void get_id();

  ClassDef(Road, 1);
};

class RoadPair : public TObject
{
 public:

  Road roads[2]; //[+,-]
  float mass;
  float x1;
  float x2;
  float xf;

  RoadPair();
  virtual ~RoadPair() {;}
  
  void clear()
  { 
    roads[0].clear(); roads[1].clear();
    mass=-9999, x1=-9999, x2=-9999, xf=-9999;
  }

  ClassDef(RoadPair, 1);
};

#endif /* _MJ_MATRIX_DATA_H__ */
