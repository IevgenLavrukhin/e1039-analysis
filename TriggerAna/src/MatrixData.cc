#include <cmath>
#include "MatrixData.h"
using namespace std;

Road::Road()
  : id(0), charge(0), pX(-9999), pY(-9999), pZ(-9999), chisq(-9999)
{
  for(int i=0; i<4; i++)
    {
      exp_xpos[i] = -9999;
      exp_ypos[i] = -9999;
      
      for(int j=0; j<3; j++)
	{
	  detectorIDs[j][i] = 0;
	  elementIDs[j][i] = 0;
	  el_pos[j][i] = -9999;
          zpos[j][i] = -9999;
	}
    }
}

void Road::get_id()
{
  int tb = 0;

  if( detectorIDs[0][0] == 32 && detectorIDs[0][1] == 38 && detectorIDs[0][2] == 40 && detectorIDs[0][3] == 46)
    {
      tb = 1;
    }
  else if( detectorIDs[0][0] == 31 && detectorIDs[0][1] == 37 && detectorIDs[0][2] == 39 && detectorIDs[0][3] == 45)
    {
      tb = -1;
    }

  id = tb*( (elementIDs[0][0]-1)*pow(16,3) + (elementIDs[0][1]-1)*pow(16,2) + (elementIDs[0][2]-1)*16 + elementIDs[0][3] );
}

RoadPair::RoadPair()
  : mass(-9999), x1(-9999), x2(-9999), xf(-9999)
{
  ;
}
