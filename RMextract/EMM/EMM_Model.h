//---------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

extern "C"{
#include "GeomagnetismHeader.h"
}



//---------------------------------------------------------------------------

/* The Geomagnetism Library is used to make a command prompt program. The program prompts
the user to enter a location, performs the computations and prints the results to the
standard output. The program expects the files GeomagnetismLibrary.c, GeomagnetismHeader.h,
EMM2010.COF, EMM2010SV.COF and EGM9615.h to be in the same directory. 

Manoj.C.Nair
Nov 23, 2009

 *  Revision Number: $Revision: 831 $
 *  Last changed by: $Author: awoods $
 *  Last changed on: $Date: 2012-04-13 14:08:09 -0600 (Fri, 13 Apr 2012) $

 */

class WMM_Model{
 public:
  ~WMM_Model();
  WMM_Model();
  WMM_Model(const char CoeffFile[256],const float date=2010.,const float lon=0.,const float lat=0.,const float h=0.);
  void setEM();
  void setDate(const float date);
  void setLonLat(const float lon,const float lat);
  void setHeight(const float h);
  float getX();
  float getY();
  float getZ();
  float getdX(){return GeoMagneticElements.Xdot;};
  float getdY(){return GeoMagneticElements.Ydot;};
  float getdZ(){return GeoMagneticElements.Zdot;};
  float getDate(){return UserDate.DecimalYear;};
  float getLon(){return CoordGeodetic.lambda ;};
  float getLat(){return CoordGeodetic.phi ;};
  float getHeight(){return CoordGeodetic.HeightAboveEllipsoid; };
  int NumTerms, LoadedEpoch, nMax, nMaxEMM;

  MAGtype_Ellipsoid Ellip;
  MAGtype_CoordSpherical CoordSpherical;
  MAGtype_CoordGeodetic CoordGeodetic;
  MAGtype_Date UserDate;
  MAGtype_GeoMagneticElements GeoMagneticElements;
  MAGtype_Geoid Geoid;
 private:
  static const int epochs = 1;
  MAGtype_MagneticModel * MagneticModels[epochs];
  MAGtype_MagneticModel * TimedMagneticModel;
};


class EMM_Model: public WMM_Model{
 public:
  EMM_Model(const char CoeffFile[256], const float date,const float lon,const float lat,const float h);
  void setEM();
 private:
  static const int epochs = 11;
  MAGtype_MagneticModel *MagneticModels[epochs+1];
  MAGtype_MagneticModel *TimedMagneticModel;

 };
