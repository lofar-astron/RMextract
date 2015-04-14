#include "EMM_Model.h"
#include "EGM9615.h"
#include <assert.h>



WMM_Model::~WMM_Model(){
  if (TimedMagneticModel){
    MAG_FreeMagneticModelMemory(TimedMagneticModel);
    TimedMagneticModel=NULL;
  }
  for(int i = 0; i < epochs + 1; i++) 
    if (MagneticModels[i]){
      MAG_FreeMagneticModelMemory(MagneticModels[i]);
      MagneticModels[i]=NULL;
    }
};

WMM_Model::WMM_Model(){  
  TimedMagneticModel=NULL;
  
  for(int i = 0; i < epochs + 1; i++) 
      MagneticModels[i]=NULL;
};

WMM_Model::WMM_Model(const char CoeffFile[256], const float date,const float lon,const float lat,const float h){
  TimedMagneticModel=NULL;
  
  for(int i = 0; i < epochs + 1; i++) 
      MagneticModels[i]=NULL;
    
  char filename[256];
  char filenameSV[256];
  
  sprintf(filename, "%s.COF", CoeffFile);
  
  FILE * filename_valid=fopen (filename, "r");
  assert(filename_valid);
  fclose(filename_valid);
  nMax = 0;
  MAG_robustReadMagModels(filename, &MagneticModels, epochs);
  /*Create Main Field Model for EMM2010*/
  if(nMax < MagneticModels[0]->nMax) nMax = MagneticModels[0]->nMax;
  NumTerms = ((nMax + 1) * (nMax + 2) / 2);
  TimedMagneticModel = MAG_AllocateModelMemory(NumTerms); /* For storing the time modified WMM Model parameters */

  MAG_SetDefaults(&Ellip, &Geoid); /* Set default values and constants */
  /* Set EGM96 Geoid parameters */
  Geoid.GeoidHeightBuffer = GeoidHeightBuffer;
  Geoid.Geoid_Initialized = 1;
  setDate(date);
  setLonLat(lon,lat);
  setHeight(h);
  setEM();
};

void WMM_Model::setEM()
  {
    int index;
    MAG_GeodeticToSpherical(Ellip, CoordGeodetic, &CoordSpherical); /*Convert from geodetic to Spherical Equations: 17-18, WMM Technical report*/
    MAG_TimelyModifyMagneticModel(UserDate, MagneticModels[0], TimedMagneticModel); /* Time adjust the coefficients, Equation 19, WMM Technical report */
    MAG_Geomag(Ellip, CoordSpherical, CoordGeodetic, TimedMagneticModel, &GeoMagneticElements); /* Computes the geoMagnetic field elements and their time change*/
    MAG_CalculateGridVariation(CoordGeodetic, &GeoMagneticElements);
  };

void WMM_Model::setDate(const float date){
    UserDate.DecimalYear=date; //decimal year
    
}

void WMM_Model::setLonLat(const float lon,const float lat){
  CoordGeodetic.phi=lat;//lattitude in degrees
  CoordGeodetic.lambda=lon; //longitude in degrees
    
}
void WMM_Model::setHeight(const float h){
    CoordGeodetic.HeightAboveEllipsoid=h; //Height above WGS-84 elipsoid
    Geoid.UseGeoid=0;
}

float WMM_Model::getX(){
    return GeoMagneticElements.X;
}

float WMM_Model::getY(){
    return GeoMagneticElements.Y;
}

float WMM_Model::getZ(){
    return GeoMagneticElements.Z;
}


EMM_Model::EMM_Model(const char CoeffFile[256], const float date,const float lon,const float lat,const float h){
  TimedMagneticModel=NULL;
  
  for(int i = 0; i < epochs + 1; i++) 
      MagneticModels[i]=NULL;
  char filename[256];
  char filenameSV[256];
  
  sprintf(filename, "%s%d.COF", CoeffFile,2010);
  sprintf(filenameSV, "%s%dSV.COF",CoeffFile,2010);
  //check if filename exists
  
  FILE * filename_valid=fopen (filename, "r");
  assert(filename_valid);
  fclose(filename_valid);
  //initialize the model for quick acces
  for(int Epoch = epochs; Epoch >= 0; Epoch--)
    {
      if(Epoch == (epochs-1))
	Epoch--;
      if(Epoch < 10)
	{
	  sprintf(filename, "%s%d.COF", CoeffFile,Epoch + 2000);
	  sprintf(filenameSV, "%s%dSV.COF",CoeffFile,Epoch + 2000);
	}

      MAG_robustReadMagneticModel_Large(filename, filenameSV, &MagneticModels[Epoch], 1);
      MagneticModels[Epoch]->CoefficientFileEndDate = MagneticModels[epochs]->CoefficientFileEndDate;
    }
  /*Create Main Field Model for EMM2010*/
  nMaxEMM = MagneticModels[0]->nMax + 1;
  NumTerms = ((nMaxEMM + 1) * (nMaxEMM + 2) / 2);
  MagneticModels[epochs - 1] = MAG_AllocateModelMemory(NumTerms);
  MagneticModels[epochs - 1]->nMax = nMaxEMM;
  MagneticModels[epochs - 1]->nMaxSecVar = nMaxEMM;
  MagneticModels[epochs - 1]->epoch = MagneticModels[0]->epoch + epochs - 1;
  MAG_AssignMagneticModelCoeffs(MagneticModels[epochs - 1], MagneticModels[epochs], MagneticModels[epochs - 1]->nMax, MagneticModels[epochs - 1]->nMaxSecVar);
  /*Allocate Memory for TimedMagneticModel*/

  nMax = MagneticModels[epochs]->nMax;
  NumTerms = ((nMax + 1) * (nMax + 2) / 2);
  TimedMagneticModel = MAG_AllocateModelMemory(NumTerms); /* For storing the time modified WMM Model parameters */
  for(int i = 0; i < epochs; i++) 
    if(MagneticModels[i] == NULL || TimedMagneticModel == NULL)
      {
	MAG_Error(2);
      }
  MAG_SetDefaults(&Ellip, &Geoid); /* Set default values and constants */
  /* Set EGM96 Geoid parameters */
  Geoid.GeoidHeightBuffer = GeoidHeightBuffer;
  Geoid.Geoid_Initialized = 1;
  setDate(date);
  setLonLat(lon,lat);
  setHeight(h);
  setEM();
};


void EMM_Model::setEM()
  {
    int index;
    int Epoch = ((int) UserDate.DecimalYear - MagneticModels[0]->epoch);
    if(Epoch < 0) Epoch = 0;
    if(Epoch > epochs - 1) Epoch = epochs - 1;
    if(LoadedEpoch != Epoch)
      {
	MagneticModels[epochs]->epoch = MagneticModels[Epoch]->epoch;
	MAG_AssignMagneticModelCoeffs(MagneticModels[epochs], MagneticModels[Epoch], MagneticModels[Epoch]->nMax, MagneticModels[Epoch]->nMaxSecVar);
	if(Epoch < epochs - 1)
	  {
	    for(int i = 0; i < 16; i++)
	      {
                        index = 16 * 17 / 2 + i;
                        MagneticModels[epochs]->Secular_Var_Coeff_G[index] = 0;
                        MagneticModels[epochs]->Secular_Var_Coeff_H[index] = 0;
	      }
	  }
                LoadedEpoch = Epoch;
      }
    MAG_GeodeticToSpherical(Ellip, CoordGeodetic, &CoordSpherical); /*Convert from geodetic to Spherical Equations: 17-18, WMM Technical report*/
    MAG_TimelyModifyMagneticModel(UserDate, MagneticModels[epochs], TimedMagneticModel); /* Time adjust the coefficients, Equation 19, WMM Technical report */
    MAG_Geomag(Ellip, CoordSpherical, CoordGeodetic, TimedMagneticModel, &GeoMagneticElements); /* Computes the geoMagnetic field elements and their time change*/
    MAG_CalculateGridVariation(CoordGeodetic, &GeoMagneticElements);
  };
