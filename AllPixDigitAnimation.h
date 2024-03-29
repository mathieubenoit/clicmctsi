// @(#)codeTDv2Final/CTDR_Animation :$Name:  $:$Id: AllPixDigitAnimation.h 49 2014-07-11 12:10:49Z mbenoit $
// Author: Mathieu Benoit   24/11/06

//////////////////////////////////////////////////////////////////////////
// Source file : class CTDR_Animation					//
// This class is used to produce a root file containing the tracks of  	//
// the charge elements simulated. It creates a geometry, defined in the //
// the constructor, and everything needed to draw image and animation   //
// of the simulated event in the geometry . The root file can be read  	//
// using cint interpreter, and the geometry and the charge tracks can	//
// then be manipulated using functions provided with root's geometry  	//
// packages , see TGeoManager and TGeoTracks classes.	 		//
//////////////////////////////////////////////////////////////////////////

#ifndef ALLPIXDIGITANIMATION_H_DEF
#define ALLPIXDIGITANIMATION_H_DEF

#include <vector>
#include "TGeoManager.h"
#include "TVirtualGeoTrack.h"
#include "TFile.h"
#include "TGeoBBox.h"
#include "TString.h"

/// This class animates the transport of all particles in the sensor
class AllPixDigitAnimation
{
      
      private :
             TGeoManager *Geo; // The geometry manager
             TFile *f; // The file where the objects are written
             int Ntracks; // The number of tracks
             TGeoMedium *medium; // the material
             TGeoVolume *top ; // the box
             TGeoBBox *PixBox; // The pixels
             TGeoVolume *PixVolume; //the pixels volume
             
             
             //double Lx,Ly,Lz;  // size of the box
	     	 int trackid;
	    	 int nx,ny;
	    	 double Lz;
	    	 double shiftx,shifty;
	   		 double pitchx,pitchy;
         	 Int_t MyPalette[100];
         	 
         	 double emax;
         	 
         	 double z_hit;
             
      public :
             AllPixDigitAnimation(int nx, int ny, double lz, double pitchx, double pitchy,int nHits, int eventid);
             void AddTrack(std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double> t,double sigma,double energy, int color);
	     	 double modulo(double a , double b);
	    	 void FirePixel(int i , int j);	
		 void SubThresholdPixel(int i , int j);
	         void SetShift(double x, double y){
	     		shiftx=x;
	     		shifty=y;
	     };
	     

	         void ViewTracks(){Geo->DrawTracks();};
             ~AllPixDigitAnimation();
      
};
      
#endif            
             
      
