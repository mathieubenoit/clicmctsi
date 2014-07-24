/*
 * ElectricField.cpp
 *
 *  Created on: Aug 1, 2012
 *      Author: mbenoit
 */
#include "ElectricField.h"
#include "TCanvas.h"
#include "Constant.h"
#include "TFile.h"
#include "TDirectory.h"
#include <fstream>
#include <iostream>
#include <math.h>
using namespace std;
/** In this class the electricfield component in the z direction is calculated. It makes sure the partcle will
 *  travel in the right direction. This will be shown in the class PixelGeometry.
 */
	ElectricField::ElectricField(double biasVoltage,double depletionVoltage,double detectorThickness,bool doTCAD, TString filename) {

        E = new double[3];
        this->biasVoltage=biasVoltage;
        this->depletionVoltage=depletionVoltage;
        this->detectorThickness=detectorThickness;
        this->doTCAD=doTCAD;
        Ex=0;
        Ey=0;
        pitch=25*um;
    	y_=0;
    	f = new TFile(filename,"open");
   /** If doTCad is false than we will use our own defined electric field. Otherwise we will get a more
    *  sophisticated electric field from the Ef.root
    */
        if(doTCAD==true)
        {
                f->GetObject("Ex",Ex);
                f->GetObject("Ey",Ey);
        }
}

double* ElectricField::GetElectricField(double x, double y, double z){

        double zp= detectorThickness-(z+detectorThickness/2);
		y_=((y+pitch/2)-pitch*floor((y+pitch/2)/pitch));

        if(doTCAD==false){
        E[0]=0;
        E[1]=0;
        if (TMath::Abs(z)>detectorThickness/2) E[2]=0;
        else
        E[2]=(biasVoltage-depletionVoltage)/detectorThickness+(1-zp/detectorThickness)*2*depletionVoltage/detectorThickness;
        }
        else {
            /** For instance the first line here I defined to make sure that the electric field is computed in
             *  all pixels. This algorithm only takes into account 1 pixel. To determine the relative position
             *  at all pixels we can use the same function!
             */
                if (TMath::Abs(z)>detectorThickness/2) {
                    E[0]=0;
                    E[1]=0;
                	E[2]=0;
                }
                else {
                	E[0]=0;
                	E[1]=Ey->Interpolate(y_/um,zp/um)*V/cm;
                	E[2]=Ex->Interpolate(y_/um,zp/um)*V/cm;
                }
        }
        return E;
}
/// In the destructor I deleted the memory. This was not done yet and it results into memory leaks.
	ElectricField::~ElectricField()
	{
        delete E;
        delete f;
        // TODO Auto-generated destructor stub
	}

