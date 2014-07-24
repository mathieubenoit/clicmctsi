/*
 * ElectricField.h
 *
 *  Created on: Aug 1, 2012
 *      Author: mbenoit
 */
#include "TH2D.h"
#include "TFile.h"
#ifndef ELECTRICFIELD_H_
#define ELECTRICFIELD_H_
/// This class sets the electric field inside the sensor
class ElectricField {
public:
	/** This is the constructor. In this field they will ask for bias Voltage, depletion Voltage, detector thickness, doTCAd and a filename.
	 *  If bool is false the first three values will be used to use our own uniform magnetic field. If true then we use the filename to point
	 *  to EF.root which is a more complicated electric field obtained from a sophisticated program. So we use an electricfield from EF.root as input
	 */
        ElectricField(double biasVoltage,double depletionVoltage,double detectorThickness, bool doTCAD,TString filename);
    /// destructor
        virtual ~ElectricField();
    /// This function returns the electric field in a certain position.
        double* GetElectricField(double x, double y, double z);
        TH2D *Ex;
        TH2D *Ey;
    /// This boolean determines whether you use the electric field from EF.root or the uniform one we made made ourself
        bool doTCAD;

private :
        double *E;
        double pitch;
        double y_;
        double biasVoltage;
        double depletionVoltage;
        double detectorThickness;
        TFile *f;

};

#endif /* ELECTRICFIELD_H_ */

