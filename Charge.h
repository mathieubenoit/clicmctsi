/*
 * Charge.h
 *
 *  Created on: Jul 30, 2012
 *      Author: tjanssen
 */

#include "Constant.h"
#include <math.h>
#include "TMath.h"
#include <cstdlib>
#include <fstream>
#include <vector>
#ifndef CHARGE_H_
#define CHARGE_H_
#include <iostream>
#include "ElectricField.h"
#include "TRandom3.h"
#include "PixelGeometry.h"

/// In this class the transport of each charge is being simulated
class Charge {
//In protected all variables and constants of a arbitrary charge are defined.
//These values can be used in all functions of the class Charge
protected:

        double x,y,z,t;
        double dx, dy, dz;

        double *E_eff;
        double rh;
        double theta_L;

        double Error;
        double Diffusion_Coefficent;                                    // Default Diffusion coefficient in a semiconductor
        double Default_Mobility;
        double Mobility;                                // Current mobility in [cm2]/[V][s]
        double *Trapping_times;
        double *DeTrapping_times;
        double sum;
        int NumberOfTraps;
        double *Probability;
        int Trap;
        double Random_number;
        double Trap_parameter;
        double Time_before_trapped;
        double Time_before_Detrapped;
        double Prob;
        double *Etemp;
        double *Ef;

        // Mobility related constants
        double alpha;
        double tnom;
        double theta;
        double beta;
        double vsat;
        double mob_0;
        double betaexp;
        double beta_0;
        double TL;

        int Q;
        double sigma;
        double t_0;
        bool TrapState;
        int *Collection_Electrode;
        bool isCollected;
        //Pointers from classes PixelGeometry and ElectricField defined to use their properties.
        PixelGeometry *Geom;
        ElectricField *E;
        //Random Generator of root included (this random generator has a lifetime of >10^200 and is extremely random
        TRandom3 *R;




//In public the behavior of Charge is determined. Functions here compute for instance its Mobility,
//Transport, electric field etc... Also functions that give access to certain values are defined here.
public:
    /// Constructor
        Charge();
    /// destructor
        virtual ~Charge();
    /// This function sets the electricfield, generated from the class ElectricField
        void SetField(ElectricField *E){this->E=E;}
    /// This functions set the geometry of the sensor, generated from the class PixelGeometry
        void SetGeometry(PixelGeometry *Geom){this->Geom=Geom;}
    /// This function set the random generator from ROOT
        void SetRandomGenerator(TRandom3 *R){this->R=R;}
    /// This function set the position of a charge
        void SetPosition(double x, double y, double z);
    /// This function computes the mobility for a given electric field
        void ComputeMobility(double E);
    /// This function returns the current mobility
        double GetMobility(){return Mobility;}
    /// This function returns the Mobility from a certain position in the sensor
        double GetMobility(double x,double y,double z);
    /// This function returns the magnitude of the electric field
        double EFieldMagnitude(double Ex,double Ey=0,double Ez=0){return sqrt(Ex*Ex+Ey*Ey+Ez*Ez); }
    /// This function set the temperature in the sensor
        void SetTemperature(double T){TL=T;}
    /// This function set the trapping times and the number of traps due to the impurities of the material of the sensor
        void SetTrap(double *Trapping_Times, int NTraps);
    /// If a charge is trapped it can also detrap. This function set the detrapping times for each corresponding trap
        void SetDeTrap(double *DeTrapping_Times, int NTraps);
    /// This function computes the time before a charge will be trapped
        void SetTimeBeforeTrapped();
    /// This function computes the time before a charge will leave the trap
        void SetTimeBeforeDeTrapped();
    /// This function checks if a charge is collected. Or in other words has left the sensor
        void Collection();
    /// This function computes the movement of a charge due to diffusion at each time step
        void DoDiffusionStep(double dt);
    /** This function simulates the transport of a charge at each time step. The simulation can be
     *  done with magnetic field (bool true) or without magnetic field (bool false)
     */
        void DoTimeStep(double dt, bool doBField);
    /// This function returns the current time
        double GetTime(){return t;}
    /// This function returns the time before a charge will be trapped
        double GetTimeToTrap(){return Time_before_trapped;}
    /// This functiong returns the time before a charge will leave the trap
        double GetTimeToDetrap(){return Time_before_Detrapped;}
    /** This function returns the trapstate of a charge (so it tells us whether
     *  a charge is trapped or not.
     */
        bool GetIsTrapped(){return this->TrapState;}
    /** This function returns isCollected. In other words it tells us whether
     *  a charge is collected or not.
     */
        bool GetIsCollected(){return isCollected;}
    /// This function returns in which trap a charge is trapped
        int GetTrap(){return Trap+1;}
    /// This function calculates the movement of a charge due to the electric field
        void RKF5Integration(double dt);
    /** When the magnetic field is switched a Lorentz force will influence the charges.
    *  We describe this Lorentz force as an effective electric field. This function
    *  takes this contribution into account and will then do the same as the previous
    *  function.
    */
        void RKF5IntegrationBField(double dt);
    /** These three functions returns the steplength in a certain time step a charge
     *  will make due to diffusion for respectively x, y and z.
     */
        double GetDiffusionStepX(){return dx;}
        double GetDiffusionStepY(){return dy;}
        double GetDiffusionStepZ(){return dz;}
    /** These three function determine the current position of a charge for
     *  respectively x, y and z.
     */
        double GetX(){return x;}
        double GetY(){return y;}
        double GetZ(){return z;}
    /** This function calculates the effective electric field (so including the magnetic field)
     *  at a certain position inside the sensor.
     */
        double* GetElectricEffectiveField(double x, double y, double z);
    /// This function calculates the effective electric field for a given electric field.
        double* GetElectricEffectiveField(double E);
    /// This function returns the Lorentz angle. This angle is the result of the magnetic field
        double GetLA(){return theta_L;}
    /// These three functions return the components of the effective electric field in respectively x, y and z
        double GetEffX(){return E_eff[0];}
        double GetEffY(){return E_eff[1];}
        double GetEffZ(){return E_eff[2];}
    /// This function returns the error in the movement of a charge.
        double GetError(){return Error;}

        std::vector<double> vX;
        std::vector<double> vY;
        std::vector<double> vZ;
        std::vector<double> vT;
};



#endif /* CHARGE_H_ */
