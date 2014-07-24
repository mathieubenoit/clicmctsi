/*
 * Interaction.h
 *
 *  Created on: Aug 1, 2012
 *      Author: tjanssen
 */

#ifndef INTERACTION_H_
#define INTERACTION_H_
#include "TH1I.h"
#include "TApplication.h"
#include "Electron.h"
#include "Hole.h"

/** This class combines all the classes and all charges. It computes the total interaction, transport for all
 *  particles.
 */
class Interaction {
  public:
   /** The constructor.
	*  Mean and RMS are set to generate the initial parameters for
	*  the Gaussian distribution to give all particles
	*  the initial position. The boolean is to generate a Gaussian
	*  position if true or make z uniform while rotating
	*  over the x axis under an angle theta if false.
	*  After the choice of the distribution one can set the angle which
	*  will only be used if the boolean is false. The second boolean is
	*  to switch the magnetic field (true) on or off (false)
	*/
        Interaction(double Mean[3], double RMS[3], bool Distr_Gaus, double theta,
		    bool doBField, int Nparticles, TRandom3 *R,ElectricField *E,
		    PixelGeometry *Geom);
	
	/// The desructor
        virtual ~Interaction();
	
	/** Generate an initial position in the sensor for a mean position with a Gaussian
	 *  distribution with a given RMS.
	 */
        double* GenerateGaussianPosition(double mean[3], double RMS[3]);

	/** Generate an initial position in the sensor for a mean position with a
	 *  uniform distribution in z and a Gaussian distribution with a given RMS in x and y.
	 */
        double* GenerateUniformPosition(double mean[3], double RMS[3]);
    /// This function set the random generator from ROOT
        void SetRandomGenerator(TRandom3 *R){this->rand=R;}
    /// This function draws graphs and histograms of all the results.
        void MakeGraphs();
    /// This function takes a step into time
        void DoATimeStep(double dt);
    /// This function returns the current time
        double GetTime(){return t;}
    /** Each time step the charges will move in the sensor. For each step
     *  in space we calculate the corresponding error. This function
     *  calculates the average error of all particles taking a step in space.
     */
        void ComputeAverageError();
    /** If the error becomes to large we want the time step to decrease to reduce the error
     *  this function makes sure that happens.
     */
        void AdjustTimeStep();
    /** We want to run the simulation till all particles are collected or the particles are trapped
     *  for more then a certain time. This function makes sure this is the case.
     */
        bool SimulationNotDone();

	/// Get electron number i from the list of electrons
        Electron* GetElectron(int i){return electrons[i];}
    /// Get hole number i from the list of hole
        Hole* GetHole(int i){return holes[i];}


private :
        Electron **electrons; //< An array of pointers to Electrons
        Hole **holes;		  //< An array of pointers to Holes
        ElectricField *E;
        PixelGeometry *Geom;
        TRandom3 *rand;
        Electron *anElectron;
        Hole *aHole;
        bool Bfield;
        int Nparticles;
        int NumberOfTraps;
        int counter;
        double thickness;
        double u, v, w;
        double u_, v_, w_;
        double alpha;
        double t;
        double Meanxe, Meanye, Meanze;
        double RMSxe, RMSye, RMSze;
        double Meanxh, Meanyh, Meanzh;
        double RMSxh, RMSyh, RMSzh;
        double sumxe, sumye, sumze;
        double sumxe2, sumye2, sumze2;
        double sumxh, sumyh, sumzh;
        double sumxh2, sumyh2, sumzh2;
        double timestep;
        double TargetError;
        double downscaling_factor;
        double upscaling_factor;
        double Error;
        double *Ef;
        double Time_stop;
        double *Coord;

        double *Hole_Trapping_Times;
        double *Electron_Trapping_Times;

        double *Hole_DeTrapping_Times;
        double *Electron_DeTrapping_Times;

};

#endif /* INTERACTION_H_ */

