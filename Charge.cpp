 /*
 * Charge.cpp
 *
 *  Created on: Jul 30, 2012
 *      Author: tjanssen
 */

#include "Charge.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFitResult.h"
#include "TPaveText.h"
#include "TF1.h"
#include "TAxis.h"
#include <math.h>
using namespace std;
/** In the constructor values of protected can be initialized. From than on they are global for all functions,
 *  In the class Charge. As one can see the pointers in protected are allocated in the free space here. In order
 *  to prevent memory leaks we have to delete this free space in the destructor.
 */
	Charge::Charge()
	{
     // TODO Auto-generated constructor stub
    	TrapState = false;
    	Time_before_trapped = 0;
    	Time_before_Detrapped = 0;
    	Geom = 0;
    	E = 0;
    	dx=dy=dz=0;
    	Collection_Electrode = new int[2];
    	isCollected = false;
    	Etemp = new double[3];
    	Ef = new double[3];
    	t = 0;
    	t_0 = 0;
    	sigma=0;
    	theta_L=0;
    	E_eff = new double[3];
    	Error=0;
	}

/** This function computes mobility. The constants used here are given in Hole and Electron.
 *  In the class interaction this function will be called through a pointer from electron or Hole.
 *  That way the right values for the constants are used in order to calculate to mobility.
 */
	void Charge::ComputeMobility(double E)
	{
    	vsat = alpha/(1.+ theta * exp(TL/tnom));
    	beta = beta_0*pow(TL/300.,betaexp);
    	Mobility =      mob_0*pow((1.0/(1. + pow(mob_0*E/vsat,beta) ) ),1./beta);
	}

/** In this function we set the current position of a charge.
 *  This->x=x tells specifically to access x from protected. Probably it will work fine without, but to
 *  be sure we force it to. The function as it says sets the position of a certain charge.
 */
	void Charge::SetPosition(double x, double y, double z)
	{
    	this->x=x;
    	this->y=y;
    	this->z=z;
	}

/** In the silicon semiconductor are certain traps. A charge might fall into a trap with a certain probability.
 *  In this function the parameters to determine in which trap a charge might fall are being set.
 */
	void Charge::SetTrap(double *Trapping_Times, int NTraps)
	{
		NumberOfTraps=NTraps;

		this->Trapping_times      =   new double[NumberOfTraps];
		this->Probability         =   new double[NumberOfTraps];
		this->Trapping_times      =   Trapping_Times;
	}

/** In here the timebeforetrapped is calculated. The name speaks for itself. In this function a charge will,
 *  fall in a certain trap with a probability. First the probability to fall in a certain trap is being calculated,
 *  If there are for instance two traps with the same trapping times than the chance to fall in trap1=trap2=0.5
 *  After we determined the chance to fall in a certain trap we generate a random number between 0 and 1.
 *  This value determines in which trap the charge will fall. After that the timebefore the charge will be
 *  trapped is being calculated using Monte Carlo. The charge will only be trapped if it will not be collected in the main time.
 *  So if timebeforecollected<timebeforetrapped. The charge will not be trapped.
 */
	void Charge::SetTimeBeforeTrapped()
	{
		/// Calculating sum of 1 over trapping_times
		sum=0;
		for(int i=0; i<NumberOfTraps; ++i)
		{
            sum += 1./Trapping_times[i];
		}

		/// Calculating trap_parameter
		Trap_parameter = 1./sum;
		/// Calculating probability to fall in certain trap j
		for(int j=0; j<NumberOfTraps; ++j)
		{
            Probability[j] = (1./Trapping_times[j])/sum;
		}

		Random_number = R->Uniform();

		for(int l=0; l<NumberOfTraps; ++l)
		{
            if(l==0)
            {
                    if(Random_number <= Probability[0])
                            Trap=0;
            }
            else
            {
                    if(Random_number <= Probability[l]+Probability[l-1] && Random_number > Probability[l-1])
                            Trap=l;
            }
		}

		Prob = R->Uniform();
		Time_before_trapped = -Trap_parameter*log(1-Prob);
	}
/** When a charge is trapped it will also detrap after some time. In here the parameters to set the
 *  timebeforedetrapped are set.
 */
	void Charge::SetDeTrap(double *DeTrapping_Times, int NTraps)
	{
        NumberOfTraps=NTraps;

        this->DeTrapping_times    =   new double[NumberOfTraps];
        this->DeTrapping_times    =   DeTrapping_Times;

	}

/** Just like in settimebeforetrapped here timebeforedetrapped is being calculated using Monte Carlo.
 *  In this function we don't have to calculate in which trap it falls. We already in which trap the charge is.
 *  So we can just just the detrapping time of the trap L.
 */
	void Charge::SetTimeBeforeDeTrapped()
	{
		Prob = R->Uniform();
		Time_before_Detrapped = -DeTrapping_times[Trap]*log(1-Prob);
	}

/** In this function a timestep is being done for each charge. In here the transport is being done for each
 *  individual particle.
 */
	void Charge::DoTimeStep(double dt, bool doBField)
	{
        	vX.push_back(x);
        	vY.push_back(y);
        	vZ.push_back(z);
        	vT.push_back(t);

        if (TrapState==false && isCollected==false)
     	 {
             	 if(doBField==false)
            	 	  RKF5Integration(dt);
             	 else RKF5IntegrationBField(dt);
             /// Here we calculate the step in space a particle will take due to diffusion.
             	 DoDiffusionStep(dt);
             /// After each time step we want to know whether the particle is collected or not with Collection
             	 Collection();
     	 }
     /// if a particle is collected than nothing has to be done with it anymore.
     	 if (isCollected==true)
     	 {

     	 }
     /// If charge is not collected it can either move in space or can be trapped/detrapped again.
     	 else
     	 {
             	 if(Time_before_trapped ==0)
            	 	 SetTimeBeforeTrapped();
             	 if (t_0>= Time_before_trapped)
             	 {
                 	 TrapState = true;
                 	 Time_before_trapped=0;
             	 }
             	 if (TrapState == true)
             	 {
                 	 if(Time_before_Detrapped==0)
                	 	 SetTimeBeforeDeTrapped();
                     	 if(t_0>= Time_before_Detrapped+Time_before_trapped)
                     	 {
                         	 TrapState = false;
                         	 Time_before_Detrapped=0;
                         	 t_0=0;
                     	 }
             	 }
     	 }

     	 t = t + dt;
     	 t_0 = t_0 + dt;
	}

/** In here we check if a charge is collected. First we check if the particle is charge 1 or -1 or in other
 *  word if the particle is an electron or hole. Than we call for the class geometry to check if the pixel
 *  is inside or outside the geometry. The particles will all start inside the geometry. When outside we
 *  define the charge as being collected. If you comment out the cout's the program will tell where the
 *  charges are being collected. To check this in detail see class PixelGeometry.
 */
	void Charge::Collection()
	{
		if (!Geom->IsInside(x,y,z) && isCollected==false)
		{
			/// The cout statements here tells if where the electron/hole is collected.
				Collection_Electrode=Geom->ClosestElectrode(x,y,z);
				if(Q==-1){
					//cout << TString::Format("Electron collected at electrode (%d,%d) \n",Collection_Electrode[0],Collection_Electrode[1]) << endl;
				}

				else {
					//cout << TString::Format("Hole collected at electrode (%d,%d) \n",Collection_Electrode[0],Collection_Electrode[1]) << endl;
				}

             isCollected=true;
		}
	}

/** Here the diffusionstep is being calculated. We consider a diffusion coefficient for a gaussian-shaped distribution.
 *  The mean is taken to be 0 and we calculate the RMS (sigma) with the diffusion coefficient. This coefficient is different for holes
 *  and electrons. This RMS is an approximation but in this algorithm the corrections will probably be neglectable
 *  If I have more time I will however also take more corrections into account.
 */
	void Charge::DoDiffusionStep(double dt)
	{
			double u1, u2, v1, v2;
			sigma = sqrt(2*Diffusion_Coefficent*dt);

            u1 = R->Uniform();
            u2 = R->Uniform();
            v1 = R->Uniform();
            v2 = R->Uniform();
        /// This expression is an easy way to calculate a Gaussion distribution. From this distribution
            dx = sigma*sin(2*pi*u1)*sqrt(-2*log(u2));
            dy = sigma*cos(2*pi*u2)*sqrt(-2*log(u1));
            dz = sigma*sin(2*pi*v1)*sqrt(-2*log(v2));

            x += dx;
            y += dy;
            z += dz;
	}

/** This function transport using Euler integration, for field (Ex,Ey,Ez),
 *  considered constant over time dt. The movement equation are those
 *  of charges in semi-conductors, sx= mu*E*dt;
 */
	void Charge::RKF5Integration(double dt)
	{
      double k1x,k2x,k3x,k4x,k5x,k6x;
      double k1y,k2y,k3y,k4y,k5y,k6y;
      double k1z,k2z,k3z,k4z,k5z,k6z;
      double dx, dy, dz;

      Ef=E->GetElectricField(x,y,z);

      k1x=-Q*GetMobility(x,y,z)*Ef[0]*dt;
      k1y=-Q*GetMobility(x,y,z)*Ef[1]*dt;
      k1z=-Q*GetMobility(x,y,z)*Ef[2]*dt;

      Ef=E->GetElectricField(x+k1x/4,y+k1y/4,z+k1z/4);

      k2x=-Q*GetMobility(x+k1x/4,y+k1y/4,z+k1z/4)*Ef[0]*dt;
      k2y=-Q*GetMobility(x+k1x/4,y+k1y/4,z+k1z/4)*Ef[1]*dt;
      k2z=-Q*GetMobility(x+k1x/4,y+k1y/4,z+k1z/4)*Ef[2]*dt;

      Ef=E->GetElectricField(x+(9./32)*k2x+(3./32)*k1x,y+(9./32)*k2y+(3./32)*k1y,z+(9./32)*k2z+(3./32)*k1z);

      k3x=-Q*GetMobility(x+(9./32)*k2x+(3./32)*k1x,y+(9./32)*k2y+(3./32)*k1y,z+(9./32)*k2z+(3./32)*k1z)*Ef[0]*dt;
      k3y=-Q*GetMobility(x+(9./32)*k2x+(3./32)*k1x,y+(9./32)*k2y+(3./32)*k1y,z+(9./32)*k2z+(3./32)*k1z)*Ef[1]*dt;
      k3z=-Q*GetMobility(x+(9./32)*k2x+(3./32)*k1x,y+(9./32)*k2y+(3./32)*k1y,z+(9./32)*k2z+(3./32)*k1z)*Ef[2]*dt;

      Ef=E->GetElectricField(x-(7200./2197)*k2x+(1932./2197)*k1x+(7296./2197)*k3x,y-(7200./2197)*k2y+(1932./2197)*k1y+(7296./2197)*k3y,z-(7200./2197)*k2z+(1932./2197)*k1z+(7296./2197)*k3z);

      k4x=-Q*GetMobility(x-(7200./2197)*k2x+(1932./2197)*k1x+(7296./2197)*k3x,y-(7200./2197)*k2y+(1932./2197)*k1y+(7296./2197)*k3y,z-(7200./2197)*k2z+(1932./2197)*k1z+(7296./2197)*k3z)*Ef[0]*dt;
      k4y=-Q*GetMobility(x-(7200./2197)*k2x+(1932./2197)*k1x+(7296./2197)*k3x,y-(7200./2197)*k2y+(1932./2197)*k1y+(7296./2197)*k3y,z-(7200./2197)*k2z+(1932./2197)*k1z+(7296./2197)*k3z)*Ef[1]*dt;
      k4z=-Q*GetMobility(x-(7200./2197)*k2x+(1932./2197)*k1x+(7296./2197)*k3x,y-(7200./2197)*k2y+(1932./2197)*k1y+(7296./2197)*k3y,z-(7200./2197)*k2z+(1932./2197)*k1z+(7296./2197)*k3z)*Ef[2]*dt;

      Ef=E->GetElectricField(x-(8)*k2x+(439./216)*k1x+(3680./513)*k3x-(845./4104)*k4x,y-(8)*k2y+(439./216)*k1y+(3680./513)*k3y-(845./4104)*k4y,z-(8)*k2z+(439./216)*k1z+(3680./513)*k3z-(845./4104)*k4z);

      k5x=-Q*GetMobility(x-(8)*k2x+(439./216)*k1x+(3680./513)*k3x-(845./4104)*k4x,y-(8)*k2y+(439./216)*k1y+(3680./513)*k3y-(845./4104)*k4y,z-(8)*k2z+(439./216)*k1z+(3680./513)*k3z-(845./4104)*k4z)*Ef[0]*dt;
      k5y=-Q*GetMobility(x-(8)*k2x+(439./216)*k1x+(3680./513)*k3x-(845./4104)*k4x,y-(8)*k2y+(439./216)*k1y+(3680./513)*k3y-(845./4104)*k4y,z-(8)*k2z+(439./216)*k1z+(3680./513)*k3z-(845./4104)*k4z)*Ef[1]*dt;
      k5z=-Q*GetMobility(x-(8)*k2x+(439./216)*k1x+(3680./513)*k3x-(845./4104)*k4x,y-(8)*k2y+(439./216)*k1y+(3680./513)*k3y-(845./4104)*k4y,z-(8)*k2z+(439./216)*k1z+(3680./513)*k3z-(845./4104)*k4z)*Ef[2]*dt;

      Ef=E->GetElectricField(x+(2)*k2x-(8./27)*k1x-(3544./2565)*k3x-(1859./4104)*k4x-(11./40)*k5x,
      y+(2)*k2y-(8./27)*k1y-(3544./2565)*k3y-(1859./4104)*k4y-(11./40)*k5y,
      z+(2)*k2z-(8./27)*k1z-(3544./2565)*k3z-(1859./4104)*k4z-(11./40)*k5z);

      k6x=-Q*GetMobility(x+(2)*k2x-(8./27)*k1x-(3544./2565)*k3x-(1859./4104)*k4x-(11./40)*k5x,
              y+(2)*k2y-(8./27)*k1y-(3544./2565)*k3y-(1859./4104)*k4y-(11./40)*k5y,
              z+(2)*k2z-(8./27)*k1z-(3544./2565)*k3z-(1859./4104)*k4z-(11./40)*k5z)*Ef[0]*dt;
      k6y=-Q*GetMobility(x+(2)*k2x-(8./27)*k1x-(3544./2565)*k3x-(1859./4104)*k4x-(11./40)*k5x,
              y+(2)*k2y-(8./27)*k1y-(3544./2565)*k3y-(1859./4104)*k4y-(11./40)*k5y,
              z+(2)*k2z-(8./27)*k1z-(3544./2565)*k3z-(1859./4104)*k4z-(11./40)*k5z)*Ef[1]*dt;
      k6z=-Q*GetMobility(x+(2)*k2x-(8./27)*k1x-(3544./2565)*k3x-(1859./4104)*k4x-(11./40)*k5x,
              y+(2)*k2y-(8./27)*k1y-(3544./2565)*k3y-(1859./4104)*k4y-(11./40)*k5y,
              z+(2)*k2z-(8./27)*k1z-(3544./2565)*k3z-(1859./4104)*k4z-(11./40)*k5z)*Ef[2]*dt;

      dx=((16./135)*k1x+(6656./12825)*k3x+(28561./56430)*k4x-(9./50)*k5x+(2./55)*k6x);
      dy=((16./135)*k1y+(6656./12825)*k3y+(28561./56430)*k4y-(9./50)*k5y+(2./55)*k6y);
      dz=((16./135)*k1z+(6656./12825)*k3z+(28561./56430)*k4z-(9./50)*k5z+(2./55)*k6z);

      double Ex,Ey,Ez;
      Ex=((1./360)*k1x-(128./4275)*k3x-(2197./75240)*k4x-(1./50)*k5x+(2./55)*k6x);
      Ey=((1./360)*k1y-(128./4275)*k3y-(2197./75240)*k4y-(1./50)*k5y+(2./55)*k6y);
      Ez=((1./360)*k1z-(128./4275)*k3z-(2197./75240)*k4z-(1./50)*k5z+(2./55)*k6z);
      Error=sqrt(Ex*Ex+Ey*Ey+Ez*Ez);

      x+=dx;
      y+=dy;
      z+=dz;
	}
/// This function calculates the electric field component in the x and y direction if the magnetic field is turned on.
	double *Charge::GetElectricEffectiveField(double x, double y, double z)
	{
        if(Geom->IsInside(x,y,z)==true)
        {
                Ef=E->GetElectricField(x,y,z);
                double B[3]={0};
                double B_length=4*T;
                B[0]=4*T;
                B[1]=sqrt(pow(B_length,2)-pow(B[0],2));
                E_eff[0]=((B[1]*sqrt(pow(B[0],2) + pow(B[1],2))*Ef[2]*sqrt(pow(Ef[0],2) + pow(Ef[1],2) + pow(Ef[2],2))*GetMobility(x,y,z)*rh)/sqrt(-2*B[0]*B[1]*Ef[0]*Ef[1] + pow(B[1],2)*(pow(Ef[0],2) + pow(Ef[2],2)) + pow(B[0],2)*(pow(Ef[1],2) + pow(Ef[2],2))))+Ef[1];
                E_eff[1]=-((B[0]*sqrt(pow(B[0],2) + pow(B[1],2))*Ef[2]*sqrt(pow(Ef[0],2) + pow(Ef[1],2) + pow(Ef[2],2))*GetMobility(x,y,z)*rh)/sqrt(-2*B[0]*B[1]*Ef[0]*Ef[1] + pow(B[1],2)*(pow(Ef[0],2) + pow(Ef[2],2)) + pow(B[0],2)*(pow(Ef[1],2) + pow(Ef[2],2))))+Ef[0];
                E_eff[2]=((sqrt(pow(B[0],2) + pow(B[1],2))*(B[1]*Ef[0] - B[0]*Ef[1])*sqrt(pow(Ef[0],2) + pow(Ef[1],2) + pow(Ef[2],2))*GetMobility(x,y,z)*rh)/sqrt(-2*B[0]*B[1]*Ef[0]*Ef[1] + pow(B[1],2)*(pow(Ef[0],2) + pow(Ef[2],2)) + pow(B[0],2)*(pow(Ef[1],2) + pow(Ef[2],2))))+Ef[2];
                }
        else
        {
                E_eff[0]=0;
                E_eff[1]=0;
                E_eff[2]=0;
        }

        return E_eff;
	}

/** This function does exectly the same thing as the previous one, but will not make use of the interaction.
 *  This one I use to make a plot of the Ex,y component vs Ez and the lorentz angle vs Ez. I made a different
 *  function for this so I would not affect the interaction with it. Besides it does not really change the performance
 *  so it was safer to do so.
 */
	double *Charge::GetElectricEffectiveField(double E)
	{
		Ef[0]=0;
		Ef[1]=0;
		Ef[2]=E;
		double B[3]={0};
		double B_length=4*T;
		B[0]=0;
		B[1]=sqrt(pow(B_length,2)-pow(B[0],2));
		this->ComputeMobility(Ef[2]);
		theta_L=atan(Mobility*rh*B_length);
		E_eff[0]=-((B[1]*sqrt(pow(B[0],2) + pow(B[1],2))*Ef[2]*sqrt(pow(Ef[0],2) + pow(Ef[1],2) + pow(Ef[2],2))*Mobility*rh)/sqrt(-2*B[0]*B[1]*Ef[0]*Ef[1] + pow(B[1],2)*(pow(Ef[0],2) + pow(Ef[2],2)) + pow(B[0],2)*(pow(Ef[1],2) + pow(Ef[2],2))))+Ef[0];
		E_eff[1]=((B[0]*sqrt(pow(B[0],2) + pow(B[1],2))*Ef[2]*sqrt(pow(Ef[0],2) + pow(Ef[1],2) + pow(Ef[2],2))*Mobility*rh)/sqrt(-2*B[0]*B[1]*Ef[0]*Ef[1] + pow(B[1],2)*(pow(Ef[0],2) + pow(Ef[2],2)) + pow(B[0],2)*(pow(Ef[1],2) + pow(Ef[2],2))))+Ef[1];
		E_eff[2]=((sqrt(pow(B[0],2) + pow(B[1],2))*(B[1]*Ef[0] - B[0]*Ef[1])*sqrt(pow(Ef[0],2) + pow(Ef[1],2) + pow(Ef[2],2))*Mobility*rh)/sqrt(-2*B[0]*B[1]*Ef[0]*Ef[1] + pow(B[1],2)*(pow(Ef[0],2) + pow(Ef[2],2)) + pow(B[0],2)*(pow(Ef[1],2) + pow(Ef[2],2))))+Ef[2];

		return E_eff;
	}


/** this integration is nearly the same as the one before. This one is being used if the magnetic field is
 *  turned on. This changes the electric field so we call a different field here. The rest is exactly the same
 *  We defined a different function for this to avoid errors. Besides it doesn't change the performance at all.
 */
	void Charge::RKF5IntegrationBField(double dt)
	{
      double k1x,k2x,k3x,k4x,k5x,k6x;
      double k1y,k2y,k3y,k4y,k5y,k6y;
      double k1z,k2z,k3z,k4z,k5z,k6z;
      double dx,dy,dz;

      E_eff=this->GetElectricEffectiveField(x,y,z);

      k1x=-Q*GetMobility(x,y,z)*E_eff[0]*dt;
      k1y=-Q*GetMobility(x,y,z)*E_eff[1]*dt;
      k1z=-Q*GetMobility(x,y,z)*E_eff[2]*dt;

      E_eff=this->GetElectricEffectiveField(x+k1x/4,y+k1y/4,z+k1z/4);

      k2x=-Q*GetMobility(x+k1x/4,y+k1y/4,z+k1z/4)*E_eff[0]*dt;
      k2y=-Q*GetMobility(x+k1x/4,y+k1y/4,z+k1z/4)*E_eff[1]*dt;
      k2z=-Q*GetMobility(x+k1x/4,y+k1y/4,z+k1z/4)*E_eff[2]*dt;

      E_eff=this->GetElectricEffectiveField(x+(9./32)*k2x+(3./32)*k1x,y+(9./32)*k2y+(3./32)*k1y,z+(9./32)*k2z+(3./32)*k1z);

      k3x=-Q*GetMobility(x+(9./32)*k2x+(3./32)*k1x,y+(9./32)*k2y+(3./32)*k1y,z+(9./32)*k2z+(3./32)*k1z)*E_eff[0]*dt;
      k3y=-Q*GetMobility(x+(9./32)*k2x+(3./32)*k1x,y+(9./32)*k2y+(3./32)*k1y,z+(9./32)*k2z+(3./32)*k1z)*E_eff[1]*dt;
      k3z=-Q*GetMobility(x+(9./32)*k2x+(3./32)*k1x,y+(9./32)*k2y+(3./32)*k1y,z+(9./32)*k2z+(3./32)*k1z)*E_eff[2]*dt;

      E_eff=this->GetElectricEffectiveField(x-(7200./2197)*k2x+(1932./2197)*k1x+(7296./2197)*k3x,y-(7200./2197)*k2y+(1932./2197)*k1y+(7296./2197)*k3y,z-(7200./2197)*k2z+(1932./2197)*k1z+(7296./2197)*k3z);

      k4x=-Q*GetMobility(x-(7200./2197)*k2x+(1932./2197)*k1x+(7296./2197)*k3x,y-(7200./2197)*k2y+(1932./2197)*k1y+(7296./2197)*k3y,z-(7200./2197)*k2z+(1932./2197)*k1z+(7296./2197)*k3z)*E_eff[0]*dt;
      k4y=-Q*GetMobility(x-(7200./2197)*k2x+(1932./2197)*k1x+(7296./2197)*k3x,y-(7200./2197)*k2y+(1932./2197)*k1y+(7296./2197)*k3y,z-(7200./2197)*k2z+(1932./2197)*k1z+(7296./2197)*k3z)*E_eff[1]*dt;
      k4z=-Q*GetMobility(x-(7200./2197)*k2x+(1932./2197)*k1x+(7296./2197)*k3x,y-(7200./2197)*k2y+(1932./2197)*k1y+(7296./2197)*k3y,z-(7200./2197)*k2z+(1932./2197)*k1z+(7296./2197)*k3z)*E_eff[2]*dt;

      E_eff=this->GetElectricEffectiveField(x-(8)*k2x+(439./216)*k1x+(3680./513)*k3x-(845./4104)*k4x,y-(8)*k2y+(439./216)*k1y+(3680./513)*k3y-(845./4104)*k4y,z-(8)*k2z+(439./216)*k1z+(3680./513)*k3z-(845./4104)*k4z);

      k5x=-Q*GetMobility(x-(8)*k2x+(439./216)*k1x+(3680./513)*k3x-(845./4104)*k4x,y-(8)*k2y+(439./216)*k1y+(3680./513)*k3y-(845./4104)*k4y,z-(8)*k2z+(439./216)*k1z+(3680./513)*k3z-(845./4104)*k4z)*E_eff[0]*dt;
      k5y=-Q*GetMobility(x-(8)*k2x+(439./216)*k1x+(3680./513)*k3x-(845./4104)*k4x,y-(8)*k2y+(439./216)*k1y+(3680./513)*k3y-(845./4104)*k4y,z-(8)*k2z+(439./216)*k1z+(3680./513)*k3z-(845./4104)*k4z)*E_eff[1]*dt;
      k5z=-Q*GetMobility(x-(8)*k2x+(439./216)*k1x+(3680./513)*k3x-(845./4104)*k4x,y-(8)*k2y+(439./216)*k1y+(3680./513)*k3y-(845./4104)*k4y,z-(8)*k2z+(439./216)*k1z+(3680./513)*k3z-(845./4104)*k4z)*E_eff[2]*dt;

      E_eff=this->GetElectricEffectiveField(x+(2)*k2x-(8./27)*k1x-(3544./2565)*k3x-(1859./4104)*k4x-(11./40)*k5x,
      y+(2)*k2y-(8./27)*k1y-(3544./2565)*k3y-(1859./4104)*k4y-(11./40)*k5y,
      z+(2)*k2z-(8./27)*k1z-(3544./2565)*k3z-(1859./4104)*k4z-(11./40)*k5z);

      k6x=-Q*GetMobility(x+(2)*k2x-(8./27)*k1x-(3544./2565)*k3x-(1859./4104)*k4x-(11./40)*k5x,
              y+(2)*k2y-(8./27)*k1y-(3544./2565)*k3y-(1859./4104)*k4y-(11./40)*k5y,
              z+(2)*k2z-(8./27)*k1z-(3544./2565)*k3z-(1859./4104)*k4z-(11./40)*k5z)*E_eff[0]*dt;
      k6y=-Q*GetMobility(x+(2)*k2x-(8./27)*k1x-(3544./2565)*k3x-(1859./4104)*k4x-(11./40)*k5x,
              y+(2)*k2y-(8./27)*k1y-(3544./2565)*k3y-(1859./4104)*k4y-(11./40)*k5y,
              z+(2)*k2z-(8./27)*k1z-(3544./2565)*k3z-(1859./4104)*k4z-(11./40)*k5z)*E_eff[1]*dt;
      k6z=-Q*GetMobility(x+(2)*k2x-(8./27)*k1x-(3544./2565)*k3x-(1859./4104)*k4x-(11./40)*k5x,
              y+(2)*k2y-(8./27)*k1y-(3544./2565)*k3y-(1859./4104)*k4y-(11./40)*k5y,
              z+(2)*k2z-(8./27)*k1z-(3544./2565)*k3z-(1859./4104)*k4z-(11./40)*k5z)*E_eff[2]*dt;

      dx=((16./135)*k1x+(6656./12825)*k3x+(28561./56430)*k4x-(9./50)*k5x+(2./55)*k6x);
      dy=((16./135)*k1y+(6656./12825)*k3y+(28561./56430)*k4y-(9./50)*k5y+(2./55)*k6y);
      dz=((16./135)*k1z+(6656./12825)*k3z+(28561./56430)*k4z-(9./50)*k5z+(2./55)*k6z);

      double Ex,Ey,Ez;
      Ex=((1./360)*k1x-(128./4275)*k3x-(2197./75240)*k4x-(1./50)*k5x+(2./55)*k6x);
      Ey=((1./360)*k1y-(128./4275)*k3y-(2197./75240)*k4y-(1./50)*k5y+(2./55)*k6y);
      Ez=((1./360)*k1z-(128./4275)*k3z-(2197./75240)*k4z-(1./50)*k5z+(2./55)*k6z);
      Error=sqrt(Ex*Ex+Ey*Ey+Ez*Ez);
      //cout << Error << endl;


      //cout << "[Transport] " << TString::Format("Q=%d x=%f y=%f z=%f ",Q,x/um,y/um,z/um) << endl;
      x+=dx;
      y+=dy;
      z+=dz;
      //cout << "[Transport after step] " << TString::Format("Q=%d x=%f y=%f z=%f ",Q,x/um,y/um,z/um) << endl;

	}
/// This function returns the Mobility computed before. Getmobility is used in the intergration.
	double Charge::GetMobility(double x,double y,double z)
	{
		Etemp=E->GetElectricField(x,y,z);
		this->ComputeMobility(EFieldMagnitude(Etemp[0],Etemp[1],Etemp[2]));
		return Mobility;
	}
/// In the destructor I delete all memory alloced in the free space. So no memory leaks will occur.
	Charge::~Charge()
	{
		delete Collection_Electrode;
		delete Ef;
		delete Trapping_times;
		delete Probability;
		delete DeTrapping_times;
		delete E_eff;
		delete Etemp;

    // TODO Auto-generated destructor stub
	}

