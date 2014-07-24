/*
 * Interaction.cpp
 *
 *  Created on: Aug 1, 2012
 *      Author: tjanssen
 */

#include "Interaction.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFitResult.h"
#include "TPaveText.h"
#include "TF1.h"
#include "TLegend.h"
#include "Constant.h"

Interaction::Interaction(double Mean[3], double RMS[3], bool Distr_Gaus, double theta, bool doBField, int Nparticles, TRandom3 *R,ElectricField *E, PixelGeometry *Geom) {
        // TODO Auto-generated constructor stub
        t=0;
        rand=R;
        counter = 0;
        Meanxe = Meanye = Meanze = 0;
        RMSxe = RMSye = RMSze = 0;
        Meanxh = Meanyh = Meanzh = 0;
        RMSxh = RMSyh = RMSzh = 0;
        sumxe = sumye = sumze = 0;
        sumxe2 = sumye2 = sumze2 = 0;
        sumxh = sumyh = sumzh = 0;
        sumxh2 = sumyh2 = sumzh2 = 0;
    /** If avarage error becomes larger then this number, we adjust
     *  time step by a downscaling facter. The same if it becomes smaller.
     *  We use an upscaling factor
     */
        TargetError=0.005;
        downscaling_factor=0.8;
        upscaling_factor=1.2;
    /// If simulation runs longer then time_stop we stop the simulation
        Time_stop=0.8*ns;
    /// Initial time step
        timestep = 1e-11*s;
        Coord = new double[3];
        thickness=Geom->GetThickness();
    /// Angles are given in radians in default. Here I convert it to degrees
        alpha = theta*(2*pi)/360;
        Bfield = doBField;

        this->E=E;
        this->Nparticles = Nparticles;
    /// Here I allocate memory on the free space for all the electrons and holes.
        electrons = new Electron * [Nparticles];
        holes = new Hole * [Nparticles];

        NumberOfTraps = 2;
        this->timestep = timestep;

        Hole_Trapping_Times = new double[NumberOfTraps];
        Electron_Trapping_Times = new double[NumberOfTraps];

        Hole_DeTrapping_Times = new double[NumberOfTraps];
        Electron_DeTrapping_Times = new double[NumberOfTraps];
    /** In the following lines the trapping and detrapping times for holes
     *  and electrons are set
     */
        Hole_Trapping_Times[0]=0.01*us;
        Hole_Trapping_Times[1]=0.02*us;

        Electron_Trapping_Times[0]=0.02*us;
        Electron_Trapping_Times[1]=0.01*us;

        Hole_DeTrapping_Times[0]=0.02*us;
        Hole_DeTrapping_Times[1]=0.01*us;

        Electron_DeTrapping_Times[0]=0.01*us;
        Electron_DeTrapping_Times[1]=0.02*us;
/** In this for loop I loop over all charges and set all initial conditions that will affect their
 *  transport. The conditions set here will be used through the entire interaction.
 */
    for(int i=0; i<Nparticles;++i)
    {
        /// Here we create electrons and holes which we will use in interaction.
			anElectron = new Electron();
			aHole = new Hole();
            electrons[i] = anElectron;
            holes[i] = aHole;

        /// Electric field for electrons and holes set.
            electrons[i]->SetField(E);
            holes[i]->SetField(E);
        /// Electrons and holes put in geometry (in other words, inside the sensor)
            electrons[i]->SetGeometry(Geom);
            holes[i]->SetGeometry(Geom);
        /// Randomgenerator set
            electrons[i]->SetRandomGenerator(rand);
            holes[i]->SetRandomGenerator(rand);
        /// Traps and corresponding trapping times set for electrons and holes
            electrons[i]->SetTrap(Electron_Trapping_Times,NumberOfTraps);
            holes[i]->SetTrap(Hole_Trapping_Times,NumberOfTraps);
        /// Detrapping times for each trap set
            electrons[i]->SetDeTrap(Electron_DeTrapping_Times,NumberOfTraps);
            holes[i]->SetDeTrap(Hole_DeTrapping_Times,NumberOfTraps);
        /// Initial position is set
            if(Distr_Gaus==true)
            {
                double *aPosition = this->GenerateGaussianPosition(Mean,RMS);
                electrons[i]->SetPosition(aPosition[0], aPosition[1] ,aPosition[2]);
                holes[i]->SetPosition(aPosition[0], aPosition[1] ,aPosition[2]);

            }
            else
            {
                double *aPosition = this->GenerateUniformPosition(Mean,RMS);
                electrons[i]->SetPosition(aPosition[0], aPosition[1] ,aPosition[2]);
                holes[i]->SetPosition(aPosition[0], aPosition[1] ,aPosition[2]);

            }


    }

/// function Makegraphs is called
    MakeGraphs();

/** This while statement makes sure the simulation will run until simulation is done or a certain time
 *  is reached (this to avoid running the program for way to long).
 */
    while (SimulationNotDone() && t<Time_stop)
    {
            this->DoATimeStep(timestep);

            for(int i=0; i<Nparticles;++i)
            {
            	/// Dotimestep will run the transport (see class charge)
                    electrons[i]->DoTimeStep(timestep, Bfield);
                    holes[i]->DoTimeStep(timestep, Bfield);

            }
        /// Adjust time step if needed
            this->AdjustTimeStep();
    }
}

/** time step adjusted. The factors are set in the constructors and we take the average error from all
 *  particles (so electrons and holes) because we don't want the particles to have different time perspectives
 */
	void Interaction::AdjustTimeStep()
	{
		ComputeAverageError();
        if(Error>TargetError && timestep>0.0001*ns)
                timestep=timestep*downscaling_factor;
        if(Error<TargetError && timestep<0.02*ns)
                timestep=timestep*upscaling_factor;
	}
/// Average error off all particles calculated
	void Interaction::ComputeAverageError()
	{
        int counte=0;
        int counth=0;
        Error=0;

        for(int i=0; i<Nparticles;++i)
        {
                if(electrons[i]->GetIsTrapped()==false && electrons[i]->GetIsCollected()==false)
                {
                        ++counte;
                        Error+=electrons[i]->GetError();
                }
                if(holes[i]->GetIsTrapped()==false && holes[i]->GetIsCollected()==false)
                {
                    ++counth;
                    Error+=holes[i]->GetError();
                }
        }
        if((counte+counth)==0) Error=TargetError;
        else Error/=(counte+counth);
	}

/** All kind of graphs drawn. For the trapping time, detrapping time, RMS, mean etc...
 *  If possible I also fit the graphs.
 */
void Interaction::MakeGraphs()
{
/// Path set to save the plots
    TString plot_path;
    plot_path.Form("/afs/cern.ch/user/t/tjanssen");

/** In the following many many lines I reserve memory for canvases and graphs, so they can be drawn
 *  and saved in the directory defined above.
 */
    TCanvas *canLA = new TCanvas();
    TCanvas *canEF = new TCanvas();
    TGraph *gr_E_eff = new TGraph();
    TGraph *gr_Angle = new TGraph();

    TCanvas *cane1 = new TCanvas();
    TCanvas *cane2 = new TCanvas();
    TCanvas *cane3 = new TCanvas();
    TCanvas *cane4 = new TCanvas();
    TH1I *histo1e = new TH1I("histo1e","Trap Type Distribution Electron",3,0,3);            //for Trap
    TH1I *histo2e = new TH1I("histo2e","Trap Time Distribution Electron",100,0,60);         //TrapTime
    TH1I *histo3e0 = new TH1I("histo3e0","Detrap Time Distribution Electron",100,0,60);      //for DeTrapTime
    TH1I *histo3e1 = new TH1I("histo3e1","Detrap Time Distribution Electron",100,0,60);      //for DeTrapTime
    TH1I *histo4xe = new TH1I("histo4e","Diffusionstep Distribution Electron",100,-3,3);    //for diffusionstep
    TH1I *histo4ye = new TH1I("histo5e","Diffusionstep Distribution Electron",100,-3,3);    //for diffusionstep
    TH1I *histo4ze = new TH1I("histo6e","Diffusionstep Distribution Electron",100,-3,3);    //for diffusionstep
    TCanvas *cane5 = new TCanvas();
    TCanvas *cane6 = new TCanvas();
    TGraph *gr_meanxe = new TGraph();
    TGraph *gr_meanye = new TGraph();
    TGraph *gr_meanze = new TGraph();
    TGraph *gr_rmsxe = new TGraph();
    TGraph *gr_rmsye = new TGraph();
    TGraph *gr_rmsze = new TGraph();

    TCanvas *canh1 = new TCanvas();
    TCanvas *canh2 = new TCanvas();
    TCanvas *canh3 = new TCanvas();
    TCanvas *canh4 = new TCanvas();
    TH1I *histo1h = new TH1I("histo1h","Trap Type Distribution Holes",3,0,3);               //for Trap
    TH1I *histo2h = new TH1I("histo2h","Trap Time Distribution Holes",100,0,60);            //TrapTime
    TH1I *histo3h0 = new TH1I("histo3h0","Detrap Time Distribution Holes",100,0,60);          //for DeTrapTime
    TH1I *histo3h1 = new TH1I("histo3h1","Detrap Time Distribution Holes",100,0,60);          //for DeTrapTime
    TH1I *histo4xh = new TH1I("histo4h","Diffusionstep Distribution Holes",100,-3,3);        //for diffusionstep
    TH1I *histo4yh = new TH1I("histo5h","Diffusionstep Distribution Holes",100,-3,3);        //for diffusionstep
    TH1I *histo4zh = new TH1I("histo6h","Diffusionstep Distribution Holes",100,-3,3);        //for diffusionstep
    TCanvas *canh5 = new TCanvas();
    TCanvas *canh6 = new TCanvas();
    TGraph *gr_meanxh = new TGraph();
    TGraph *gr_meanyh = new TGraph();
    TGraph *gr_meanzh = new TGraph();
    TGraph *gr_rmsxh = new TGraph();
    TGraph *gr_rmsyh = new TGraph();
    TGraph *gr_rmszh = new TGraph();

    TCanvas *candt = new TCanvas();
    TCanvas *canErr = new TCanvas();
    TGraph *gr_dt = new TGraph();
    TGraph *gr_Err = new TGraph();

    TCanvas *canme = new TCanvas();
    TCanvas *canmh = new TCanvas();
    TGraph *gr_Mob430e = new TGraph();
    TGraph *gr_Mob300e = new TGraph();
    TGraph *gr_Mob200e = new TGraph();
    TGraph *gr_Mob430h = new TGraph();
    TGraph *gr_Mob300h = new TGraph();
    TGraph *gr_Mob200h = new TGraph();

    TCanvas *cande = new TCanvas();
    TCanvas *candh = new TCanvas();
    TGraph *gr_Drift430e = new TGraph();
    TGraph *gr_Drift300e = new TGraph();
    TGraph *gr_Drift200e = new TGraph();
    TGraph *gr_Drift430h = new TGraph();
    TGraph *gr_Drift300h = new TGraph();
    TGraph *gr_Drift200h = new TGraph();

    int count_mob=0;
/// In this for loop I draw the graphs of the Mobility and the drift speed vs the electric field
    for(double E=300e3*V/cm;E>0;E-=1000*V/cm)
    {

            electrons[0]->SetTemperature(430);
            electrons[0]->ComputeMobility(E);
            double mobe430=electrons[0]->GetMobility()/(cm*cm/(V*s));

            electrons[0]->SetTemperature(300);
            electrons[0]->ComputeMobility(E);
            double mobe300=electrons[0]->GetMobility()/(cm*cm/(V*s));

            electrons[0]->SetTemperature(200);
            electrons[0]->ComputeMobility(E);
            double mobe200=electrons[0]->GetMobility()/(cm*cm/(V*s));

            holes[0]->SetTemperature(430);
            holes[0]->ComputeMobility(E);
            double mobh430=holes[0]->GetMobility()/(cm*cm/(V*s));

            holes[0]->SetTemperature(300);
            holes[0]->ComputeMobility(E);
            double mobh300=holes[0]->GetMobility()/(cm*cm/(V*s));

            holes[0]->SetTemperature(200);
            holes[0]->ComputeMobility(E);
            double mobh200=holes[0]->GetMobility()/(cm*cm/(V*s));

            gr_Drift430e->SetPoint(count_mob, E*cm/V, E*mobe430*cm/V);
            gr_Drift430e->SetLineColor(8);
            gr_Drift300e->SetPoint(count_mob, E*cm/V, E*mobe300*cm/V);
            gr_Drift300e->SetLineColor(4);
            gr_Drift200e->SetPoint(count_mob, E*cm/V, E*mobe200*cm/V);
            gr_Drift200e->SetLineColor(2);
            gr_Drift430h->SetPoint(count_mob, E*cm/V, E*mobh430*cm/V);
            gr_Drift430h->SetLineColor(8);
            gr_Drift300h->SetPoint(count_mob, E*cm/V, E*mobh300*cm/V);
            gr_Drift300h->SetLineColor(4);
            gr_Drift200h->SetPoint(count_mob, E*cm/V, E*mobh200*cm/V);
            gr_Drift200h->SetLineColor(2);

            gr_Mob430e->SetPoint(count_mob, E*cm/V, mobe430);
            gr_Mob430e->SetLineColor(8);
            gr_Mob300e->SetPoint(count_mob, E*cm/V, mobe300);
            gr_Mob300e->SetLineColor(4);
            gr_Mob200e->SetPoint(count_mob, E*cm/V, mobe200);
            gr_Mob200e->SetLineColor(2);
            gr_Mob430h->SetPoint(count_mob, E*cm/V, mobh430);
            gr_Mob430h->SetLineColor(8);
            gr_Mob300h->SetPoint(count_mob, E*cm/V, mobh300);
            gr_Mob300h->SetLineColor(4);
            gr_Mob200h->SetPoint(count_mob, E*cm/V, mobh200);
            gr_Mob200h->SetLineColor(2);

            ++count_mob;

    }

    cande->cd();
    gr_Drift200e->Draw("AL");
    gr_Drift300e->Draw("same");
    gr_Drift430e->Draw("same");
    gr_Drift430e->SetTitle("430K");
    gr_Drift430e->SetFillColor(0);
    gr_Drift300e->SetTitle("300K");
    gr_Drift300e->SetFillColor(0);
    gr_Drift200e->SetTitle("200K");
    gr_Drift200e->SetFillColor(0);
    gr_Drift200e->GetXaxis()->SetTitle("E (V/cm)");
    gr_Drift200e->GetYaxis()->SetTitle("Drift Speed (cm/s)");
    TLegend *legde = cande->BuildLegend(0.6,0.5,0.89,0.69);
    legde->SetFillColor(0);
    cande->SaveAs(TString::Format("%s/Drift_Speed_Electrons.png",plot_path.Data()));
    delete gr_Drift430e;
    delete gr_Drift300e;
    delete gr_Drift200e;
    delete cande;

    candh->cd();
    gr_Drift200h->Draw("AL");
    gr_Drift300h->Draw("same");
    gr_Drift430h->Draw("same");
    gr_Drift200h->SetTitle("200K");
    gr_Drift200h->SetFillColor(0);
    gr_Drift300h->SetTitle("300K");
    gr_Drift300h->SetFillColor(0);
    gr_Drift430h->SetTitle("430K");
    gr_Drift430h->SetFillColor(0);
    gr_Drift200h->GetXaxis()->SetTitle("E (V/cm)");
    gr_Drift200h->GetYaxis()->SetTitle("Drift Speed (cm/s)");
    TLegend *legdh = candh->BuildLegend(0.6,0.5,0.89,0.69);
    legdh->SetFillColor(0);
    candh->SaveAs(TString::Format("%s/Drift_Speed_Holes.png",plot_path.Data()));
    delete gr_Drift430h;
    delete gr_Drift300h;
    delete gr_Drift200h;
    delete candh;

    canme->cd();
    gr_Mob200e->Draw("AL");
    gr_Mob300e->Draw("same");
    gr_Mob430e->Draw("same");
    gr_Mob430e->SetTitle("430K");
    gr_Mob430e->SetFillColor(0);
    gr_Mob300e->SetTitle("300K");
    gr_Mob300e->SetFillColor(0);
    gr_Mob200e->SetTitle("200K");
    gr_Mob200e->SetFillColor(0);
    gr_Mob200e->GetXaxis()->SetTitle("E (V/cm)");
    gr_Mob200e->GetYaxis()->SetTitle("Mobility (cm^{2}/Vs)");
    TLegend *legme = canme->BuildLegend(0.6,0.5,0.89,0.69);
    legme->SetFillColor(0);
    canme->SaveAs(TString::Format("%s/Mobility_Electrons.png",plot_path.Data()));
    delete gr_Mob430e;
    delete gr_Mob300e;
    delete gr_Mob200e;
    delete canme;




    canmh->cd();
    gr_Mob200h->Draw("AL");
    gr_Mob300h->Draw("same");
    gr_Mob430h->Draw("same");
    gr_Mob200h->SetTitle("200K");
    gr_Mob200h->SetFillColor(0);
    gr_Mob300h->SetTitle("300K");
    gr_Mob300h->SetFillColor(0);
    gr_Mob430h->SetTitle("430K");
    gr_Mob430h->SetFillColor(0);
    gr_Mob200h->GetXaxis()->SetTitle("E (V/cm)");
    gr_Mob200h->GetYaxis()->SetTitle("Mobility (cm^{2}/Vs)");
    TLegend *legmh = canmh->BuildLegend(0.6,0.5,0.89,0.69);
    legmh->SetFillColor(0);
    canmh->SaveAs(TString::Format("%s/Mobility_Holes.png",plot_path.Data()));
    delete gr_Mob430h;
    delete gr_Mob300h;
    delete gr_Mob200h;
    delete canmh;

/** If magnetic field is turned on the drift of electrons and holes is tilted by an angle called the Lorentz
 *  angle. This magnetic field results in a Lorentz force that will influence the mobility
 *  We therefore write the Lorentz force in terms of an effective electric field
 *  In this function a plot of the Lorentz angle vs Ez is made and a plot of
 *  Ex vs Ez. These plots are made with ROOT.
 */
    for(double E=100000*V/cm;E>0;E-=100*V/cm)
    {
    electrons[0]->GetElectricEffectiveField(E);
    gr_E_eff->SetPoint(E/(100*V/cm), E*cm/V, electrons[0]->GetEffX()*cm/V);
    gr_Angle->SetPoint(E/(100*V/cm), E*cm/V, electrons[0]->GetLA()*360/(2*pi));
    }

    canEF->cd();
    gr_E_eff->Draw("AL");
    gStyle->SetOptStat(0);
    gr_E_eff->SetTitle("Electric field for Electrons");
    gr_E_eff->GetXaxis()->SetTitle("Ez (V/cm)");
    gr_E_eff->GetYaxis()->SetTitle("Ex (V/cm)");
    canEF->SaveAs(TString::Format("%s/Electric_Field.png",plot_path.Data()));

    canLA->cd();
    gr_Angle->Draw("AL");
    gStyle->SetOptStat(0);
    gr_Angle->SetTitle("Lorentz angle for Electrons");
    gr_Angle->GetXaxis()->SetTitle("Ez (V/cm)");
    gr_Angle->GetYaxis()->SetTitle("Angle (degree)");

    canLA->SaveAs(TString::Format("%s/Angle.png",plot_path.Data()));
    delete gr_E_eff;
    delete canLA;
    delete gr_Angle;
    delete canEF;

/** In this for loop I make histograms of the diffusionstep, trapping time,
 *  detrapping time and trap type.
 */
    for(int i=0; i<Nparticles;++i)
    {
            electrons[i]->SetTimeBeforeTrapped();
            holes[i]->SetTimeBeforeTrapped();

            electrons[i]->SetTimeBeforeDeTrapped();
            holes[i]->SetTimeBeforeDeTrapped();

            electrons[i]->DoDiffusionStep(timestep);
            holes[i]->DoDiffusionStep(timestep);

            histo1e->Fill(electrons[i]->GetTrap());
            histo2e->Fill(electrons[i]->GetTimeToTrap());

            for(int j=0;j<NumberOfTraps;++j)
            {
            	if(electrons[i]->GetTrap()==1){
            		histo3e0->Fill(electrons[i]->GetTimeToDetrap());
            		histo3e0->SetLineColor(8);
            	}
            	if(electrons[i]->GetTrap()==2){
            		histo3e1->Fill(electrons[i]->GetTimeToDetrap());
            		histo3e1->SetLineColor(6);
            	}

            	if(holes[i]->GetTrap()==1){
            		histo3h0->Fill(holes[i]->GetTimeToDetrap());
            		histo3h0->SetLineColor(8);
            	}
            	if(holes[i]->GetTrap()==2){
            		histo3h1->Fill(holes[i]->GetTimeToDetrap());
            		histo3h1->SetLineColor(6);
            	}
            }

            histo4xe->Fill(electrons[i]->GetDiffusionStepX());
            histo4xe->SetLineColor(8);
            histo4ye->Fill(electrons[i]->GetDiffusionStepY());
            histo4ye->SetLineColor(6);
            histo4ze->Fill(electrons[i]->GetDiffusionStepZ());
            histo4ze->SetLineColor(2);



            histo1h->Fill(holes[i]->GetTrap());
            histo2h->Fill(holes[i]->GetTimeToTrap());
            histo4xh->Fill(holes[i]->GetDiffusionStepX());
            histo4xh->SetLineColor(8);
            histo4yh->Fill(holes[i]->GetDiffusionStepY());
            histo4yh->SetLineColor(6);
            histo4zh->Fill(holes[i]->GetDiffusionStepZ());
            histo4zh->SetLineColor(2);
    }

/** In this while loop I make the graphs of the mean and RMS of the carrier distribution
 *  in the sensor during the simulation. For the RMS I have also the option to make a fit.
 *  These fits will only look well if the trapping time is very large.
 */
    while (SimulationNotDone() && t<Time_stop)
    {
            this->DoATimeStep(timestep);
            sumxe=0;
            sumxh=0;
            sumye=0;
            sumyh=0;
            sumze=0;
            sumzh=0;
            sumxe2=0;
            sumxh2=0;
            sumye2=0;
            sumyh2=0;
            sumze2=0;
            sumzh2=0;

            for(int i=0; i<Nparticles;++i)
            {
                    sumxe+=electrons[i]->GetX();
                    sumxe2+=pow(electrons[i]->GetX(),2);
                    sumye+=electrons[i]->GetY();
                    sumye2+=pow(electrons[i]->GetY(),2);
                    sumze+=electrons[i]->GetZ();
                    sumze2+=pow(electrons[i]->GetZ(),2);

                    sumxh+=holes[i]->GetX();
                    sumxh2+=pow(holes[i]->GetX(),2);
                    sumyh+=holes[i]->GetY();
                    sumyh2+=pow(holes[i]->GetY(),2);
                    sumzh+=holes[i]->GetZ();
                    sumzh2+=pow(holes[i]->GetZ(),2);

                    electrons[i]->DoTimeStep(timestep, Bfield);
                    holes[i]->DoTimeStep(timestep, Bfield);

            }
            this->AdjustTimeStep();

            Meanxe = sumxe/Nparticles;
            Meanye = sumye/Nparticles;
            Meanze = sumze/Nparticles;
            RMSxe = sqrt(sumxe2/Nparticles);
            RMSye = sqrt(sumye2/Nparticles);
            RMSze = sqrt(sumze2/Nparticles);

            gr_meanxe->SetPoint(counter, t, Meanxe);
            gr_meanxe->SetLineColor(8);
            gr_meanye->SetPoint(counter, t, Meanye);
            gr_meanye->SetLineColor(4);
            gr_meanze->SetPoint(counter, t, Meanze);
            gr_meanze->SetLineColor(2);
            gr_rmsxe->SetPoint(counter, t, RMSxe);
            gr_rmsxe->SetLineColor(8);
            gr_rmsye->SetPoint(counter, t, RMSye);
            gr_rmsye->SetLineColor(4);
            gr_rmsze->SetPoint(counter, t, RMSze);
            gr_rmsze->SetLineColor(2);
            gr_dt->SetPoint(counter, t, timestep);
            gr_Err->SetPoint(counter, t, Error);

            Meanxh = sumxh/Nparticles;
            Meanyh = sumyh/Nparticles;
            Meanzh = sumzh/Nparticles;
            RMSxh = sqrt(sumxh2/Nparticles);
            RMSyh = sqrt(sumyh2/Nparticles);
            RMSzh = sqrt(sumzh2/Nparticles);

            gr_meanxh->SetPoint(counter, t, Meanxh);
            gr_meanxh->SetLineColor(8);
            gr_meanyh->SetPoint(counter, t, Meanyh);
            gr_meanyh->SetLineColor(4);
            gr_meanzh->SetPoint(counter, t, Meanzh);
            gr_meanzh->SetLineColor(2);
            gr_rmsxh->SetPoint(counter, t, RMSxh);
            gr_rmsxh->SetLineColor(8);
            gr_rmsyh->SetPoint(counter, t, RMSyh);
            gr_rmsyh->SetLineColor(4);
            gr_rmszh->SetPoint(counter, t, RMSzh);
            gr_rmszh->SetLineColor(2);
            ++counter;

    }

    //Timestep vs time is drawn to check how much the timestep changes over time.
    candt->cd();
    gr_dt->Draw("AL");
    gStyle->SetOptStat(0);
    gr_dt->SetTitle("Timestep vs Time");
    gr_dt->GetXaxis()->SetTitle("Time (ns)");
    gr_dt->GetYaxis()->SetTitle("Timestep (ns)");
    candt->SaveAs(TString::Format("%s/TimeStep.png",plot_path.Data()));
    delete gr_dt;
    delete candt;

    //Graph of the error vs time made to check if adjust time is doing its work properly
    canErr->cd();
    gr_Err->Draw("AL");
    gStyle->SetOptStat(0);
    gr_Err->SetTitle("Error vs Time");
    gr_Err->GetXaxis()->SetTitle("Time (ns)");
    gr_Err->GetYaxis()->SetTitle("Average Error");
    canErr->SaveAs(TString::Format("%s/Error.png",plot_path.Data()));
    delete gr_Err;
    delete canErr;

    //Histogram for trap type. This one doesnt need a fit. We just want to check how many times a charge
    //fals in trap l
    cane1->cd();
    histo1e->Draw();
    gStyle->SetOptStat(0);
    histo1e->GetXaxis()->SetTitle("Trap Type");
    histo1e->GetYaxis()->SetTitle("A.U.");
    cane1->SaveAs(TString::Format("%s/TrapDistrElectrons.png",plot_path.Data()));
    delete histo1e;
    delete cane1;

    //Histogram of trap time. This should be an exponential. So that is what we check here.
    TF1 *f1e= new TF1("histo2efit","[0]*TMath::Exp(-[1]*x)",0,100);
    TFitResultPtr r1e =  histo2e->Fit(f1e,"S","",0,60);
    cane2->cd();
    histo2e->Draw();
    TPaveText *par1e = new TPaveText(0.5,0.7,0.89,0.89,"NDC");
    par1e->SetFillColor(0);
    par1e->SetTextAlign(12);
    par1e->SetTextSize(0.04);
    par1e->AddText(TString::Format("#tau = %f",1.0/r1e->Parameter(1)));
    par1e->AddText(TString::Format("Amplitude = %f",r1e->Parameter(0)));
    par1e->Draw();
    histo2e->GetXaxis()->SetTitle("Trap Time (ns)");
    histo2e->GetYaxis()->SetTitle("A.U.");
    cane2->SaveAs(TString::Format("%s/TrapTimeElectrons.png",plot_path.Data()));
    delete histo2e;
    delete f1e;
    delete par1e;
    delete cane2;

    //Histogram of the detrap time for trap 0. This also has to be an exponential
    TF1 *f2e0= new TF1("histo3e0fit","[0]*TMath::Exp(-[1]*x)",0,100);
    f2e0->SetLineColor(1);
    TF1 *f2e1= new TF1("histo3e1fit","[0]*TMath::Exp(-[1]*x)",0,100);
    f2e1->SetLineColor(3);
    TFitResultPtr r2e0 =  histo3e0->Fit(f2e0,"S","",0,60);
    TFitResultPtr r2e1 =  histo3e1->Fit(f2e1,"S","",0,60);
    cane3->cd();
    histo3e0->Draw();
    histo3e1->Draw("same");
    TPaveText *par2e0 = new TPaveText(0.5,0.7,0.89,0.89,"NDC");
    par2e0->SetFillColor(0);
    par2e0->SetTextAlign(12);
    par2e0->SetTextSize(0.04);
    par2e0->AddText(TString::Format("#tau = %f",1.0/r2e0->Parameter(1)));
    par2e0->AddText(TString::Format("Amplitude = %f",r2e0->Parameter(0)));
    par2e0->Draw("");
    TPaveText *par2e1 = new TPaveText(0.5,0.5,0.89,0.69,"NDC");
    par2e1->SetFillColor(0);
    par2e1->SetTextAlign(12);
    par2e1->SetTextSize(0.04);
    par2e1->AddText(TString::Format("#tau = %f",1.0/r2e1->Parameter(1)));
    par2e1->AddText(TString::Format("Amplitude = %f",r2e1->Parameter(0)));
    par2e1->Draw("");
    histo3e0->SetTitle("Trap 1");
    histo3e0->SetFillColor(0);
    histo3e1->SetTitle("Trap 2");
    histo3e1->SetFillColor(0);
    histo3e0->GetXaxis()->SetTitle("Detrap Time (ns)");
    histo3e0->GetYaxis()->SetTitle("A.U.");
    TLegend *lege3 = cane3->BuildLegend(0.2,0.7,0.39,0.89);
    lege3->SetFillColor(0);
    cane3->SaveAs(TString::Format("%s/DeTrapTimeElectrons.png",plot_path.Data()));
    delete histo3e1;
    delete histo3e0;
    delete f2e0;
    delete f2e1;
    delete par2e0;
    delete par2e1;
    delete cane3;

    //Histrogram of the diffusionstep (in this case in the Z direction). This should be a gaussian.
    TF1 *f3xe= new TF1("histo4xefit","gaus",0,2);
    f3xe->SetLineColor(1);
    TF1 *f3ye= new TF1("histo4yefit","gaus",0,2);
    f3ye->SetLineColor(3);
    TF1 *f3ze= new TF1("histo4zefit","gaus",0,2);
    f3ze->SetLineColor(5);
    TFitResultPtr r3xe =  histo4xe->Fit(f3xe,"S","",-3,3);
    TFitResultPtr r3ye =  histo4ye->Fit(f3ye,"S","",-3,3);
    TFitResultPtr r3ze =  histo4ze->Fit(f3ze,"S","",-3,3);
    cane4->cd();
    histo4xe->Draw();
    histo4ye->Draw("same");
    histo4ze->Draw("same");
    TPaveText *par3xe = new TPaveText(0.62,0.2,0.89,0.39,"NDC");
    par3xe->SetFillColor(0);
    par3xe->SetTextAlign(12);
    par3xe->SetTextSize(0.04);
    par3xe->AddText(TString::Format("x:Const = %f",r3xe->Parameter(0)));
    par3xe->AddText(TString::Format("x:Mean = %f",r3xe->Parameter(1)));
    par3xe->AddText(TString::Format("x:RMS = %f",r3xe->Parameter(2)));
    par3xe->Draw();
    TPaveText *par3ye = new TPaveText(0.62,0.4,0.89,0.59,"NDC");
    par3ye->SetFillColor(0);
    par3ye->SetTextAlign(12);
    par3ye->SetTextSize(0.04);
    par3ye->AddText(TString::Format("y:Const = %f",r3ye->Parameter(0)));
    par3ye->AddText(TString::Format("y:Mean = %f",r3ye->Parameter(1)));
    par3ye->AddText(TString::Format("z:RMS = %f",r3ye->Parameter(2)));
    par3ye->Draw();
    TPaveText *par3ze = new TPaveText(0.62,0.6,0.89,0.79,"NDC");
    par3ze->SetFillColor(0);
    par3ze->SetTextAlign(12);
    par3ze->SetTextSize(0.04);
    par3ze->AddText(TString::Format("z:Const = %f",r3ze->Parameter(0)));
    par3ze->AddText(TString::Format("z:Mean = %f",r3ze->Parameter(1)));
    par3ze->AddText(TString::Format("z:RMS = %f",r3ze->Parameter(2)));
    par3ze->Draw();
    histo4xe->SetTitle("x");
    histo4xe->SetFillColor(0);
    histo4ye->SetTitle("y");
    histo4ye->SetFillColor(0);
    histo4ze->SetTitle("z");
    histo4ze->SetFillColor(0);
    histo4xe->GetXaxis()->SetTitle("Diffusionstep (#mum)");
    histo4xe->GetYaxis()->SetTitle("A.U.");
    TLegend *lege4 = cane4->BuildLegend(0.2,0.7,0.39,0.89);
    lege4->SetFillColor(0);
    cane4->SaveAs(TString::Format("%s/DiffusionStepElectrons.png",plot_path.Data()));
    delete histo4xe;
    delete histo4ye;
    delete histo4ze;
    delete f3xe;
    delete f3ye;
    delete f3ze;
    delete par3xe;
    delete par3ye;
    delete par3ze;
    delete cane4;

    cane5->cd();
    gr_meanze->Draw("AL");
    gr_meanye->Draw("same");
    gr_meanxe->Draw("same");
    gStyle->SetOptStat(0);
    gr_meanxe->SetTitle("x");
    gr_meanxe->SetFillColor(0);
    gr_meanye->SetTitle("y");
    gr_meanye->SetFillColor(0);
    gr_meanze->SetTitle("z");
    gr_meanze->SetFillColor(0);
    gr_meanze->GetYaxis()->SetRangeUser(-10,25);
    gr_meanze->GetXaxis()->SetTitle("Time (ns)");
    gr_meanze->GetYaxis()->SetTitle("Mean (#mum)");
    TLegend *lege5 = cane5->BuildLegend(0.7,0.6,0.89,0.79);
    lege5->SetFillColor(0);
    cane5->SaveAs(TString::Format("%s/MeanElectrons.png",plot_path.Data()));
    delete gr_meanxe;
    delete gr_meanye;
    delete gr_meanze;
    delete cane5;

//    TF1 *f4xe= new TF1("gr_rmsefit","TMath::Sqrt([0]*x)+[1]",0,10);
//    f4xe->SetLineColor(1);
//    TF1 *f4ye= new TF1("gr_rmsefit","TMath::Sqrt([0]*x)+[1]",0,10);
//    f4ye->SetLineColor(3);
//    TF1 *f4ze= new TF1("gr_rmsefit","TMath::Sqrt([0]*x)+[1]",0,10);
//    f4ze->SetLineColor(5);
//    f4xe->SetParameter(0,36*cm*cm/s);
//    f4ye->SetParameter(0,36*cm*cm/s);
//    f4ze->SetParameter(0,36*cm*cm/s);
//    TFitResultPtr r4xe =  gr_rmsxe->Fit(f4xe,"S","",0,0.4);
//    TFitResultPtr r4ye =  gr_rmsye->Fit(f4ye,"S","",0,0.4);
//    TFitResultPtr r4ze =  gr_rmsze->Fit(f4ze,"S","",0,0.4);
    cane6->cd();
    gr_rmsze->Draw("AL");
    gr_rmsye->Draw("same");
    gr_rmsxe->Draw("same");
//    TPaveText *par4xe = new TPaveText(0.6,0.2,0.89,0.29,"NDC");
//    par4xe->SetFillColor(0);
//    par4xe->SetTextAlign(12);
//    par4xe->SetTextSize(0.04);
//    par4xe->AddText(TString::Format("x:D = %f",0.5*r4xe->Parameter(0)));
//    par4xe->Draw();
//    TPaveText *par4ye = new TPaveText(0.6,0.3,0.89,0.39,"NDC");
//    par4ye->SetFillColor(0);
//    par4ye->SetTextAlign(12);
//    par4ye->SetTextSize(0.04);
//    par4ye->AddText(TString::Format("y:D = %f",0.5*r4ye->Parameter(0)));
//    par4ye->Draw();
//    TPaveText *par4ze = new TPaveText(0.6,0.4,0.89,0.49,"NDC");
//    par4ze->SetFillColor(0);
//    par4ze->SetTextAlign(12);
//    par4ze->SetTextSize(0.04);
//    par4ze->AddText(TString::Format("z:D = %f",0.5*r4ze->Parameter(0)));
//    par4ze->Draw();
    gr_rmsxe->SetTitle("x");
    gr_rmsxe->SetFillColor(0);
    gr_rmsye->SetTitle("y");
    gr_rmsye->SetFillColor(0);
    gr_rmsze->SetTitle("z");
    gr_rmsze->SetFillColor(0);
    gr_rmsze->GetYaxis()->SetRangeUser(0,25);
    gr_rmsze->GetXaxis()->SetTitle("Time (ns)");
    gr_rmsze->GetYaxis()->SetTitle("RMS (#mum)");
    TLegend *lege6 = cane6->BuildLegend(0.7,0.5,0.89,0.69);
    lege6->SetFillColor(0);
    cane6->SaveAs(TString::Format("%s/RMSElectrons.png",plot_path.Data()));
    delete gr_rmsxe;
    delete gr_rmsye;
    delete gr_rmsze;
//    delete f4xe;
//    delete f4ye;
//    delete f4ze;
//    delete par4xe;
//    delete par4ye;
//    delete par4ze;
    delete cane6;

    canh1->cd();
    histo1h->Draw();
    gStyle->SetOptStat(0);
    histo1h->GetXaxis()->SetTitle("Trap Type");
    histo1h->GetYaxis()->SetTitle("A.U.");
    canh1->SaveAs(TString::Format("%s/TrapDistrHoles.png",plot_path.Data()));
    delete histo1h;
    delete canh1;

    TF1 *f1h= new TF1("histo2hfit","[0]*TMath::Exp(-[1]*x)",0,100);
    TFitResultPtr r1h =  histo2h->Fit(f1h,"S","",0,60);
    canh2->cd();
    histo2h->Draw();
    TPaveText *par1h = new TPaveText(0.5,0.7,0.89,0.89,"NDC");
    par1h->SetFillColor(0);
    par1h->SetTextAlign(12);
    par1h->SetTextSize(0.04);
    par1h->AddText(TString::Format("#tau = %f",1.0/r1h->Parameter(1)));
    par1h->AddText(TString::Format("Amplitude = %f",r1h->Parameter(0)));
    par1h->Draw();
    histo2h->GetXaxis()->SetTitle("Trap Time (ns)");
    histo2h->GetYaxis()->SetTitle("A.U.");
    canh2->SaveAs(TString::Format("%s/TrapTimeHoles.png",plot_path.Data()));
    delete histo2h;
    delete f1h;
    delete par1h;
    delete canh2;

    TF1 *f2h0= new TF1("histo3h0fit","[0]*TMath::Exp(-[1]*x)",0,100);
    f2h0->SetLineColor(1);
    TF1 *f2h1= new TF1("histo3h1fit","[0]*TMath::Exp(-[1]*x)",0,100);
    f2h1->SetLineColor(3);
    f2h0->SetParameter(1,0.05);
    f2h1->SetParameter(1,0.01);
    TFitResultPtr r2h0 =  histo3h0->Fit(f2h0,"S","",0,60);
    TFitResultPtr r2h1 =  histo3h1->Fit(f2h1,"S","",0,60);
    canh3->cd();
    histo3h0->Draw();
    histo3h1->Draw("same");
    TPaveText *par2h0 = new TPaveText(0.5,0.7,0.89,0.89,"NDC");
    par2h0->SetFillColor(0);
    par2h0->SetTextAlign(12);
    par2h0->SetTextSize(0.04);
    par2h0->AddText(TString::Format("#tau = %f",1.0/r2h0->Parameter(1)));
    par2h0->AddText(TString::Format("Amplitude = %f",r2h0->Parameter(0)));
    par2h0->Draw();
    TPaveText *par2h1 = new TPaveText(0.5,0.5,0.89,0.69,"NDC");
    par2h1->SetFillColor(0);
    par2h1->SetTextAlign(12);
    par2h1->SetTextSize(0.04);
    par2h1->AddText(TString::Format("#tau = %f",1.0/r2h1->Parameter(1)));
    par2h1->AddText(TString::Format("Amplitude = %f",r2h1->Parameter(0)));
    par2h1->Draw();
    histo3h0->SetTitle("Trap 1");
    histo3h0->SetFillColor(0);
    histo3h1->SetTitle("Trap 2");
    histo3h1->SetFillColor(0);
    histo3h0->GetXaxis()->SetTitle("Detrap Time (ns)");
    histo3h0->GetYaxis()->SetTitle("A.U.");
    TLegend *legh3 = canh3->BuildLegend(0.2,0.7,0.39,0.89);
    legh3->SetFillColor(0);
    canh3->SaveAs(TString::Format("%s/DeTrapTimeHoles.png",plot_path.Data()));
    delete histo3h1;
    delete histo3h0;
    delete f2h0;
    delete f2h1;
    delete par2h0;
    delete par2h1;
    delete canh3;

    TF1 *f3xh= new TF1("histo4xhfit","gaus",0,100);
    f3xh->SetLineColor(1);
    TF1 *f3yh= new TF1("histo4yhfit","gaus",0,100);
    f3yh->SetLineColor(3);
    TF1 *f3zh= new TF1("histo4zhfit","gaus",0,100);
    f3zh->SetLineColor(5);
    TFitResultPtr r3xh =  histo4xh->Fit(f3xh,"S","",-3,3);
    TFitResultPtr r3yh =  histo4yh->Fit(f3yh,"S","",-3,3);
    TFitResultPtr r3zh =  histo4zh->Fit(f3zh,"S","",-3,3);
    canh4->cd();
    histo4xh->Draw();
    histo4yh->Draw("same");
    histo4zh->Draw("same");
    TPaveText *par3xh = new TPaveText(0.62,0.2,0.89,0.39,"NDC");
    par3xh->SetFillColor(0);
    par3xh->SetTextAlign(12);
    par3xh->SetTextSize(0.04);
    par3xh->AddText(TString::Format("x:Const = %f",r3xh->Parameter(0)));
    par3xh->AddText(TString::Format("x:Mean = %f",r3xh->Parameter(1)));
    par3xh->AddText(TString::Format("x:RMS = %f",r3xh->Parameter(2)));
    par3xh->Draw();
    TPaveText *par3yh = new TPaveText(0.62,0.4,0.89,0.59,"NDC");
    par3yh->SetFillColor(0);
    par3yh->SetTextAlign(12);
    par3yh->SetTextSize(0.04);
    par3yh->AddText(TString::Format("y:Const = %f",r3xh->Parameter(0)));
    par3yh->AddText(TString::Format("y:Mean = %f",r3xh->Parameter(1)));
    par3yh->AddText(TString::Format("y:RMS = %f",r3xh->Parameter(2)));
    par3yh->Draw();
    TPaveText *par3zh = new TPaveText(0.62,0.6,0.89,0.79,"NDC");
    par3zh->SetFillColor(0);
    par3zh->SetTextAlign(12);
    par3zh->SetTextSize(0.04);
    par3zh->AddText(TString::Format("z:Const = %f",r3xh->Parameter(0)));
    par3zh->AddText(TString::Format("z:Mean = %f",r3xh->Parameter(1)));
    par3zh->AddText(TString::Format("z:RMS = %f",r3xh->Parameter(2)));
    par3zh->Draw();
    histo4xh->SetTitle("x");
    histo4xh->SetFillColor(0);
    histo4yh->SetTitle("y");
    histo4yh->SetFillColor(0);
    histo4zh->SetTitle("z");
    histo4zh->SetFillColor(0);
    histo4xh->GetXaxis()->SetTitle("Diffusionstep (#mum)");
    histo4xh->GetYaxis()->SetTitle("A.U.");
    TLegend *legh4 = canh4->BuildLegend(0.2,0.7,0.39,0.89);
    legh4->SetFillColor(0);
    canh4->SaveAs(TString::Format("%s/DiffusionStepHoles.png",plot_path.Data()));
    delete histo4xh;
    delete histo4yh;
    delete histo4zh;
    delete f3xh;
    delete f3yh;
    delete f3zh;
    delete par3xh;
    delete par3yh;
    delete par3zh;
    delete canh4;

    canh5->cd();
    gr_meanzh->Draw("AL");
    gr_meanyh->Draw("same");
    gr_meanxh->Draw("same");
    gStyle->SetOptStat(0);
    gr_meanzh->GetYaxis()->SetRangeUser(-25,5);
    gr_meanzh->GetXaxis()->SetTitle("Time (ns)");
    gr_meanzh->GetYaxis()->SetTitle("Mean (#mum)");
    gr_meanxh->SetTitle("x");
    gr_meanxh->SetFillColor(0);
    gr_meanyh->SetTitle("y");
    gr_meanyh->SetFillColor(0);
    gr_meanzh->SetTitle("z");
    gr_meanzh->SetFillColor(0);
    TLegend *legh5 = canh5->BuildLegend(0.2,0.2,0.39,0.39);
    legh5->SetFillColor(0);
    canh5->SaveAs(TString::Format("%s/MeanHoles.png",plot_path.Data()));
    delete gr_meanxh;
    delete gr_meanyh;
    delete gr_meanzh;
    delete canh5;

//    TF1 *f4xh= new TF1("gr_rmshfit","TMath::Sqrt([0]*x)+[1]",0,0.8);
//    f4xh->SetLineColor(1);
//    TF1 *f4yh= new TF1("gr_rmshfit","TMath::Sqrt([0]*x)+[1]",0,0.8);
//    f4yh->SetLineColor(3);
//    TF1 *f4zh= new TF1("gr_rmshfit","TMath::Sqrt([0]*x)+[1]",0,0.8);
//    f4zh->SetLineColor(5);
//    f4xh->SetParameter(0,12*cm*cm/s);
//    f4yh->SetParameter(0,12*cm*cm/s);
//    f4zh->SetParameter(0,12*cm*cm/s);
//    TFitResultPtr r4xh =  gr_rmsxh->Fit(f4xh,"S","",0,9);
//    TFitResultPtr r4yh =  gr_rmsyh->Fit(f4yh,"S","",0,9);
//    TFitResultPtr r4zh =  gr_rmszh->Fit(f4zh,"S","",0,9);
    canh6->cd();
    gr_rmszh->Draw("AL");
    gr_rmsyh->Draw("same");
    gr_rmsxh->Draw("same");
//    TPaveText *par4xh = new TPaveText(0.6,0.2,0.89,0.29,"NDC");
//    par4xh->SetFillColor(0);
//    par4xh->SetTextAlign(12);
//    par4xh->SetTextSize(0.04);
//    par4xh->AddText(TString::Format("x:D = %f",0.5*r4xh->Parameter(0)));
//    par4xh->Draw();
//    TPaveText *par4yh = new TPaveText(0.6,0.3,0.89,0.39,"NDC");
//    par4yh->SetFillColor(0);
//    par4yh->SetTextAlign(12);
//    par4yh->SetTextSize(0.04);
//    par4yh->AddText(TString::Format("y:D = %f",0.5*r4yh->Parameter(0)));
//    par4yh->Draw();
//    TPaveText *par4zh = new TPaveText(0.6,0.4,0.89,0.49,"NDC");
//    par4zh->SetFillColor(0);
//    par4zh->SetTextAlign(12);
//    par4zh->SetTextSize(0.04);
//    par4zh->AddText(TString::Format("z:D = %f",0.5*r4zh->Parameter(0)));
//    par4zh->Draw();
    gr_rmsxh->SetTitle("x");
    gr_rmsxh->SetFillColor(0);
    gr_rmsyh->SetTitle("y");
    gr_rmsyh->SetFillColor(0);
    gr_rmszh->SetTitle("z");
    gr_rmszh->SetFillColor(0);
    gr_rmszh->GetYaxis()->SetRangeUser(0,25);
    gr_rmszh->GetXaxis()->SetTitle("Time (ns)");
    gr_rmszh->GetYaxis()->SetTitle("RMS (#mum)");
    TLegend *legh6 = canh6->BuildLegend(0.7,0.5,0.89,0.69);
    legh6->SetFillColor(0);
    canh6->SaveAs(TString::Format("%s/RMSHoles.png",plot_path.Data()));
    delete gr_rmsxh;
    delete gr_rmsyh;
    delete gr_rmszh;
//    delete f4xh;
//    delete f4yh;
//    delete f4zh;
//    delete par4xh;
//    delete par4yh;
//    delete par4zh;
    delete canh6;

}

/** Checking if simulation is done or not done. We want to run the simulation till all particles are either
 *  collected or trapped for more than a certain time.
 */
	bool Interaction::SimulationNotDone()
	{

        bool notDone = false;
        for(int i=0; i<Nparticles;++i)
        {
                if(!electrons[i]->GetIsCollected()) {

                        if (electrons[i]->GetIsTrapped()){
                                if (electrons[i]->GetTimeToDetrap() < 0.025*us) notDone=true;
                                }
                        else notDone=true;
                };


                if(!holes[i]->GetIsCollected()) {

                        if (holes[i]->GetIsTrapped()){
                                if (holes[i]->GetTimeToDetrap() < 0.025*us) notDone=true;
                                }
                        else notDone=true;
                };
        }
        //cout << "[Collection] " <<  notDone << endl;
        return notDone;
	}

/** Initial position of charge set. This position is uniform in z and Gaussian in x and y.
 *  I also rotate over the x axis under an angle alpha (theta in the constructor).
 */
	double* Interaction::GenerateUniformPosition(double mean[3], double RMS[3])
	{
        double u1, u2, v1;
        u1 = rand->Uniform();
        u2 = rand->Uniform();
        v1 = rand->Uniform();

        u = RMS[0]*sin(2*pi*u1)*sqrt(-2*log(u2))+mean[0];
        v = RMS[1]*cos(2*pi*u2)*sqrt(-2*log(u1))+mean[1];
        w = v1*thickness-thickness/2;

        u_ = u;
        v_ = cos(alpha)*v-sin(alpha)*w;
        w_ = sin(alpha)*v+cos(alpha)*w;

        Coord[0] = u_;
        Coord[1] = v_;
        Coord[2] = w_;

        return Coord;

	}
/// Initial position of charge set using a gaussian distribution.
double* Interaction::GenerateGaussianPosition(double mean[3], double RMS[3])
{
                double u1, u2, v1, v2;
                u1 = rand->Uniform();
                u2 = rand->Uniform();
                v1 = rand->Uniform();
                v2 = rand->Uniform();

                u = RMS[0]*sin(2*pi*u1)*sqrt(-2*log(u2))+mean[0];
                v = RMS[1]*cos(2*pi*u2)*sqrt(-2*log(u1))+mean[1];
                w = RMS[2]*sin(2*pi*v1)*sqrt(-2*log(v2))+mean[2];

                Coord[0]=u;
                Coord[1]=v;
                Coord[2]=w;

        return Coord;
}

/// Timestep defined in interaction.
	void Interaction::DoATimeStep(double dt)
	{
        t = t + dt;
	}

/// Delete all memory allocated in the free space.
	Interaction::~Interaction()
	{
        delete Hole_Trapping_Times;
        delete Electron_Trapping_Times;
        delete Hole_DeTrapping_Times;
        delete Electron_DeTrapping_Times;
        delete Coord;
        delete electrons;
        delete holes;
        delete aHole;
        delete anElectron;
        // TODO Auto-generated destructor stub
	}

