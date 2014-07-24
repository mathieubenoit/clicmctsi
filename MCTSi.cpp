//============================================================================
// Name        : MCTSi.cpp
// Author      : Tim Janssen
// Version     :
// Copyright   : copyleft
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "Constant.h"
#include "Electron.h"
#include "Hole.h"
#include "TH1I.h"
#include "TApplication.h"
#include "TRandom3.h"
#include "Interaction.h"
#include "PixelGeometry.h"
#include "AllPixDigitAnimation.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFitResult.h"
#include "TPaveText.h"
#include "TF1.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TLegend.h"
using namespace std;

/** \mainpage
 *  This Monte Carlo program simulates the charge collection behavior of semiconductor silicon tracking detectors.
 */
int main(int argc, char* argv[]) {


        double SensorThickness=50*um;
    /** This appliction makes sure I can use root in C++. This is the only application I wont delete, because
     *  otherwise the graphs will not be drawn properly. So one has to shutdown the program himself after
     *  all the graphs are drawn.
     */
        TApplication* rootapp = new TApplication("example",&argc, argv);
        gROOT->SetStyle("Pub");
    /** random number generator defined (there are 3 random number generators in root, this one is the best
     *  in this situation, because of the long lifetime and its randomness)
     */
        TRandom3 *rand = new TRandom3();

    /// ElectricField set
        ElectricField *E = new ElectricField(40,10,SensorThickness,true,"/afs/cern.ch/user/t/tjanssen/Ef.root");
    /// If uncomment this line, it will draw the strength of the electric field, this will however greatly reduce the performance
        //E->Ey->Draw("colz");
    /// Geometry set
        PixelGeometry *Geom = new PixelGeometry(SensorThickness,256,256,25*um,25*um);

    /// Set random numbers between 0 and 1
        rand->SetSeed(1);

        cout << "Monte-Carlo Charge Transport in Silicon" << endl; // prints Monte-Carlo Charge Transport in Silicon

    /// Here we define the RMS and the mean for the initial distribution
        double RMS[3]={0.1*um,0.1*um,1*um};
        double Mean[3]={0*um,0*um,0};

        int ntracks =100;
    /** Mean and RMS are set to generate the initial parameters for
     *  the Gaussian distribution to give all particles
     *  the initial position. The boolean is to generate a Gaussian
     *  position if true or make z uniform while rotating
     *  over the x axis under an angle theta if false.
     *  After the choice of the distribution one can set the angle which
     *  will only be used if the boolean is false. The second boolean is
     *  to switch the magnetic field (true) on or off (false)
     */
        Interaction I(Mean,RMS, false, 20, true, ntracks,rand,E,Geom);

    /** The following 8 lines will draw the track of the charges. It is recommended to comment it out, if you use more then 100 charges, since
     *  it will  greately reduce the performance and absorbs memmory.
     */
        AllPixDigitAnimation aDisplay(3,3,50*um,25*um,25*um,2*ntracks,0);

        for(int i=0;i<ntracks;++i){
                aDisplay.AddTrack(I.GetElectron(i)->vX,I.GetElectron(i)->vY,I.GetElectron(i)->vZ,I.GetElectron(i)->vT,1,1,kBlue);
                aDisplay.AddTrack(I.GetHole(i)->vX,I.GetHole(i)->vY,I.GetHole(i)->vZ,I.GetHole(i)->vT,1,1,kRed);
        }

        aDisplay.ViewTracks();



    /** If you want to keep the picture running you need this on. If you comment allpixdigitanimation out. Don't use it. It will only
     *  slow down the program
     */

        rootapp->Run();

        delete Geom;
        delete E;
        delete rand;

        return 0;
  }



