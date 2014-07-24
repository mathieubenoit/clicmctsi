/*
 * Electron.cpp
 *
 *  Created on: Jul 30, 2012
 *      Author: tjanssen
 */

#include "Electron.h"
/// All properties of Electron are given in Charge. Here its constants are defined.
	Electron::Electron()
	{
	 // TODO Auto-generated constructor stub
		alpha=2.4e7*cm/s;
		tnom=600.0;
    	theta=0.8;
		betaexp = 0.17;
		beta_0 = 1.213;
    	beta = 1.0;
    	mob_0 = 1400.0*cm*cm/(V*s);
    	TL = 300;
    	Diffusion_Coefficent = 36*cm*cm/s;
    	Q=-1;
    	rh=1.15;
	}



