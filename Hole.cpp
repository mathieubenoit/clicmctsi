/*
 * Hole.cpp
 *
 *  Created on: Jul 30, 2012
 *      Author: tjanssen
 */

#include "Hole.h"
/// All properties of Hole are given in Charge. Here its constants are
	Hole::Hole()
	{
	 // TODO Auto-generated constructor stub
		alpha=2.4e7*cm/s;
		tnom=600.0;
		theta=0.8;
		betaexp = 0.66;
		beta_0 = 1.109;
		beta = 2.0;
		mob_0 = 450.0*cm*cm/(V*s);
		TL = 300;
		Diffusion_Coefficent = 12*cm*cm/s;
		Q=1;
		rh=0.7;
	}

