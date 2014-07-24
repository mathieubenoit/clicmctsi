#include <iostream>
#include "TString.h"
#include "TMath.h"
/** ifndef and endif CONSTANT_H_ makes sure that if I include libraries twice. Meaning in here and a class
 *  that inherits CONSTANT. The library will only be copied ones. So I can be lazy and not check this without
 *  performance loss. If I call it twice than a large algorithm will be copied while compiling the code after all
 *  this construction is made in each class.
 */
	#ifndef CONSTANT_H_
	#define CONSTANT_H_
/// in here I define the unit system and the constant Pi. Without it one easily makes unit errors
	const double cm=10000.;
	const double m =100*cm;
	const double um=1.;
	const double s = 1.e9;
	const double us= 1.e3;
	const double ns= 1.;

	const double keV =1.;
	const double MeV=1000.;

	const double V=1;
	const double uA=1.;
	const double mA=1000.;
	const double A=1.e6;

	const double T=V*s/m/m;
/// Here I define pi from the root directory since its very precise. The one from C++ is only 3.14
	const double pi = TMath::Pi();

	#endif


