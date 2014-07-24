/*
 * PixelGeometry.cpp
 *
 *  Created on: Aug 9, 2012
 *      Author: tjanssen
 */

#include "PixelGeometry.h"
#include <iostream>
using namespace std;

PixelGeometry::PixelGeometry(double Thickness,int nPixX,int nPixY,double pitchX,double pitchY) {
        // TODO Auto-generated constructor stub
        Lx=pitchX*nPixX;
        Ly=pitchY*nPixY;
        nPx=nPixX;
        nPy=nPixY;
        pitchx=pitchX;
        pitchy=pitchY;
        thickness=Thickness;
        address = new int[2];
        close = new int[2];

}
/** This function returns the address of a particle in the geometry. The z component is unimportant here.
 *  The address is stored in a two dimensional array, so it will return a Cartesian coordinate.
 */
	int* PixelGeometry::WhereAmI(double x, double y , double /*z=0*/)
		{
        double x_, y_;
        x_ = x + Lx/2;
        y_ = y + Ly/2;
        address[0] = (int)floor(x_/pitchx)+1;
        address[1] = (int)floor(y_/pitchy)+1;

        return address;
		}
/** Well the name speaks for itself. The if statement checks if the charge is positioned inside or outside
 *  the sensor. If it is the case it returns true and otherwise it is set false by default.
 */
	bool PixelGeometry::IsInside(double x, double y , double z)
	{
        bool inside = false;
        if (x<=Lx/2 && x>=-Lx/2 && y<=Ly/2 && y>=-Ly/2 && z<=thickness/2 && z>=-thickness/2)
                inside = true;
        return inside;
	}
/** When a particle has left the detector we want to know through what electrode. Usually electrons
 *  should leave at pixels and holes at the back. This is due to the electric field. There are 6 electrodes
 *  we just give an address for each of these electrodes. For pixels for instance I have chosen 1,1.
 */
	int* PixelGeometry::ClosestElectrode(double x, double y , double z)
	{
        double left, right, bot, top, back, pixels;
    /// Here all distanced to each electrode is calculated.
        left   = abs(x + Lx/2);
        right  = abs(Lx/2 - x);
        bot    = abs(y + Ly/2);
        top    = abs(Ly/2 - y);
        back   = abs(z + thickness/2);
        pixels = abs(thickness/2 - z);
    /// At the if statements we check for each distance if it is the smallest compared to all others.
        if(left<right && left<bot && left<top && left<back && left<pixels)
        {
                close[0] = -1;
                close[1] = 0;
        }
        if(right<left && right<bot && right<top && right<back && right<pixels)
        {
                close[0] = 257;
                close[1] = 0;
        }
        if(bot<left && bot<right && bot<top && bot<back && bot<pixels)
        {
                close[0] = 0;
                close[1] = -1;
        }
        if(top<left && top<right && top<bot && top<back && top<pixels)
        {
                close[0] = 0;
                close[1] = 257;
        }
        if(back<left && back<right && back<bot && back<top && back<pixels)
        {
                close[0] = 0;
                close[1] = 0;
        }
        if(pixels<left && pixels<right && pixels<bot && pixels<top && pixels<back)
        {
                close[0] = 1;
                close[1] = 1;
        }
    /** There is of course a small chance that there are two smallest distances that are the same. Therefore this
     *  else statement. In that case we always choose either pixels or back. This will also always be the
     *  case because of the geometry (distance from back to pixels is smaller then from left to right and from bottom to top)
     */
        else
        {
                if(back<pixels)
                {
                close[0] = 0;
                close[1] = 0;
                }
                if(pixels<=back)
                {
                close[0] = 1;
                close[1] = 1;
                }
        }

        return close;
	}

/// Here I delete the memory allocated in the free space.
	PixelGeometry::~PixelGeometry() {
        delete address;
        delete close;
        // TODO Auto-generated destructor stub
	}
