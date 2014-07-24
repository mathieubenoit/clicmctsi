/*
 * PixelGeometry.h
 *
 *  Created on: Aug 9, 2012
 *      Author: tjanssen
 */

#ifndef PIXELGEOMETRY_H_
#define PIXELGEOMETRY_H_

#include <math.h>
#include "TMath.h"
/** The class PixelGeometry defines the edges of the sensor. It determines if the particle is inside or outside
 *  the sensor and at what the nearest electrode to a particle is at all given times.
 */
class PixelGeometry {
public:
	/** In the constructor of pixel geometry we determine the pixel sizes and the number of pixels.
	 *  first we determine the pixel thickness, then the number of pixels in x, then the number of pixels in y
	 *  In the end the pixel size in x and then the pixel size in y.
	 */
        PixelGeometry(double Thickness,int nPixX,int nPixY,double pitchX,double pitcxhY);
    /// Destructor
        virtual ~PixelGeometry();
    /// Return (i,j) address of the pixel volume (x,y,z) is in
        int* WhereAmI(double x, double y , double /*z=0*/);
    /// Returns the thickness of the sensor
        double GetThickness(){return thickness;}
    /// Return true if (x,y,z) is inside the sensor
        bool IsInside(double x, double y , double z);
    /** Return a (i,j) value of the closest electrode. (0,0) for back side, (-1,0) for left edge
     *  (257,0) for right edge, (0,-1) for bottom edge and (0,257) for top edger, (1,1) for pixels.
     */
        int* ClosestElectrode(double x, double y , double z);

private :

        double Lx,Ly;
        double pitchx,pitchy;
        double thickness;
        int nPx, nPy;
        int *address;
        int *close;


};

#endif /* PIXELGEOMETRY_H_ */

