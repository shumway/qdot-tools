// $Id: Lens.h,v 1.2 2004/04/22 20:09:44 jshumwa Exp $
/*

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/
#ifndef Ellipse_H
#define Ellipse_H

#include <cmath>
#include "Vec3.h"
#include "Structure.h"

class Material;

/**The Ellipse can be seen as a circle with a radius of radius in making a 
cycling movement around the axes Z=z (paralelling the zaxe) with a radius of 
radius out. And then the Ellipse is the part that Pt.z>z) */
class Ellipse : public Structure {
protected:
  const double x;
  const double y;
  const double z;
  const double phi;
  const double height;
  const double a;
  const double a2;
  const double r2;
  const double r;
public:
  Ellipse(const double x, const double y, const double z,const double phi,
		const double height, const double a, const double a2, const double r2, const double r,
           	const Material* material);
  bool isPointInStruct(const Vec3& pt) const;
};
#endif
