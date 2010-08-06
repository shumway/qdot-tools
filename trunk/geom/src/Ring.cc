// $Id: Lens.h,v 1.2 2004/04/22 20:09:44 jshumwa Exp $
/*
    Copyright (C) 2008 Lixin He

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
#include "Ring.h"
#include "Vec3.h"
#include <cmath>

Ring::Ring(const double x, const double y, const double z,
  const double height, const double radiusin, const double radiusout, const double s11,
  const Material *material)
  : Structure(material), x(x), y(y), z(z), height(height), 
    radiusin(radiusin), radiusout(radiusout), s11(s11) {
}

bool Ring::isPointInStruct(const Vec3 &pt) const {
  double a, b, c, d, r, x1, y1, x2, y2;
  x1=0.5*pt.x*((1/s11)+1)+0.5*pt.y*((1/s11)-1);
  y1=0.5*pt.x*((1/s11)-1)+0.5*pt.y*((1/s11)+1);
  y1=0.5*pt.x*((1/s11)-1)+0.5*pt.y*((1/s11)+1);

  x2=0.5*x*((1/s11)+1)+0.5*y*((1/s11)-1);
  y2=0.5*x*((1/s11)-1)+0.5*y*((1/s11)+1);

  a = sqrt((x1-x2) * (x1-x2) + (y1-y2) * (y1-y2));

  b = (radiusout - radiusin)/2;
  d = (radiusout + radiusin)/2;
  r = 0.5 * (height * height + b * b)/height; 
  c = r - height;   

  if (pt.z-z > height || pt.z < z 
                      || r * r  < (a-d) * (a-d) + (pt.z+c-z) * (pt.z+c-z))
    return false;
  else
    return true;
}

