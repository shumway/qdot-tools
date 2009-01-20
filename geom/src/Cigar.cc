// $Id: Cigar.cc,v 1.1 2004/05/10 21:56:26 jshumwa Exp $
/*
    Copyright (C) 2004 Gabriel Bester

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
#include "Cigar.h"
#include "Vec3.h"
#include <cmath>
Cigar::Cigar(const double x, const double y, const double z, 
             const double radius, const double length, 
             const double vx, const double vy, const double vz,
             const Material* material) :
  Structure(material), x(x), y(y), z(z), radius(radius), length(length), vx(vx), vy(vy), vz(vz){}

bool Cigar::isPointInStruct(const Vec3& pt) const {
  double v2 = vx*vx+vy*vy+vz*vz;
  double r2 = radius*radius;
  double d = (length-2*radius)/2;
  double rlarge2 = d*d + radius*radius;
  double part1 = pt.x-x-(vx*(pt.x-x)+vy*(pt.y-y)+vz*(pt.z-z))*vx/v2;
  double part2 = pt.y-y-(vx*(pt.x-x)+vy*(pt.y-y)+vz*(pt.z-z))*vy/v2;
  double part3 = pt.z-z-(vx*(pt.x-x)+vy*(pt.y-y)+vz*(pt.z-z))*vz/v2;
  double c1x = x + d*vx/sqrt(v2);
  double c1y = y + d*vy/sqrt(v2);
  double c1z = z + d*vz/sqrt(v2);
  double c2x = x - d*vx/sqrt(v2);
  double c2y = y - d*vy/sqrt(v2);
  double c2z = z - d*vz/sqrt(v2);
  return ((pt.x-x)*(pt.x-x)+(pt.y-y)*(pt.y-y)+(pt.z-z)*(pt.z-z) < rlarge2 &&
  part1*part1 + part2*part2 + part3*part3 < r2) ||
  (pt.x-c1x)*(pt.x-c1x)+(pt.y-c1y)*(pt.y-c1y)+(pt.z-c1z)*(pt.z-c1z) < r2 ||
  (pt.x-c2x)*(pt.x-c2x)+(pt.y-c2y)*(pt.y-c2y)+(pt.z-c2z)*(pt.z-c2z) < r2;
}
