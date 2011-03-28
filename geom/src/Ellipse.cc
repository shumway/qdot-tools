/*
    Peter G McDonald, 2010, Heriot Watt University

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
#include "Ellipse.h"
#include "Vec3.h"
#include <cmath>

Ellipse::Ellipse(const double x, const double y, const double z, const double phi,
const double height, const double a, const double a2, const double r2, const double r,
const Material* material)
  : Structure(material), x(x), y(y), z(z), height(height), phi(phi),
    a(a), a2(a2), r2(r2), r(r) {
}

bool Ellipse::isPointInStruct(const Vec3 &pt) const {
double xi,yi,xi_r,yi_r,r_ijk,r_min,r_max,zeta;
double pi;
pi=3.14159265;
xi=pt.x-x;
yi=pt.y-y;

xi_r=xi*cos(phi) - yi*sin(phi);
yi_r=xi*sin(phi) + yi*cos(phi);

r_ijk=sqrt((xi_r)*(xi_r)+(yi_r)*(yi_r));

zeta=atan2(yi_r,xi_r);
                
r_max=((a+r)*(a2+r))/sqrt((pow((a2+r)*cos(zeta),2)+(pow((a+r)*sin(zeta),2))));
r_min=((a-r)*(a2-r))/sqrt((pow((a2-r)*cos(zeta),2)+(pow((a-r)*sin(zeta),2))));

  if (pt.z-z > height || pt.z < z
                      || r_ijk>r_max || r_ijk<r_min
                      || ((pt.z-z)>=height*sin(((r_ijk-r_min)/(r_max-r_min))*pi)))
    return false;
  else
    return true;
}

