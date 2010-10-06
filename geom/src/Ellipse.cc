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
 double xi,yi,xi_r,yi_r,a_grid_test,y_grid_mint,x_grid_mint,x1_min,y1_min,x_grid_maxt,y_grid_maxt,x1_max,y1_max,theta_grid_min,theta_grid_max;
double pi;
pi=3.14159265;
xi=pt.x-x;
yi=pt.y-y;

xi_r=xi*cos(phi) - yi*sin(phi);
yi_r=xi*sin(phi) + yi*cos(phi);

if (xi_r/(a-r2)>1 || xi_r/(a-r2)<-1)
{
theta_grid_min=0;
}
else
{
theta_grid_min=acos(xi_r/(a-r2));
}


theta_grid_max=acos(xi_r/(a+r2));

a_grid_test=sqrt((xi_r)*(xi_r)+(yi_r)*(yi_r));

x_grid_mint=(a - r2) * cos(theta_grid_min);
y_grid_mint=(a2 - r) * sin(theta_grid_min);

x1_min=x_grid_mint*cos(phi) - y_grid_mint*sin(phi);
y1_min=x_grid_mint*sin(phi) + y_grid_mint*cos(phi);

x_grid_maxt=(a + r2) * cos(theta_grid_max);
y_grid_maxt=(a2 + r) * sin(theta_grid_max);

x1_max=x_grid_maxt*cos(phi) - y_grid_maxt*sin(phi);
y1_max=x_grid_maxt*sin(phi) + y_grid_maxt*cos(phi);


  if (pt.z-z > height || pt.z < z || a_grid_test>a+r2 || a_grid_test<a2-r
                      || (a_grid_test>sqrt(x1_max*x1_max+y1_max*y1_max) || a_grid_test<sqrt(x1_min*x1_min+y1_min*y1_min))
                      || ((pt.z-z)>height*sin(((a_grid_test-sqrt(x1_min*x1_min+y1_min*y1_min))/(sqrt(x1_max*x1_max+y1_max*y1_max)-sqrt(x1_min*x1_min+y1_min*y1_min)))*pi)))
    return false;
  else
    return true;
}

