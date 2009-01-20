// $Id: Dome.cc,v 1.2 2004/06/14 23:33:31 jshumwa Exp $
/*
    Copyright (C) 2004 John B. Shumway, Jr.

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
#include "Dome.h"
#include "Vec3.h"
#include <cmath>

Dome::Dome(const double x, const double y, const double z, 
           const FArray& facet, const Material* material) :
  Structure(material), x(x), y(y), z(z), height(1000),
  nfacet(facet.size()), boundingPlane(nfacet) {
  for (int i=0; i<nfacet; ++i) {
    if (height>facet[i].height) height=facet[i].height;
    boundingPlane(i)[0]=facet[i].height;
    boundingPlane(i)[1]=-(1.0*facet[i].i)/(1.0*facet[i].k);
    boundingPlane(i)[2]=-(1.0*facet[i].j)/(1.0*facet[i].k);
  }
}

bool Dome::isPointInStruct(const Vec3& pt) const {
  double deltaz=pt.z-z;
  if (deltaz<0 || deltaz>height) return false;
  double deltax=std::abs(pt.x-x);
  double deltay=std::abs(pt.y-y);
  if (deltax<deltay) {double temp=deltax; deltax=deltay; deltay=temp;}
  for (int i=0; i<nfacet; ++i) {
    if (deltaz > dot(boundingPlane(i),Plane(1.0,deltax,deltay))) return false;
  }
  return true;  
}
