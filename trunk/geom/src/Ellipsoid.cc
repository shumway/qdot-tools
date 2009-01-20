// $Id: Ellipsoid.cc,v 1.1 2004/05/10 21:56:26 jshumwa Exp $
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
#include "Ellipsoid.h"
#include "Vec3.h"

Ellipsoid::Ellipsoid(const double x, const double y, const double z, const double xr, const double yr, const double zr, const Material* material) :
  Structure(material), x(x), y(y), z(z), xr(xr), yr(yr), zr(zr) {
}

bool Ellipsoid::isPointInStruct(const Vec3& pt) const {
    return (pt.x-x)*(pt.x-x)/(xr*xr) + 
           (pt.y-y)*(pt.y-y)/(yr*yr) +
           (pt.z-z)*(pt.z-z)/(zr*zr) < 1.0;
}
