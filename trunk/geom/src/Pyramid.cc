// $Id: Pyramid.cc,v 1.2 2004/04/22 20:09:44 jshumwa Exp $
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
#include "Pyramid.h"
#include "Vec3.h"
#include <cmath>

Pyramid::Pyramid(const double x, const double y, const double z, 
                 const double height, const double base, const double top,
                 const Material* material) :
  Structure(material), x(x), y(y), z(z), height(height), base(base), top(top) {
}

bool Pyramid::isPointInStruct(const Vec3& pt) const {
    double deltax=std::abs(pt.x-x);
    double deltay=std::abs(pt.y-y);
    double deltaz=pt.z-z;
    return deltaz >= 0 && deltaz <= height &&
         2*((deltax>deltay)?deltax:deltay)
         <= base+deltaz*(top-base)/height;
}
