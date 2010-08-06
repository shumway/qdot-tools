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
#ifndef Ring_H
#define Ring_H

#include <cmath>
#include "Vec3.h"
#include "Structure.h"

class Material;

/**The ring can be seen as a circle with a radius of radius in making a 
cycling movement around the axes Z=z (paralelling the zaxe) with a radius of 
radius out. And then the Ring is the part that Pt.z>z) */
class Ring : public Structure {
protected:
  const double x;
  const double y;
  const double z;
  const double height;
  const double radiusin;
  const double radiusout;
  const double s11;
public:
  Ring(const double x, const double y, const double z,const double height,
           const double radiusin, const double radiusout, const double s11,
           const Material* material);
  bool isPointInStruct(const Vec3& pt) const;
};
#endif
