// $Id: Lens.h,v 1.2 2004/04/22 20:09:44 jshumwa Exp $
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
#ifndef __Lens_h_
#define __Lens_h_
#include "Structure.h"

class Material;

/** Representation for a lens shape quantum dot structure.  A lens
    is a slice off of a solid sphere, and is uniquely defined by
    a base diameter and height.  Points (x,y,z) on the surface satisfy
    the equation:
    \f[ (x-x_0)^2 + (y-y_0)^2  = [h-(z-z_0)][(z-z_0)+d^2/4h]\f]

    @author John Shumway */
class Lens : public Structure {
public:
  /** Constructor.  Provide the x,y,z coordinates for the center
      of the base, the height,  and the base diameter.  All coords are in 
      eight-atom cell lattice units. */
  Lens(const double x, const double y, const double z, 
         const double height, const double diameter, const Material* material);

  /** Is a given point inside the lens? */
  bool isPointInStruct(const Vec3& pt) const;

protected:
  /** The x-coordinate of the center of the base, measured in eight 
      atom lattice units. */
  const double x;
  /** The y-coordinate of the center of the base, measured in eight 
      atom lattice units. */
  const double y;
  /** The z-coordinate of the center of the base, measured in eight 
      atom lattice units. */
  const double z;
  /** The height of the lens, measured in eight atom lattice units. */
  const double height;
  /** The diameter of the base, measured in eight atom lattice units. */
  const double diameter;
  /** The quantity $d^2/4h$, used for determining the surface of the lens. */
  const double d2over4h;
};
#endif
