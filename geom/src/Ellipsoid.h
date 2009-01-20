// $Id: Ellipsoid.h,v 1.1 2004/05/10 21:56:26 jshumwa Exp $
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
#ifndef __Ellipsoid_h_
#define __Ellipsoid_h_
#include "Structure.h"

class Material;

/** Representation for a cone shape quantum dot structure.  
    @author Gabriel Bester */
class Ellipsoid : public Structure {
public:
  /** Constructor.  Provide the x,y,z coordinates for the center
      and rx,yr,zr for the radii in x,y,z direction.  All coords are in 
      eight-atom cell lattice units. */
  Ellipsoid(const double x,  const double y,  const double z, 
            const double xr, const double yr, const double zr,
            const Material* material);

  /** Is a given point inside the cone? */
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
  /** The radius in x direction. */
  const double xr;
  /** The radius in y direction. */
  const double yr;
  /** The radius in z direction. */
  const double zr;
};
#endif
