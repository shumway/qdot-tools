// $Id: Cylinder.h,v 1.1 2004/05/10 21:56:26 jshumwa Exp $
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
#ifndef __Cylinder_h_
#define __Cylinder_h_
#include "Structure.h"

class Material;

/** Representation for a cylindrical quantum wire structure.
    @author Gabriel Bester */
class Cylinder : public Structure {
public:
  /** Representatation for a Cylinder (quantum wire). The only
      Implemented groth direction of the wire is the z-direction. 
      The wire is infinitely long. REM: composition flucutuation
      should be possible using composition gradients (not tested)*/
  Cylinder(const double x, const double y, const double radius, 
           const Material* material);

  /** Is a given point inside the cylinder? */
  bool isPointInStruct(const Vec3& pt) const;

protected:
  /** The x-coordinate of the center of the cylinder, measured in eight 
      atom lattice units. */
  const double x;
  /** The y-coordinate of the center of the cylinder, measured in eight 
      atom lattice units. */
  const double y;
  /** The radius of the cylinder, measured in eight atom lattice units. */
  const double radius;
  /** Radius squared */
  const double r2;
};
#endif
