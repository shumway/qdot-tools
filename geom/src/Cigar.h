// $Id: Cigar.h,v 1.1 2004/05/10 21:56:26 jshumwa Exp $
/*
    Copyright (C) 2004 Gabriel Bester and John B. Shumway, Jr.

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
#ifndef __Cigar_h_
#define __Cigar_h_
#include "Structure.h"

class Material;

/** Representation for a cigar shaped nanostructure (two spheres connected
    by a cylinder, arbitrary orientation). 
    @author Gabriel Bester */
class Cigar : public Structure {
public:
  /** Constructor.  Provide the x,y,z coordinates for the center
      of the cigar, the radius, the lenght, the orientation (vector). 
      All coords are in eight-atom cell lattice units. */
  Cigar(const double x, const double y, const double z, 
        const double radius, const double length, 
        const double vx, const double vy, const double vz,
        const Material* material);

  /** Is a given point inside the sphere? */
  bool isPointInStruct(const Vec3& pt) const;

protected:
  /** The x-coordinate of the center of the cigar, measured in eight 
      atom lattice units. */
  const double x;
  /** The y-coordinate of the center of the cigar, measured in eight 
      atom lattice units. */
  const double y;
  /** The z-coordinate of the center of the cigar, measured in eight 
      atom lattice units. */
  const double z;
  /** The radius of the cigar, measured in eight atom lattice units. */
  const double radius;
  /** The lenght of the cigar, measured in eight atom lattice units. */
  const double length;
  /** The x coordinate of the vector along the cigar axis. */
  const double vx;
  /** The y coordinate of the vector along the cigar axis. */
  const double vy;
  /** The z coordinate of the vector along the cigar axis. */
  const double vz;
};
#endif
