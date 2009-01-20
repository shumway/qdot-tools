// $Id: Hut.h,v 1.3 2004/05/10 21:56:26 jshumwa Exp $
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
#ifndef __Hut_h_
#define __Hut_h_
#include "Structure.h"

class Material;

/** Representation for a rectangluar Hut shaped quantum dot structure. 
    @author John Shumway, Paul Logan */
class Hut : public Structure {
public:
  /** Constructor.  Provide the x,y,z coordinates for the center of the base, 
      the height, the base-x length, the base-y length, the x-top edge length,
      and the y-top edge length.  Set top=0 for a non-truncated hut.  All 
      hcoords are in eight-atom cell lattice units. */
  Hut(const double x, const double y, const double z, const double height, 
      const double xbase, const double ybase, const double xtop, 
      const double ytop, const Material* material);

  /** Is a given point inside the hut? */
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
  /** The height of the hut, measured in eight atom lattice units. */
  const double height;
  /** The length of the x-base edge, measured in eight atom lattice units. */
  const double xbase;
  /** The length of the y-base edge, measured in eight atom lattice units. */
  const double ybase;
  /** The length of the top edge, in the x direction, measured in eight atom 
      lattice units.  This is used for truncated square huts, use top=0 for 
      non-truncated huts. */
  const double xtop;
  /** The length of the top edge, in the y direction, measured in eight atom 
      lattice units.  This is used for truncated square huts, use top=0 for 
      non-truncated huts. */
  const double ytop;
    
};
#endif
