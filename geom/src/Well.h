// $Id: Well.h,v 1.2 2004/04/22 20:09:44 jshumwa Exp $
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
#ifndef __Well_h_
#define __Well_h_
#include "Structure.h"

class Material;

/** Representation for a quantum well structure.
    @author John Shumway */
class Well : public Structure {
public:
  /** Constructor.  Provide the z coordinate for the center
      of the well and the thickness.  All coords are in 
      eight-atom cell lattice units. */
  Well(const double z, const double thickness, const Material* material);

  /** Is a given point inside the well? */
  bool isPointInStruct(const Vec3& pt) const;

protected:
  /** The center of the well, measured in eight atom lattice units. */
  const double z;
  /** The thickness of the well, measured in eight atom lattice units. */
  const double thickness;
};
#endif
