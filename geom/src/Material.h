// $Id: Material.h,v 1.2 2004/04/22 20:09:44 jshumwa Exp $
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
#ifndef __Material_h_
#define __Material_h_

#include <string>
#include <valarray>
class Vec3;

/** Base class of materials (pure binary compounds, random alloys, ordered 
    alloys, etc).  Construct the material, then set the
    integer IDs before calling getAnionIndex or getCationIndex.
    @author John Shumway */
class Material {
public:
  /** Constructor. */
  Material (const std::string& name, const int nspecies);
  /** Return the index for an anion at point p. */
  virtual int getAnionIndex(const Vec3&) const=0;
  /** Return the index for a cation at point p. */
  virtual int getCationIndex(const Vec3&) const=0;
  /** The species names for this material. */
  std::valarray<std::string> speciesName;
  /** The species integer IDs for this material. */
  mutable std::valarray<int> speciesID;
  /** Get the material name. */
  const std::string& getName() const {return name;}
protected:
  /** The material name. */
  const std::string name;
};
#endif
