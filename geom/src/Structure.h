// $Id: Structure.h,v 1.2 2004/04/22 20:09:44 jshumwa Exp $
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
#ifndef __Structure_h_
#define __Structure_h_

class SpeciesH5;
class LatticeIndexingH5;
class Material;
class Vec3;

/** Base class of all structures.  Structures may be combined to
    generate a nanostructure.  When structures overlap, each call
    to ``build'' overwrites any overlapping structure.

    Structures have pointers materials.
    @author John Shumway */
class Structure {
public: 
  /** Constructor. */
  Structure(const Material* material);
  /** Lay down the structure in the cell. */
  void build(SpeciesH5&, const LatticeIndexingH5&) const;
  /** Is the given x,y,z point in the structure? */
  virtual bool isPointInStruct(const Vec3& pt)const=0;
protected:
  /** The material for the strucutre. */
  const Material* material;
};
#endif
