// $Id: BinaryAlloy.h,v 1.2 2004/04/22 20:09:44 jshumwa Exp $
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
#ifndef __BinaryAlloy_h_
#define __BinaryAlloy_h_

#include "Material.h"

/** Material class for a binary alloy.  Construct the material, 
    then set the integer IDs before calling getAnionIndex or getCationIndex.  
    Composition is randomly generated, and repeated calls to the same atom
    site will generate random results. 
    @author John Shumway */
class BinaryAlloy : public Material {
public:
  /** Constructor.  The last argument x1 = (1-x2) is the fraction of 
      atom type one.*/
  BinaryAlloy(const std::string& name, const std::string& atom1, 
               const std::string& atom2, const double x1);
  /** Return the index for an anion at point p. */
  virtual int getAnionIndex(const Vec3& p) const;
  /** Return the index for a cation at point p.  Not reproducable,
      and repeated calls will generate different results. */
  virtual int getCationIndex(const Vec3& p) const;
private:
  /// The percentage of atoms of type one.
  double x1;
};
#endif
