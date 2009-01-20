// $Id: Binary.h,v 1.2 2004/04/22 20:09:44 jshumwa Exp $
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
#ifndef __Binary_h_
#define __Binary_h_

#include "Material.h"

/** Material class for a pure binary compond.  Construct the material, 
    then set the integer IDs before calling getAnionIndex or getCationIndex.  
    @author John Shumway */
class Binary : public Material {
public:
  /** Constructor. */
  Binary(const std::string& name, const std::string& anionName,
                                  const std::string& cationName);
  /** Return the index for an anion at point p. */
  virtual int getAnionIndex(const Vec3& p) const;
  /** Return the index for a cation at point p. */
  virtual int getCationIndex(const Vec3& p) const;
};
#endif
