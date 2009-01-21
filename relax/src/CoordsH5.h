// $Id: CoordsH5.h,v 1.3 2007/03/14 19:47:50 jshumwa Exp $
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
#ifndef __CoordsH5_h_
#define __CoordsH5_h_

#include <valarray>
#include <string>
#include "StructH5.h"
#include <blitz/tinyvec.h>

/** The coordinates.
 
    Layout in the H5 file:
    \begin{itemize}
    \item coords \begin{itemize}
      \item nAtom (attribute, integer) - Number of atoms.
      \item coord (dataset, real[nAtom,3]) - Positions of the atoms.
    \end{itemize}\end{itemize}
    @author John Shumway */
class CoordsH5 : public StructH5 {
public:
#ifdef ENABLE_FLOAT
  typedef blitz::TinyVector<float,3> Vec3;
#else
  typedef blitz::TinyVector<double,3> Vec3;
#endif
  /** Constructor */
  CoordsH5(const int natom);
  /** Construct from a struct.h5 file. */
  CoordsH5(const std::string& filename);
  /** The coordinates. */
  mutable std::valarray<Vec3> coords;
  /** Get coordinates for atom i. */
  Vec3& operator[](const int i) {return coords[i];}
  /** Get coordinates for atom i. */
  const Vec3& operator[](const int i) const {return coords[i];}
  /** Write the coordinates to a struct.h5 file. */
  void h5Write(const std::string& filename, const int mode=APPEND) const;
  /** Read the coordinates from a struct.h5 file. */
  void h5Read(const std::string& filename);
};
#endif
