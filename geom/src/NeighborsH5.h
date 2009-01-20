// $Id: NeighborsH5.h,v 1.3 2004/06/10 19:41:32 jshumwa Exp $
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
#ifndef __NeighborsH5_h_
#define __NeighborsH5_h_

#include <valarray>
#include <string>
#include "StructH5.h"

/** Nearest neighbors.
 
    Layout in the H5 file:
    <ul><li>neighbors
      - neighbor (dataset, integer[nAtom,maxNeighbors]) - 
          Neighbor tables, with atoms enumerated starting from 1.
          Missing neighbor indicated by "0". </ul>
    @author John Shumway */
class NeighborsH5 : public StructH5 {
public:
  /** Constructor */
  NeighborsH5(const int nmax, const int natoms);
  /** The nearest neighbor lists. */
  mutable std::valarray<int> nn;
  /** The maximum number of neighbors per atoms. */
  int nmax;
  /** Write the supercell to a struct.h5 file. */
  void h5Write(const std::string& filename, const int mode=APPEND) const;
};
#endif
