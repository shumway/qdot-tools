// $Id: SuperCellH5.h,v 1.3 2004/06/10 19:41:32 jshumwa Exp $
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
#ifndef __SuperCellH5_h_
#define __SuperCellH5_h_

#include <string>
#include "StructH5.h"

/** The supercell.

    Layout in the H5 file:
    <ul><li>superCell
      - a1 (attribute, real[3]) - 1st supercell lattice vector
      - a2 (attribute, real[3]) - 2nd supercell lattice vector
      - a3 (attribute, real[3]) - 3rd supercell lattice vector
    </li></ul>
    @author John Shumway */
class SuperCellH5 : public StructH5 {
public:
  /** Constructor */
  SuperCellH5();
  /** Supercell lattice vector $a_1$. */
  mutable double a1[3];
  /** Supercell lattice vector $a_2$. */
  mutable double a2[3];
  /** Supercell lattice vector $a_3$. */
  mutable double a3[3];
  /** Write the supercell to a struct.h5 file. */
  void h5Write(const std::string& filename, const int mode=APPEND) const;
};
#endif
