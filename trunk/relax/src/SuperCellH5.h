// $Id: SuperCellH5.h,v 1.3 2007/03/14 19:47:50 jshumwa Exp $
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
#include <blitz/tinyvec.h>
#include "StructH5.h"

/** The supercell.

    Layout in the H5 file:
    \begin{itemize}
    \item superCell \begin{itemize}
      \item a1 (attribute, real[3]) - 1st supercell lattice vector
      \item a2 (attribute, real[3]) - 2nd supercell lattice vector
      \item a3 (attribute, real[3]) - 3rd supercell lattice vector
    \end{itemize}\end{itemize}
    @author John Shumway */
class SuperCellH5 : public StructH5 {
public:
#ifdef ENABLE_FLOAT
  typedef blitz::TinyVector<float,3> Vec3;
#else
  typedef blitz::TinyVector<double,3> Vec3;
#endif
  /** Constructor */
  SuperCellH5();
  /** Construct from file. */
  SuperCellH5(const std::string& filename);
  /** Supercell lattice vector $a_1$. */
  mutable Vec3 a1;
  /** Supercell lattice vector $a_2$. */
  mutable Vec3 a2;
  /** Supercell lattice vector $a_3$. */
  mutable Vec3 a3;
  /** Read the supercell from a struct.h5 file. */
  void h5Read(const std::string& filename);
  /** Write the supercell to a struct.h5 file. */
  void h5Write(const std::string& filename, const int mode=APPEND) const;
  /** Compute reciprical lattice vectors. */
  void computeRecipricalVectors();
  /** Set to the smallest displacement with PBC.
   * @bug Only projects to primative unit cell, not-neccesarily
   * the smallest vector. */
  void pbc(Vec3&) const;
private:
  /** Reciprical supercell lattice vector $b_1,b_2,b_3$. */
  Vec3 b1,b2,b3;
};
#endif
