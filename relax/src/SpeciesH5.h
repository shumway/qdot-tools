// $Id: SpeciesH5.h,v 1.1.1.1 2004/04/22 23:44:54 jshumwa Exp $
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
#ifndef __SpeciesH5_h_
#define __SpeciesH5_h_

#include <valarray>
#include <string>
#include "StructH5.h"

/** The species (atom types).

    Layout in the H5 file:
    \begin{itemize}
    \item species \begin{itemize}
      \item name (attribute, string[nSpecies]) - Species names.
      \item nSpecies (attribute, integer) - Number of species.
      \item species (attribute, integer[nPart]) - Species type of each atom
              (enumerated as integers starting from 1).
    \end{itemize}\end{itemize}
    @author John Shumway */
class SpeciesH5 : public StructH5 {
public:
  /** Constructor */
  SpeciesH5(const int nspecies, const int natoms);
  /** Construct from a file */
  SpeciesH5(const std::string& filename);
  /** The species names. */
  std::valarray<std::string> name;
  /** The species indexing. */
  mutable std::valarray<int> species;
  /** Read the supercell from a struct.h5 file. */
  void h5Read(const std::string& filename);
  /** Write the supercell to a struct.h5 file. */
  void h5Write(const std::string& filename, const int mode=APPEND) const;
};
#endif
