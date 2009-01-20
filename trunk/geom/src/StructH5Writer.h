// $Id: StructH5Writer.h,v 1.4 2005/01/07 17:47:28 jshumwa Exp $
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
#ifndef __StructH5Writer_h_
#define __StructH5Writer_h_

#include <string>
class SuperCell;
class Nanostructure;
class LatticeUtil;

/** Class for writing a struct.h5 file.
    @author John Shumway */
class StructH5Writer{
public:
  /** Constructor requires a supercell and the structural info. */
  StructH5Writer(const SuperCell&, const Nanostructure&);
  /** Destructor. */
  ~StructH5Writer();
  /// Etch an atom species to create vacuum or surface, constructs an
  /// EtchedLattice object to replace current lattice.
  void etch(const int etchIndex);
  /// Reconstuct the surface.
  void reconstruct();
  /** Write the structure to a struct.h5 file. */
  void write(const std::string& filename);
protected:
  /** The supercell. */
  const SuperCell& cell;
  /** The structure. */
  const Nanostructure& structs;
  /** The atomic structural data for the nanostructure. */
  LatticeUtil* lattice;
};
#endif
