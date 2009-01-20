// $Id: LatticeUtil.h,v 1.3 2004/06/11 00:48:52 jshumwa Exp $
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
#ifndef __LatticeUtil_h_
#define __LatticeUtil_h_

class CoordsH5;
class LatticeIndexingH5;
class NeighborsH5;
class SpeciesH5;
class SuperCellH5;

/** Base class of lattice generation utilities.  This has a protected
    constructor, so you can only instantiate subclass objects. 
    @author John Shumway */
class LatticeUtil {
protected:
  LatticeUtil();
public:
  /** Destructor deletes supercell. */
  ~LatticeUtil();
  /** Return a reference to the supercell. */
  SuperCellH5& getSuperCell() const {return *cell;}
  /** Return a reference to the species info. */
  SpeciesH5& getSpecies() const {return *species;}
  /** Return a reference to the lattice indexing info. */
  LatticeIndexingH5& getLatticeIndexing() const {return *indexing;}
  /** Return a reference to the atomic coordinates. */
  CoordsH5& getCoords() const {return *coords;}
  /** Return a reference to the nearest neighbor lists. */
  NeighborsH5& getNeighbors() const {return *neighbor;}
protected:
  /** Pointer to the supercell. */
  SuperCellH5* cell;
  /** Pointer to the species info. */
  SpeciesH5* species;
  /** Pointer to the lattice indexing info. */
  LatticeIndexingH5* indexing;
  /** Pointer to the atomic coordinates. */
  CoordsH5* coords;
  /** Pointer to the nearest-neighbor lists. */
  NeighborsH5* neighbor;
};
#endif
