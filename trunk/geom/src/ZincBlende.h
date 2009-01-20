// $Id: ZincBlende.h,v 1.3 2004/06/10 19:41:32 jshumwa Exp $
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
#ifndef __ZincBlende_h_
#define __ZincBlende_h_

#include <string>
#include "LatticeUtil.h"
class SuperCell;
class Nanostructure;
class LatticeIndexingH5;

/** Zinc-blende lattice generation utility.  The anion (i.e. As) is
    located at (0,0,0) and the cation (i.e. Ga) is located at (0.25,0.25,0.25).
    @author John Shumway */
class ZincBlende : public LatticeUtil {
public:
  /** Constructor. */
  ZincBlende(const SuperCell& cell, const Nanostructure& ns);

  /** Directly convert the zinc-blend lattice indexing into neighbor lists.
  This is done by working in cubic, eight-atom base cells.  
  There are 23 atoms that are nieghbors to these eight atoms
  (Atoms 0-7 are in the base cell, and 8-22 are in neighboring cells):
  @verbatim
   Atom       Coords     Anion/Cation
    0       0    0    0       A
    1      1/2  1/2   0       A
    2      1/4  1/4  1/4      C
    3      3/4  3/4  1/4      A
    4      1/2   0   1/2      A
    5       0   1/2  1/2      A
    6      3/4  1/4  3/4      C
    7      1/4  3/4  3/4      C
    8     -1/4 -1/4  1/4      C
    9     -1/4 -1/4 -1/4      C
   10      1/4 -1/4 -1/4      C
   11      3/4  1/4 -1/4      C
   12      1/4  3/4 -1/4      C
   13       1   1/2  1/2      A
   14      1/2   1   1/2      A
   15       1    1    0       A
   16      1/4 -1/4  3/4      C
   17     -1/4  1/4  3/4      C
   18      3/4 -1/4  1/4      C
   19     -1/4  3/4  1/4      C
   20       1    0    1       A
   21      1/2  1/2   1       A
   22       0    1    1       A
  @endverbatim
   Using this indexing, the neighbors of atoms 0-7 are:
  @verbatim
   Atom     Neighbors
    0       2, 8, 9,10
    1       3, 2,11,12
    2       4, 5, 1, 0
    3      13,14,15, 1
    4       6,16,17, 2
    5       7,18, 2,19
    6      20,21,13, 4
    7      21,22,14, 5
  @endverbatim
   The subroutine generates the neigbors by looping over the atoms twice
   times, to do the following steps.:

   <ol>
   <li>(First loop) Compute base-cell indexing from primative-cell indexing.
   <li>(First loop) Place each atom in up to four 23-atom cell-lists.
   <li>(Second loop) Compute base-cell indexing from primative-cell indexing.
   <li> (Second loop) For each atom, read the four neigbors directly 
         of its cell-list.
   </ol>
   */
  NeighborsH5* findNeighbors(const LatticeIndexingH5& indexing) const; 

protected:
  /** Offsets of two atoms in the primative cell. */
  static const double offset[];
};
#endif
