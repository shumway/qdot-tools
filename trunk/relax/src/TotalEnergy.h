// $Id: TotalEnergy.h,v 1.4 2007/03/12 20:07:11 jshumwa Exp $
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
#ifndef __TotalEnergy_h_
#define __TotalEnergy_h_

class CoordsH5;
class SpeciesH5;
class NeighborsH5;
class SuperCellH5;
class AtomicPotential;

/** Class for computing the total energy.
 * @author John Shumway */
class TotalEnergy {
public:
  /// Constructor.
  TotalEnergy(const CoordsH5& coords, const NeighborsH5& nbr,
    const SpeciesH5& species, const SuperCellH5& cell, 
    const AtomicPotential& pot);
  /// Compute the total energy.
  void compute();
  /// Report details of the total energy.
  void report();
  /// Get energy.
  double getEnergy() {return energy;}
private:
  /// Coordinates.
  const CoordsH5& coords;
  /// Neighbor table.
  const NeighborsH5& nbr;
  /// Species index.
  const SpeciesH5& species;
  /// Supercell info.
  const SuperCellH5& cell;
  /// Atomic potential.
  const AtomicPotential& pot;
  /// Energy.
  double energy;
  /// Stretch energy.
  double stretchEnergy;
  /// Bend energy.
  double bendEnergy;
  /// Binding energy.
  double bindingEnergy;
  /// Number of threads (if openmp).
  //int numThreads;
  /// Displacements to neighbors.
  //blitz::Array<Vec3,2> delta;
  /// Distances to neighbors.
  //blitz::Array<double,2> d;
  /// Species of neighbors.
  //blitz::Array<int,2> spec;
  /// Is i greater than j for the atom (used to avoid double
  /// counting bonds).
  //blitz::Array<bool,2> igej;
};
#endif
