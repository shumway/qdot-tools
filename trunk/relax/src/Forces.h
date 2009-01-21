// $Id: Forces.h,v 1.3 2007/03/14 19:47:50 jshumwa Exp $
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
#ifndef __Forces_h_
#define __Forces_h_

#include <valarray>
#include <blitz/tinyvec.h>
class CoordsH5;
class SpeciesH5;
class NeighborsH5;
class SuperCellH5;
class AtomicPotential;

/** Class for computing forces on each atom.
 * @author John Shumway */
class Forces {
public:
#ifdef ENABLE_FLOAT
  typedef blitz::TinyVector<float,3> Vec3;
#else
  typedef blitz::TinyVector<double,3> Vec3;
#endif
  /// Constructor.
  Forces(const CoordsH5& coords, const NeighborsH5& nbr,
    const SpeciesH5& species, const SuperCellH5& cell, 
    const AtomicPotential& pot);
  /// Compute the forces.
  void compute();
  /// Access force on atom i.
  Vec3& operator[](const int i) {return force[i];}
  /// Access force on atom i.
  const Vec3& operator[](const int i) const {return force[i];}
  /// Access the force vector.
  const std::valarray<Vec3>& getForces() const {return force;}
private:
  /// Coordinates.
  const CoordsH5& coords;
  /// Neighbor table.
  const NeighborsH5& nbr;
  /// Species index.
  const SpeciesH5& species;
  /// Supercell info.
  const SuperCellH5& cell;
  /// Stillenger-Weber potential.
  const AtomicPotential& pot;
  /// Forces.
  std::valarray<Vec3> force;
};
#endif
