// $Id: Stress.h,v 1.2 2004/06/24 18:28:41 jshumwa Exp $
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
#ifndef __Stress_h_
#define __Stress_h_

#include <blitz/array.h>
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>

class CoordsH5;
class SpeciesH5;
class NeighborsH5;
class SuperCellH5;
class AtomicPotential;

/** Class for computing stresses on each atom.
@version $Revision: 1.2 $
@author John Shumway */
class Stress {
public:
  typedef blitz::TinyVector<float,3> Vec;
  typedef blitz::TinyVector<float,6> SymMat3;
  static double trace(const SymMat3& m) {return m[0]+m[1]+m[2];}
  static double biaxial(const SymMat3& m) {return 2*m[2]-m[0]-m[1];}
  static void outerProduct(const Vec& v1, const Vec& v2, SymMat3& m) {
    m[0]=v1[0]*v2[0]; m[1]=v1[1]*v2[1]; m[2]=v1[2]*v2[2];
    m[3]=0.5*(v1[1]*v2[2]+v1[2]*v2[1]);
    m[4]=0.5*(v1[0]*v2[2]+v1[2]*v2[0]);
    m[5]=0.5*(v1[0]*v2[1]+v1[1]*v2[0]);
  }
  typedef blitz::Array<SymMat3,1> SMArray1;
  typedef blitz::Array<float,1> Array1;
  /// Constructor.
//  Stress(const CoordsH5& coords, const NeighborsH5& nbr,
//    const SpeciesH5& species, const SuperCellH5& cell, 
//    const AtomicPotential& pot);
  Stress() {}
  /// Compute the stresses.
//  void compute();
  /// Access stresses on atom i.
  SymMat3& operator()(const int i) {return stress(i);}
  /// Access stresses on atom i.
  const SymMat3 operator()(const int i) const {return stress(i);}
  /// Access the stresses vector.
  const SMArray1& getStress() const {return stress;}
  /// Allow I/O access by StressH5.
  friend class StressH5;
  /// Get the number of particles.
  int getNPart() const {return stress.size();}
private:
/*  /// Coordinates.
  const CoordsH5& coords;
  /// Neighbor table.
  const NeighborsH5& nbr;
  /// Species index.
  const SpeciesH5& species;
  /// Supercell info.
  const SuperCellH5& cell;
  /// Potential.
  const AtomicPotential& pot;
  /// Stress. */
  SMArray1 stress;
  /// Atomic volumes.
//  Array1 volume;
};
#endif
