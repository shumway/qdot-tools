// $Id: SurfaceStress.h,v 1.2 2007/03/14 19:47:50 jshumwa Exp $
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
#ifndef __SurfaceStress_h_
#define __SurfaceStress_h_

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
class SurfaceStress {
public:
#ifdef ENABLE_FLOAT
  typedef blitz::TinyVector<float,3> Vec3;
#else
  typedef blitz::TinyVector<double,3> Vec3;
#endif
  static const int XX=0,YY=1,ZZ=2,YZ=3,XZ=4,XY=5;
  /// Internal class for storing six components of stress.
#ifdef ENABLE_FLOAT
  class SymMat : public blitz::TinyVector<float,6> {
#else
  class SymMat : public blitz::TinyVector<double,6> {
#endif
  public:
    /// Set matrix value to symmetrized outer product of two vectors.
    SymMat& tensorProduct(const Vec3& v1, const Vec3& v2) {
      (*this)[XX]=v1[0]*v2[0]; (*this)[YY]=v1[1]*v2[1];
      (*this)[ZZ]=v1[2]*v2[2];
      (*this)[YZ]=0.5*(v1[1]*v2[2]+v1[2]*v2[1]);
      (*this)[XZ]=0.5*(v1[0]*v2[2]+v1[2]*v2[0]);
      (*this)[XY]=0.5*(v1[0]*v2[1]+v1[1]*v2[0]);
      return *this;
    }
  };
  typedef blitz::Array<SymMat,1> SMArray1;
  typedef blitz::Array<double,1> Array1;
  /// Constructor.
  SurfaceStress(const CoordsH5& coords, const NeighborsH5& nbr,
    const SpeciesH5& species, const SuperCellH5& cell, 
    const AtomicPotential& pot);
  /// Compute the stresses.
  void compute();
  /// Access stresses on atom i.
  SymMat& operator()(const int i) {return stress(i);}
  /// Access stresses on atom i.
  const SymMat operator()(const int i) const {return stress(i);}
  /// Access the stresses vector.
  const SMArray1& getStress() const {return stress;}
  /// Allow I/O access by StressH5.
  friend class StressH5;
private:
  /// Coordinates.
  const CoordsH5& coords;
  /// Neighbor table.
  const NeighborsH5& nbr;
  /// Species index.
  const SpeciesH5& species;
  /// Supercell info.
  const SuperCellH5& cell;
  /// Potential.
  const AtomicPotential& pot;
  /// Stress.
  SMArray1 stress;
  /// Atomic volumes.
  Array1 volume;
};
#endif
