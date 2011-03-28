// $Id: Strain.h,v 1.3 2007/12/06 23:59:59 jshumwa Exp $
/*
    Copyright (C) 2007 John B. Shumway, Jr.

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
#ifndef __Strain_h_
#define __Strain_h_
class GridFactory;
class CoordsH5;
class SuperCellH5;
class NeighborsH5;
class CompositionGrids;
#include <blitz/tinyvec.h>
#include <blitz/tinymat.h>
#include "Grid.h"
#include "StrainGrid.h"

class Strain: public StrainGrid {
public:
  /// Typedefs.
  typedef blitz::TinyVector<int,3> IVec;
  typedef blitz::TinyVector<float,3> Vec;
  typedef blitz::TinyVector<float,6> SymMat3;
  typedef Grid<SymMat3> SMGrid;
  static double trace(const SymMat3& m) {return m[0]+m[1]+m[2];}
  static double biaxial(const SymMat3& m) {return 2*m[2]-m[0]-m[1];}
  /// Constructor.
  Strain(const GridFactory&, const CoordsH5&, 
         const CompositionGrids&, const NeighborsH5&, const SuperCellH5& cell, bool isInGaAs);
  /// Virtual destructor.
  virtual ~Strain();
  /// Get reference to the grid.
  virtual SMGrid& getGrid() {return *grid;}
  /// Get const reference to the grid.
  virtual const SMGrid& getGrid() const {return *grid;}
  /// Get reference to the trace grid.
  virtual Grid<float>& getTraceGrid() {return *traceGrid;}
  /// Get const reference to the trace grid.
  virtual const Grid<float>& getTraceGrid() const {return *traceGrid;}
  /// Get reference to the trace grid.
  virtual Grid<float>& getBiaxialGrid() {return *biaxialGrid;}
  /// Get const reference to the trace grid.
  virtual const Grid<float>& getBiaxialGrid() const {return *biaxialGrid;}
  /// Get reference to a grid element.
  virtual SymMat3& operator()(const int i, const int j, const int k) {
                              return (*grid)(i,j,k);}
  /// Get const reference to a grid element.
  virtual const SymMat3& operator()(const int i, const int j, const int k) 
                   const {return (*grid)(i,j,k);}
private:
  /// The grid.
  SMGrid* grid;
  /// The grid for the trace of stress.
  Grid<float>* traceGrid;
  /// The grid for the biaxial stress.
  Grid<float>* biaxialGrid;
  /// Coordinates.
  const CoordsH5 &coords;
  /// Coordinates.
  const NeighborsH5 &neighbors;
};
#endif
