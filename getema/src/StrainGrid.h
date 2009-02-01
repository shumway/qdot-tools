// $Id: StrainGrid.h,v 1.2 2004/06/24 18:28:41 jshumwa Exp $
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
#ifndef __StrainGrid_h_
#define __StrainGrid_h_
class GridFactory;
class ElasticConstants;
class StressGrid;
#include <blitz/tinyvec.h>
#include <blitz/tinymat.h>
#include "Grid.h"

class StrainGrid {
public:
  /// Typedefs.
  typedef blitz::TinyVector<int,3> IVec;
  typedef blitz::TinyVector<float,6> SymMat3;
  typedef Grid<SymMat3> SMGrid;
  static double trace(const SymMat3& m) {return m[0]+m[1]+m[2];}
  static double biaxial(const SymMat3& m) {return 2*m[2]-m[0]-m[1];}
  /// Constructor.
  StrainGrid() {}
  /// Virtual destructor.
  virtual ~StrainGrid() {}
  /// Get reference to the grid.
  virtual SMGrid& getGrid()=0;
  /// Get const reference to the grid.
  virtual const SMGrid& getGrid() const=0;
  /// Get reference to the trace grid.
  virtual Grid<float>& getTraceGrid()=0;
  /// Get const reference to the trace grid.
  virtual const Grid<float>& getTraceGrid() const=0;
  /// Get reference to the trace grid.
  virtual Grid<float>& getBiaxialGrid()=0;
  /// Get const reference to the trace grid.
  virtual const Grid<float>& getBiaxialGrid() const=0;
  /// Get reference to a grid element.
  virtual SymMat3& operator()(const int, const int, const int)=0;
  /// Get const reference to a grid element.
  virtual const SymMat3& operator()(const int, const int, const int) const=0;
};
#endif
