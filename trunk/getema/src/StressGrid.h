// $Id: StressGrid.h,v 1.3 2004/06/24 18:28:41 jshumwa Exp $
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
#ifndef __StressGrid_h_
#define __StressGrid_h_

class GridFactory;
class Stress;
class CoordsH5;
#include "Grid.h"
#include <vector>

/** Class for gridding stress.
 @bug Assumes grid is eight-atom cublic cell.
 @todo Need to store volume with atomic stress data, then use volume ratios.
 @author John Shumway
 @version $Revision: 1.3 $ */
class StressGrid {
public:
  /// Typedefs.
  typedef blitz::TinyVector<int,3> IVec;
  typedef blitz::TinyVector<float,6> SymMat3;
  typedef Grid<SymMat3> SMGrid;
  static double trace(const SymMat3& m) {return m[0]+m[1]+m[2];}
  static double biaxial(const SymMat3& m) {return 2*m[2]-m[0]-m[1];}
  /// Constructor.
  StressGrid(const GridFactory&, const Stress&, const CoordsH5&);
  /// Destructor.
  ~StressGrid();
  /// Get reference to the grid.
  SMGrid& getGrid() {return *grid;}
  /// Get const reference to the grid.
  const SMGrid& getGrid() const {return *grid;}
  /// Get reference to the trace grid.
  Grid<float>& getTraceGrid() {return *traceGrid;}
  /// Get const reference to the trace grid.
  const Grid<float>& getTraceGrid() const {return *traceGrid;}
  /// Get reference to the trace grid.
  Grid<float>& getBiaxialGrid() {return *biaxialGrid;}
  /// Get const reference to the trace grid.
  const Grid<float>& getBiaxialGrid() const {return *biaxialGrid;}
  /// Calculate gridded stress.
  void calculate();
  /// Get the stress tensor at a grid point.
  const SymMat3 operator()(const int i, const int j, const int k) const { 
    return (*grid)(i,j,k);}
private:
  /// The grid.
  SMGrid* grid;
  /// The grid for the trace of stress.
  Grid<float>* traceGrid;
  /// The grid for the biaxial stress.
  Grid<float>* biaxialGrid;
  /// The grid for normalization in average of stress.
  Grid<float>* normGrid;
  /// The species data.
  const Stress& stress;
  /// The coordinate data.
  const CoordsH5& coordsH5;
};
#endif
