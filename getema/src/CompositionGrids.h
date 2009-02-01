// $Id: CompositionGrids.h,v 1.2 2006/07/02 06:27:58 jshumwa Exp $
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
#ifndef __CompositionGrids_h_
#define __CompositionGrids_h_
#include "Grid.h"
class GridFactory;
class SpeciesH5;
class CoordsH5;
#include <vector>

/// Class for projecting and storing the composition on a grid.
/// Grid point has number of atoms of that type in the cell.
/// @author John Shumway
/// @version $Revision: 1.2 $
class CompositionGrids {
public:
  /// Constructor.
  CompositionGrids(const GridFactory&, const SpeciesH5&, const CoordsH5&);
  /// Destructor.
  ~CompositionGrids();
  /// Get the number of composition grids.
  int getNGrid() const {return grid.size();}
  /// Get grid by index.
  Grid<float>& getGrid(const int i) const {return *(grid[i]);}
  /// Caclulate gridded compositon.
  void calculate();
private:
  /// The grids.
  std::vector<Grid<float>*> grid;
  /// The species data.
  const SpeciesH5& species;
  /// The coordinate data.
  const CoordsH5& coords;
};
#endif
