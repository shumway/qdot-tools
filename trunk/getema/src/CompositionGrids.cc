// $Id: CompositionGrids.cc,v 1.1.1.1 2004/05/03 16:49:21 jshumwa Exp $
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
#include "CompositionGrids.h"
#include "SpeciesH5.h"
#include "CoordsH5.h"
#include "GridFactory.h"

CompositionGrids::CompositionGrids(const GridFactory& factory,
    const SpeciesH5& species, const CoordsH5& coords)
  : grid(species.getNSpecies()), species(species), coords(coords) {
  std::cout << "Setting composition grids for " << grid.size()
            << " species." << std::endl;
  for (int i=0; i<grid.size(); ++i) grid[i]=factory.getNewGrid<float>();
} 

CompositionGrids::~CompositionGrids() {
  for (int i=0; i<grid.size(); ++i) delete grid[i];
}

void CompositionGrids::calculate() {
  for (int i=0; i<species.getNPart(); ++i) {
    grid[species(i)-1]->addData(coords[i],1.);
  }
}
