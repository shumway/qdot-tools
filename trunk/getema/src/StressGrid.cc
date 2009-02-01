// $Id: StressGrid.cc,v 1.3 2006/07/01 23:07:17 jshumwa Exp $
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
#include "StressGrid.h"
#include "GridFactory.h"
#include "Stress.h"
#include "CoordsH5.h"
#include <blitz/matrix.h>

StressGrid::StressGrid(const GridFactory &factory, 
    const Stress &stress, const CoordsH5 &coordsH5) 
  : grid(factory.getNewGrid<SymMat3>()),
    traceGrid(factory.getNewGrid<float>()), 
    biaxialGrid(factory.getNewGrid<float>()), 
    normGrid(factory.getNewGrid<float>()), 
    stress(stress), coordsH5(coordsH5) {
}

StressGrid::~StressGrid() {
  delete grid;
  delete traceGrid;
}

void StressGrid::calculate() {
  grid->getData()=0.0; traceGrid->getData()=0.0;
  biaxialGrid->getData()=0.0; normGrid->getData()=0.0;
  for (int i=0; i<stress.getNPart(); ++i) {
    const SymMat3& m=stress(i);
    grid->addData(coordsH5[i],m);
    traceGrid->addData(coordsH5[i],trace(m));
    biaxialGrid->addData(coordsH5[i],biaxial(m));
    normGrid->addData(coordsH5[i],1);
  }
  grid->getData()/=normGrid->getData();
  traceGrid->getData()/=normGrid->getData();
  biaxialGrid->getData()/=normGrid->getData();
}
