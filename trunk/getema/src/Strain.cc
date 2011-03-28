// $Id: Strain.cc,v 1.3 2007/12/06 23:59:59 jshumwa Exp $
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
#include "Strain.h"
#include "GridFactory.h"
#include "CoordsH5.h"
#include "SuperCellH5.h"
#include "CompositionGrids.h"
#include "NeighborsH5.h"

extern "C" void sgetrf_(const int*, const int*, float*, const int*, const int*, int*);
extern "C" void sgetri_(const int*, float*, const int*, const int*, float*, const int*, int*); 

Strain::Strain(const GridFactory &factory, const CoordsH5 &coords,
    const CompositionGrids &composition, const NeighborsH5 &neighbors,
    const SuperCellH5 &cell, bool isInGaAs) 
  : grid(factory.getNewGrid<SymMat3>()),
    traceGrid(factory.getNewGrid<float>()), 
    biaxialGrid(factory.getNewGrid<float>()),
    coords(coords), neighbors(neighbors) {
  // Allocate work grids to store vector displacements.
  //Grid<Vec> *deltax = factory.getNewGrid<Vec>();
  //Grid<Vec> *deltay = factory.getNewGrid<Vec>();
  //Grid<Vec> *deltaz = factory.getNewGrid<Vec>();
  int natom=coords.getNAtom();
  std::cout << "Calculating stain for " << natom << " atoms." << std::endl;
  // First build x, y, and z vectors from each atom and store in each cell.
  for (int i=0; i<natom; ++i) {
    Vec pos = coords[i];
    Vec x=Vec(0.,0.,0.);
    Vec y=Vec(0.,0.,0.);
    Vec z=Vec(0.,0.,0.);
    for (int jnbr=0; jnbr<neighbors.nmax; ++jnbr) {
      int j = neighbors.getNeighbor(i,jnbr)-1; 
      if (j==-1) break;
      Vec delta = coords[j] - pos;
      cell.pbc(delta);
      x += (delta[0]>0) ? delta : -delta;
      y += (delta[1]>0) ? delta : -delta;
      z += (delta[2]>0) ? delta : -delta;
    }
    SymMat3 symDelta(0.,0.,0.,0.,0.,0.);
    symDelta *= 0;
    symDelta(0) = x[0];
    symDelta(1) = y[1];
    symDelta(2) = z[2];
    symDelta(3) = 0.5*(y[2]+z[1]);
    symDelta(4) = 0.5*(z[0]+x[2]);
    symDelta(5) = 0.5*(x[1]+y[0]);
    grid->addData(pos,symDelta);
  }
  // Next calculate strain tensor in each cell.
  std::cout << "Normalizing strain tensors" << std::endl;

  const Grid<float> &numAs(composition.getGrid(0));
  const Grid<float> &numGa(composition.getGrid(1));
  const Grid<float> &numIn(composition.getGrid(2));
  const Grid<float> &numSb(composition.getGrid(3));

  SymMat3 one(1.,1.,1.,0.,0.,0.);
  IVec extent=grid->getExtent();
  const double aGaAs = 10.683;
  const double aInAs = 11.449;
  const double aGaSb = 11.5196; //http://www.ioffe.ru/SVA/NSM/Semicond/GaSb/basic.html
  for (int i=0; i<extent[0]; ++i ) {
    for (int j=0; j<extent[1]; ++j ) {
      for (int k=0; k<extent[2]; ++k ) {
        float nAs=numAs(i,j,k);
        float nGa=numGa(i,j,k);
        float nIn=numIn(i,j,k);
	float nSb=numSb(i,j,k);
        float ntot = nAs + nGa + nIn + nSb;
        float xIn = nIn / (nGa + nIn);
        float xSb = nSb / (nAs + nSb + 1E-6);
	double a=0;	

	if (isInGaAs) {
         a = xIn*aInAs + (1.-xIn)*aGaAs;
	} else { 
         a = xSb*aGaSb + (1.-xSb)*aGaAs;
	}

//  std::cout << a << std::endl;

        blitz::TinyVector<float,6> strain = (*grid)(i,j,k);
        strain = (strain/(a*ntot)) - one;
        (*grid)(i,j,k) = strain;
        (*traceGrid)(i,j,k) = trace(strain);
        (*biaxialGrid)(i,j,k) = biaxial(strain);
      }
    }
  }
  std::cout << "Done calculating strain tensors" << std::endl;
}

Strain::~Strain() {
  delete grid, traceGrid, biaxialGrid;
}
