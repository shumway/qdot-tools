// $Id: StrainFromStress.cc,v 1.1.1.1 2004/05/03 16:49:21 jshumwa Exp $
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
#include "StrainFromStress.h"
#include "GridFactory.h"
#include "ElasticConstants.h"
#include "StressGrid.h"

extern "C" void sgetrf_(const int*, const int*, float*, const int*, const int*, int*);
extern "C" void sgetri_(const int*, float*, const int*, const int*, float*, const int*, int*); 

StrainFromStress::StrainFromStress(const GridFactory &factory,
    const ElasticConstants &elasticity, const StressGrid& stress) 
  : grid(factory.getNewGrid<SymMat3>()),
    traceGrid(factory.getNewGrid<float>()), 
    biaxialGrid(factory.getNewGrid<float>()) {
  ElasticConstants::Mat mat;
  static const int SIX=6;
  blitz::TinyVector<int,SIX> ipiv;
  static const int lwork=SIX*SIX;
  blitz::Array<float,1> work(lwork);
  int info;
  IVec extent=grid->getExtent();
  for (int i=0; i<extent[0]; ++i ) {
    for (int j=0; j<extent[1]; ++j ) {
      for (int k=0; k<extent[2]; ++k ) {
        mat=elasticity(i,j,k);
        // Find inverse matrix with LU decomposition.
        sgetrf_(&SIX,&SIX,mat.data(),&SIX,ipiv.data(),&info);
        if (info!=0) std::cout << "BAD RETURN FROM SGETRF!!!!" << std::endl;
        sgetri_(&SIX,mat.data(),&SIX,ipiv.data(),work.data(),&lwork,&info);
        if (info!=0) std::cout << "BAD RETURN FROM SGETRI!!!!" << std::endl;
        // Compute strain tensor using elastic constants.
        (*grid)(i,j,k)=blitz::product(mat,stress(i,j,k));
        (*traceGrid)(i,j,k)=trace((*grid)(i,j,k));
        (*biaxialGrid)(i,j,k)=biaxial((*grid)(i,j,k));
      }
    }
  }
}

StrainFromStress::~StrainFromStress() {
  delete grid, traceGrid, biaxialGrid;
}
