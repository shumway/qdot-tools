// $Id: ZincBlendeElasticity.cc,v 1.3 2006/07/02 06:27:58 jshumwa Exp $
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
#include "ZincBlendeElasticity.h"
#include "GridFactory.h"
#include "CompositionGrids.h"

ZincBlendeElasticity::ZincBlendeElasticity(const GridFactory &factory, 
    const coef& mat1, const coef& mat2,
    const int ind1, const int ind2)
  : factory(factory), mat1(mat2), mat2(mat2), ind1(ind1), ind2(ind2),
    c11(factory.getNewGrid<float>()), c12(factory.getNewGrid<float>()),
    c44(factory.getNewGrid<float>()) {
}

ZincBlendeElasticity::~ZincBlendeElasticity() {
  delete c11; delete c12; delete c44;
}

void ZincBlendeElasticity::calculate(const CompositionGrids& composition) {
  // Use linear interpolation of elastic constants (Vegard's law).
  c11->getData()=(mat1.c11*composition.getGrid(ind1).getData()
                 +mat2.c11*composition.getGrid(ind2).getData())
                /(composition.getGrid(ind1).getData()
                 +composition.getGrid(ind2).getData())*tenMbar;
  c12->getData()=(mat1.c12*composition.getGrid(ind1).getData()
                 +mat2.c12*composition.getGrid(ind2).getData())
                /(composition.getGrid(ind1).getData()
                 +composition.getGrid(ind2).getData())*tenMbar;
  c44->getData()=(mat1.c44*composition.getGrid(ind1).getData()
                 +mat2.c44*composition.getGrid(ind2).getData())
                /(composition.getGrid(ind1).getData()
                 +composition.getGrid(ind2).getData())*tenMbar;
}

const ElasticConstants::Mat& ZincBlendeElasticity::operator()(
    const int i, const int j, const int k) const {
  mat=0;
  mat(0,0)=mat(1,1)=mat(2,2)=c11->getData()(i,j,k);
  mat(0,1)=mat(1,0)=mat(0,2)=mat(2,0)=mat(1,2)=mat(2,1)=c12->getData()(i,j,k);
  mat(3,3)=mat(4,4)=mat(5,5)=c44->getData()(i,j,k);
  return mat;
}



// Elastic constants in 10^12 dyne/cm2 = 1 Mbar (300K values).
// From "Semiconductors - Basic Data", O. Madelung (ed) (Springer, 1996).
//
// GaAs (c11,c12,c44)
const ZincBlendeElasticity::coef
ZincBlendeElasticity::GAAS(1.19,0.538,0.595);
// InAs (c11,c12,c44)
const ZincBlendeElasticity::coef
ZincBlendeElasticity::INAS(0.8329,0.4526,0.3959);
// GaSb (c11,c12,c44)
//http://www.ioffe.ru/SVA/NSM/Semicond/GaSb/mechanic.html
const ZincBlendeElasticity::coef
ZincBlendeElasticity::GASB(0.883,0.402,0.432);
// Si (c11,c12,c44)
const ZincBlendeElasticity::coef
ZincBlendeElasticity::SI(1.675,0.650,0.801);
// Ge (c11,c12,c44)
const ZincBlendeElasticity::coef
ZincBlendeElasticity::GE(1.315,0.494,0.684);

// Conversion from 1 Mbar to Ha/a0^3 (1bar=10^5 J/m^3).
const double
ZincBlendeElasticity::tenMbar=1e11*1e-30*(.529*.529*.529)/(27.211*1.6e-19);
