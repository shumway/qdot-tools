// $Id: ZincBlendeEMA.cc,v 1.7 2007/10/29 19:23:24 jshumwa Exp $
/*
    Copyright (C) 2004-7 John B. Shumway, Jr.

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
#include "ZincBlendeEMA.h"
#include "GridFactory.h"
#include "CompositionGrids.h"
#include "StrainGrid.h"
#include <complex>

extern "C" void zhpev_(const char *jobz, const char *uplo, const int *n,
  const std::complex<double>* ap, const double *w, 
  std::complex<double>* z, const int *ldz, std::complex<double> *work,
  const double *rwork, const int *info);

ZincBlendeEMA::ZincBlendeEMA(const GridFactory &factory, 
    const param& mat1, const param& mat2, const int ind1, const int ind2)
  : factory(factory), ind1(ind1), ind2(ind2), mat1(mat1), mat2(mat2),
    ve(*factory.getNewGrid<float>()), vh(*factory.getNewGrid<float>()) {
}

ZincBlendeEMA::~ZincBlendeEMA() {
  delete &ve,&vh;
}

void ZincBlendeEMA::calculate(const CompositionGrids& composition,
    const StrainGrid& strain) {
  const int N=6;
  // upper packed array: h_i,j=hamil(i+(j-1)*j/2)
  CArray1 hamil(21),work(2*N-1);
  Array1 en(N), rwork(3*N-2);
  CArray2 vec(N,N);
  const char jobz='V', uplo='U';
  const int XX=0,YY=1,ZZ=2,YZ=3,XZ=4,XY=5;
  const double sqrt3=sqrt(3.);
  IVec extent=ve.getExtent();
  for (int i=0; i<extent[0]; ++i ) {
    for (int j=0; j<extent[1]; ++j ) {
      for (int k=0; k<extent[2]; ++k ) {
        double trStrain=strain.getTraceGrid()(i,j,k);
        StrainGrid::SymMat3 eps=strain(i,j,k);
        double biaxStrain=2*eps(ZZ)-eps(XX)-eps(YY);
        float n1=composition.getGrid(ind1)(i,j,k);
        float n2=composition.getGrid(ind2)(i,j,k);
        double x=n2/(n1+n2+1e-6);
        const double oneThird=1.0/3.0;
        // Electron offset is linear in composition and strain.
        ve(i,j,k)=x*mat2.ec+(1.-x)*mat1.ec+trStrain*(
                  x*mat2.ac+(1.-x)*mat1.ac)+oneThird*biaxStrain*(
                  x*mat2.xi+(1.-x)*mat1.xi);
        // Diagonalize matrix to find hole offset.
        hamil=0.0;
        double b=x*mat2.b+(1-x)*mat1.b;
        double d=x*mat2.d+(1-x)*mat1.d;
        hamil(0)=hamil(2)=hamil(5)=trStrain*(x*mat2.av+(1-x)*mat1.av);
        hamil(0)-=b*(-2*eps(XX)+eps(YY)+eps(ZZ));
        hamil(1)=sqrt3*d*eps(XY);
        hamil(2)-=b*(eps(XX)-2*eps(YY)+eps(ZZ));
        hamil(3)=sqrt3*d*eps(XZ);
        hamil(4)=sqrt3*d*eps(YZ);
        // HACK: push pz energy up by 1000 Ha to remove LH contribution at 
        // edge of dot.
        hamil(5)-=b*(eps(XX)+eps(YY)-2*eps(ZZ))+1000;
        double SO=(x*mat2.so+(1.-x)*mat1.so)/3.0;
//        SO = 0.;
        std::complex<double> I=std::complex<double>(0.0,1.0);
        std::complex<double> zSO=I*SO;
        hamil(9)=hamil(0);
        hamil(13)=hamil(1);
        hamil(14)=hamil(2);
        hamil(18)=hamil(3);
        hamil(19)=hamil(4);
        hamil(20)=hamil(5);
        hamil(1)-=zSO;
        hamil(8)-=SO;
        hamil(12)+=zSO;
        hamil(13)+=zSO;
        hamil(15)+=SO;
        hamil(16)-=zSO;
        // Diagonalize the hamiltonian.
        int info=0;
        zhpev_(&jobz, &uplo, &N, hamil.data(), en.data(),
               vec.data(), &N, work.data(), rwork.data(), &info);
        if (info!=0) std::cout << "ERROR DIAGONALIZING IN ZHPEV" << std::endl;
        // HH and LH are states 3,4 and 5,6, but order can change.
        // Define heavy hole band to be the one with the largest lz * sz.
        double lz3=std::real(I*vec(3,0)*std::conj(vec(3,1)))
                  -std::real(I*vec(3,3)*std::conj(vec(3,4)));
        double lz5=std::real(I*vec(5,0)*std::conj(vec(5,1)))
                  -std::real(I*vec(5,3)*std::conj(vec(5,4)));
        int indexHH=(lz3 > lz5) ? 3 : 5;
        vh(i,j,k)=x*mat2.ev+(1-x)*mat1.ev+en(indexHH);
      }
    }
  }
}

/// Conversion eVtoHa.
const double ZincBlendeEMA::eVtoHa=1./27.211396;

// Data from Tables II, III and IV of Chris g. Van de Walle, "Band lineups and 
// deformation potentials in model-solid theory," PRB 39, 1871-1883, (1989).

// Si            (    ec ,  ev ,  ac , av ,   b ,  d  , so , xi ) in eV.
const ZincBlendeEMA::param 
ZincBlendeEMA::SI(  -5.85,-7.03, 4.18,2.46,-2.35,-5.32,0.04,9.16);

// Ge            (    ec ,  ev ,  ac , av ,   b ,  d  , so , xi ) in eV.
const ZincBlendeEMA::param 
ZincBlendeEMA::GE(  -5.51,-6.35, 2.55,1.24,-2.55,-5.50,0.30,9.42);

// GaAs            (  ec ,  ev ,  ac , av ,   b ,  d  , so , xi ) in eV.
const ZincBlendeEMA::param 
ZincBlendeEMA::GAAS(-5.29,-6.92,-7.17,1.16,-1.90,-4.23,0.34,0.00);

// InAs            (  ec ,  ev ,  ac , av ,   b ,  d  , so , xi ) in eV.
const ZincBlendeEMA::param
ZincBlendeEMA::INAS(-6.13,-6.67,-5.08,1.00,-1.55,-3.10,0.38,0.00);

// GaSb            (  ec ,  ev ,  ac , av ,   b ,  d  , so , xi ) in eV.
const ZincBlendeEMA::param
ZincBlendeEMA::GASB(-5.23,-6.25,-6.85,0.79,-2.0,-4.8,0.82,0.00);

