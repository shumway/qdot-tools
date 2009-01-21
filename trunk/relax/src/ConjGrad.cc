// $Id: ConjGrad.cc,v 1.11 2007/12/28 02:19:40 jshumwa Exp $
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "ConjGrad.h"
#include "Forces.h"
#include "TotalEnergy.h"
#include "CoordsH5.h"
#include "Checkpoint.h"
#include "CellRelaxer.h"
#include <iostream>

#ifdef HAVE_BLAS
#define DAXPY_F77 F77_FUNC(daxpy,DAXPY)
extern "C" void DAXPY_F77(const int *N, const double *ALPHA, const double *X, 
                       const int *INCX, double *Y, const int *INCY);
#define SAXPY_F77 F77_FUNC(saxpy,SAXPY)
extern "C" void SAXPY_F77(const int *N, const float *ALPHA, const float *X, 
                       const int *INCX, float *Y, const int *INCY);
#endif

ConjGrad::ConjGrad(Forces& f, TotalEnergy& en, CoordsH5& coords, 
    const Checkpoint* checkpt, const double ftol, const int itermax,
    CellRelaxer* cellRelaxer, const int cellRelaxSteps, const double scale)
  : minimizer(*this),
    f(f), en(en), coords(coords), direction(f.getForces().size()),
    checkpt(checkpt), ftol(ftol), itermax(itermax),
    cellRelaxer(cellRelaxer), cellRelaxSteps(cellRelaxSteps),
    scale(scale) {
}

ConjGrad::~ConjGrad() {
  delete cellRelaxer;
}

void ConjGrad::minimize() {
  // Initialization.
  std::cout << "*** Initializing conjugant gradients ***" << std::endl;
  std::cout << "ftol=" << ftol << ", maxiter=" << itermax << std::endl;
  int natom=direction.size();
  en.compute(); en.report();
  f.compute();
  double fold2=0, fnew2=0;
  // Compute maximum force on any atom.
  double fmax = 0;
  for (int iatom = 0; iatom<natom; ++iatom) {
    double fi=dot(f[iatom],f[iatom]);
    fold2+=fi;
    if (fi>fmax) fmax=fi;
  }
  fmax=sqrt(fmax);
  std::cout << "Maximum force on any atom is " << fmax << std::endl;
  direction=f.getForces();
  // CG loop.
  for (int iter=0;iter<itermax;++iter) {
    std::cout << "*** CG iteration " << iter << " ***"<< std::endl;
    // Relax the supercell.
    if (cellRelaxer && iter%cellRelaxSteps==0) cellRelaxer->relax();
    // Line minimize.
    lineMinimize();
    f.compute();
    // Compute maximum force on any atom.
    fmax=0;
    for (int iatom = 0; iatom<natom; ++iatom) {
      double fi=dot(f[iatom],f[iatom]);
      if (fi>fmax) fmax=fi;
    }
    fmax=sqrt(fmax);
    std::cout << "Maximum force on any atom is " << fmax << std::endl;
    if (fmax<ftol) break;
    //if (fmax<2e-6) break;
    // Compute new conjugate direction.
    fnew2=0;
    //double foldFnew=0;
    for (int iatom=0; iatom<natom;++iatom) {
      fnew2+=dot(f[iatom],f[iatom]);
      //foldFnew+=(*fnew)[iatom].dot((*fold)[iatom]);
    }
    double gamma = fnew2/fold2;
    //double gamma = (fnew2-foldFnew)/fold2;
    for (int iatom=0; iatom<natom;++iatom) {
      (direction[iatom]*=gamma)+=f[iatom];
      //direction[iatom]=f[iatom];
      //direction[iatom]=(*fnew)[iatom];
    }
    fold2=fnew2;
    if (checkpt) checkpt->write(iter);
  }
  // Finished
  std::cout << "*** Finished conjugant gradients ***" << std::endl;
  if (fmax>ftol) {
    std::cout << "Warning: Failed to reach force tolerance " << ftol
              << " in " << itermax << " steps." << std::endl;
  }
  en.compute(); en.report();
}

double ConjGrad::operator()(const double x) {
  int natom=direction.size();
  if (x!=xLast) {
#ifdef HAVE_BLAS
    // BLAS call for coord+=(x-xLast)*direction.
    const int N=3*natom, ONE=1;
#ifdef ENABLE_FLOAT
    float alpha=x-xLast;
    SAXPY_F77(&N,&alpha,(const float*)&direction[0],&ONE,
                              (float*)&coords[0],&ONE);
#else
    double alpha=x-xLast;
    DAXPY_F77(&N,&alpha,(const double*)&direction[0],&ONE,
                              (double*)&coords[0],&ONE);
#endif
#else
    for (int iatom=0; iatom<natom;++iatom) {
      coords[iatom]+=direction[iatom]*(x-xLast);
    }
#endif
    xLast=x;
    en.compute();
  }
  return en.getEnergy();
}

void ConjGrad::lineMinimize() {
  // Start with points ax and bx.
  double ax=xLast=0;
  double ena=(*this)(ax);
  //double bx=10.0;
  double bx=2.0*scale;
  double enb=(*this)(bx);
  std::cout << "ax="<<ax<<", bx="<<bx<< std::endl;
  std::cout << "ena="<<27.211*ena<<", enb="<<27.211*enb<< std::endl;
  // Bracket the minimum.
  double cx=0,enc=0;
  std::cout << "Bracketing minimum." << std::endl;
  minimizer.bracketMin(ax,bx,cx,ena,enb,enc);
  std::cout << "ax="<<ax<<", bx="<<bx<<", cx="<<cx << std::endl;
  std::cout << "ena="<<27.211*ena<<", enb="<<27.211*enb<<", enc="<<27.211*enc << std::endl;
  // Find solution using golden section search.
  double xmin, emin;
  std::cout << "Searching for minimum." << std::endl;
  //emin = minimizer.goldenSearch(ax,bx,cx,xmin);
  emin = minimizer.brentSearch(ax,bx,cx,xmin);
  (*this)(xmin);
  en.compute();
  en.report();
}
