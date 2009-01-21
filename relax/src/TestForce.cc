// $Id: TestForce.cc,v 1.3 2004/08/14 13:53:21 jshumwa Exp $
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
#include "TestForce.h"
#include "Forces.h"
#include "TotalEnergy.h"
#include "CoordsH5.h"
#include <iostream>
#include <blitz/tinyvec-et.h>

TestForce::TestForce(Forces& f, TotalEnergy& en, CoordsH5& coords)
  : f(f), en(en), coords(coords) {
}

TestForce::~TestForce() {
}

void TestForce::test() {
  std::cout << "Testing forces" << std::endl;
  f.compute();
  int natom=coords.coords.size(); 
  en.compute();
  for (int iatom=0; iatom<20; ++iatom) {
    Vec3 deltae;
    const double delta=3e-4;
    Vec3 orig=coords[iatom];
    // Compute force in x direction.
    coords[iatom][0]+=delta;
    en.compute(); deltae[0]=en.getEnergy();
    (coords[iatom]=orig)[0]-=delta;
    en.compute(); deltae[0]-=en.getEnergy();
    // Compute force in y direction.
    (coords[iatom]=orig)[1]+=delta;
    en.compute(); deltae[1]=en.getEnergy();
    (coords[iatom]=orig)[1]-=delta;
    en.compute(); deltae[1]-=en.getEnergy();
    // Compute force in z direction.
    (coords[iatom]=orig)[2]+=delta;
    en.compute(); deltae[2]=en.getEnergy();
    (coords[iatom]=orig)[2]-=delta;
    en.compute(); deltae[2]-=en.getEnergy();
    Vec3 deltaf=deltae; (deltaf*=(-0.5/delta))-=f[iatom];
    if (dot(deltaf,deltaf)>1e-15) {
      std::cout << "atom:" << iatom << std::endl;
      std::cout << "  force:" << f[iatom] << std::endl;
      std::cout << "  -dedr:" << deltae*(-0.5/delta) << std::endl;
    }
    coords[iatom]=orig;
  }
}
/*
void TestForce::minimize() {
  // Initialization.
  std::cout << "*** Initializing conjugant gradients ***" << std::endl;
  int natom=direction.size();
  Forces *fold(&f1), *fnew(&f2);
  en.compute(); en.report();
  fold->compute();
  double fold2=0, fnew2=0;
  for (int iatom=0; iatom<natom;++iatom) fold2+=(*fold)[iatom].length2();
  direction=fold->getForces();
  // CG loop.
  for (int iter=0;;++iter) {
    std::cout << "CG iteration " << iter << std::endl;
    // Line minimize.
    lineMinimize();
    fnew->compute();
    // Compute maximum force on any atom.
    double fmax = 0;
    for (int iatom = 0; iatom<natom; ++iatom) {
      double f=(*fnew)[iatom].length2();
      if (f>fmax) fmax=f;
    }
    fmax=sqrt(fmax);
    std::cout << "Maximum force on any atom is " << fmax << std::endl;
    if (fmax<0.003) break;
    // Compute new conjugate direction.
    fnew2=0;
    double foldFnew=0;
    for (int iatom=0; iatom<natom;++iatom) {
      fnew2+=(*fnew)[iatom].length2();
      foldFnew+=(*fnew)[iatom].dot((*fold)[iatom]);
    }
    double lambda = (fnew2-foldFnew)/fold2;
    for (int iatom=0; iatom<natom;++iatom) {
      //(direction[iatom]*=lambda)+=(*fnew)[iatom];
      direction[iatom]=(*fnew)[iatom];
    }
    // Move fnew to fold.
    if (fnew=&f1) {
      fnew=&f2; fold=&f1; fold2=fnew2;
    } else {
      fnew=&f1; fold=&f2; fold2=fnew2;
    }
  }
  // Finished
  std::cout << "*** Finished conjugant gradients ***" << std::endl;
  en.compute(); en.report();
}

double TestForce::operator()(const double x) {
  int natom=direction.size();
  if (x!=xLast) {
    for (int iatom=0; iatom<natom;++iatom) {
      coords[iatom]+=direction[iatom]*(x-xLast);
    }
    xLast=x;
  }
  en.compute();
  return en.getEnergy();
}

void TestForce::lineMinimize() {
  // Start with points ax and bx.
  double ax=xLast=0;
  double ena=(*this)(ax);
  en.report();
  double bx=0.02;
  double enb=(*this)(bx);
  en.report();
  // Bracket the minimum.
  double cx=0,enc=0;
  minimizer.bracketMin(ax,bx,cx,ena,enb,enc);
  std::cout << "ax=" << ax << ", bx=" << bx << ", cx=" << cx << std::endl;
  std::cout << std::scientific;
  std::cout << "ena=" << ena << ", enb=" << enb << ", enc=" << enc << std::endl;
  // Find solution using golden section search.
  double xmin, emin;
  emin = minimizer.goldenSearch(ax,bx,cx,xmin);
  std::cout << "Minimum energy " << emin << " at " << xmin << std::endl;
  (*this)(xmin);
  en.compute();
  en.report();
  //int natom=direction.size();
  //for (int iatom=0; iatom<natom;++iatom) coords[iatom]+=direction[iatom]*0.2;
  en.compute();
//  std::cout << "Energy is " << en.getEnergy()*27.211 << std::endl;
  double a=0, ae=en.getEnergy(); 
  double b=0.1;
  for (int iatom=0; iatom<natom;++iatom) coords[iatom]+=direction[iatom]*(b-a);
  en.compute();
  double be=en.getEnergy(); 
  double c=-0.1;
  for (int iatom=0; iatom<natom;++iatom) coords[iatom]+=direction[iatom]*(c-b);
  en.compute();
  double ce=en.getEnergy(); 
//  std::cout << std::scientific << "a= " << a << " ae=" << ae;
//  std::cout << " b= " << b << " be=" << be;
//  std::cout << " c= " << c << " ce=" << ce << std::endl;
  double xold(c),xeold(ce),x,xe;
  // loop 
  while (true) {
//  std::cout << std::scientific << "a= " << a << " ae=" << ae;
//  std::cout << " b= " << b << " be=" << be;
//  std::cout << " c= " << c << " ce=" << ce << std::endl;
    // extrapolate x from parabola abc
    x=b-0.5*((b-a)*(b-a)*(be-ce)-(b-c)*(b-c)*(be-ae))
         /((b-a)*(be-ce)-(b-c)*(be-ae));
    for (int iatom=0; iatom<natom;++iatom) {
      coords[iatom]+=direction[iatom]*(x-xold);
    }
    en.compute();
    xe=en.getEnergy();
    std::cout << std::scientific << "x= " << x << " xe=" << xe << std::endl;;
    if (fabs(xe-xeold)<1e-3) break;
    c=b;b=a;a=x; ce=be;be=ae;ae=xe;
    xold=x; xeold=xe;
  }
//  std::cout.unsetf(std::ios_base::scientific);
//  std::cout << "Energy is " << en.getEnergy()*27.211 << std::endl;
  std::cout.unsetf(std::ios_base::scientific);
}
*/
