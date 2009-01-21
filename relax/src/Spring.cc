// $Id: Spring.cc,v 1.3 2006/05/29 15:48:45 jshumwa Exp $
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
#include "Spring.h"
#include <cmath>
#include <iostream>

Spring::Spring() {}

double Spring::f2(const double r, const int i, const int j) const {
  return (r<rcut) ? 0.5*k*pow(r-r0,2):0;
}

double Spring::drF2(const double r, const int i, const int j) const {
  return  (r<rcut) ? k*(r-r0):0;
}

double Spring::h(const double r1, const double r2, const double costheta,
                   const int i, const int j, const int k ) const {
  return (r1<rcut && r2<rcut) ? pow(costheta+1./3.,2) : 0;
}

void Spring::dh(const double r1, const double r2, const double costheta,
                  const int i, const int j, const int k,
                  double& dr1H, double& dr2H, double& dthetaH) const {
  dr1H=dr2H=0;
  if (r1<rcut && r2<rcut) {
    dthetaH=2*(costheta+1./3.);
  }
}

const double Spring::rcut=5;
const double Spring::k=1;
const double Spring::r0=4.4366;
