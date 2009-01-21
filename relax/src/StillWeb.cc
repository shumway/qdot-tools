//$Id: StillWeb.cc,v 1.5 2005/01/21 17:29:27 jshumwa Exp $
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
#include "StillWeb.h"
#include <cmath>
#include <iostream>

StillWeb::StillWeb()
  : A(7.049556277), B(0.6022245584), p(4), q(0), a(1.80), gamma(1.20), param(3) {
  // Si-Si:
  param[0].lambda=21.0;
  param[0].sigmainv=1./((2.347772155/0.529177)/pow(2,1./6));
  param[0].epsilon=(2.17/27.211);
  // Si-Ge, Ge-Si:
  param[1].lambda=sqrt(21.0*31.0);
  param[1].sigmainv=1./((2.396885/0.529177)/pow(2,1./6));
  param[1].epsilon=(2.0427/27.211);
  // Ge-Ge
  param[2].lambda=31.0;
  param[2].sigmainv=1./((2.44598/0.529177)/pow(2,1./6));
  param[2].epsilon=(1.93/27.211);
}

double StillWeb::f2(const double r, const int i, const int j) const {
  const SWParam& p1=param[ind(i,j)];
  double y=r*p1.sigmainv;
  //return (y<a) ? p1.epsilon*A*(B*pow(y,-p)-1)*exp(1/(y-a)) : 0;
  return (y<a) ? p1.epsilon*(A*(B/(y*y*y*y)-1)*exp(1/(y-a))+1) : p1.epsilon;
}

double StillWeb::drF2(const double r, const int i, const int j) const {
  const SWParam& p1=param[ind(i,j)];
  double y=r*p1.sigmainv;
  double y4=y*y; y4*=y4;
  return (y<a) ? -p1.epsilon*p1.sigmainv*A*exp(1/(y-a))
                //*((B*pow(y,-p)-1)*pow((y-a),-2)+B*p*pow(y,-p-1)) : 0;
                *((B/y4-1)*1./((y-a)*(y-a))+B*p/(y4*y)) : 0;
}

double StillWeb::h(const double r1, const double r2, const double costheta,
                   const int i, const int j, const int k ) const {
  const SWParam& p1=param[ind(i,j)];
  const SWParam& p2=param[ind(i,k)];
  const double onethird=1./3.;
  double y1=r1*p1.sigmainv;
  double y2=r2*p2.sigmainv;
  return (y1<a && y2<a)
	  ? sqrt(p1.epsilon*p2.epsilon*p1.lambda*p2.lambda)
//            *exp(gamma/(y1-a)+gamma/(y2-a))*pow((costheta+onethird),2)
            *exp(gamma/(y1-a)+gamma/(y2-a))
	    *(costheta+onethird)*(costheta+onethird)
	  : 0;
}

void StillWeb::dh(const double r1, const double r2, const double costheta,
                  const int i, const int j, const int k,
                  double& dr1H, double& dr2H, double& dthetaH) const {
  const SWParam& p1=param[ind(i,j)];
  const SWParam& p2=param[ind(i,k)];
  dr1H=dr2H=dthetaH=0;
  double y1=r1*p1.sigmainv;
  double y2=r2*p2.sigmainv;
  if (y1<a && y2<a) {
     double hr=sqrt(p1.epsilon*p2.epsilon*p1.lambda*p2.lambda)
              *exp(gamma/(y1-a)+gamma/(y2-a));
     //double h=hr*pow((costheta+1./3),2);
     double h=hr*(costheta+1./3)*(costheta+1./3);
     //dr1H=-gamma*h*pow(y1-a,-2)*p1.sigmainv;
     //dr2H=-gamma*h*pow(y2-a,-2)*p2.sigmainv;
     dr1H=-gamma*h*p1.sigmainv/((y1-a)*(y1-a));
     dr2H=-gamma*h*p2.sigmainv/((y2-a)*(y2-a));
     dthetaH=2*hr*(costheta+1./3);
  } else {
    dr1H=dr2H=dthetaH=0;
  }
}

double StillWeb::bind(const int i, const int j) const {
  return -param[ind(i,j)].epsilon;
}
