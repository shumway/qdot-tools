// $Id: VFF.cc,v 1.2 2004/08/14 01:42:07 jshumwa Exp $
/*
    Copyright (C) 2004 Matt Harowitz and John Shumway

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
#include "VFF.h"
#include <cmath>
#include <iostream>

VFF::VFF()
  :param(10) {
  //set parameters
  //assume that As, Ga, In, Sb are indexed 0,1,2,3,4 respectively
  //(GaAs,InAs,InGaAs,GaSb,GaAsSb,error)
//http://www.tbrri.com/~orubel/publications/2008/PhysRevB_v78_p075207_y2008.pdf GaSb
  cosphi[0] = -0.33333;	//-0.32557;
  cosphi[1] = -0.33333;	//-0.32557;
  cosphi[2] = -0.35021; //InGaAs
  cosphi[3] = -0.33333; //GaSb
  cosphi[4] = -0.33333; //GaAsSb
  cosphi[5] = 0.0; //error
  beta[0] = 6.1122E-3;//6.018E-3;
  beta[1] = 3.9129E-3;//3.700E-3;
  beta[2] = 5.0125E-3;//4.859E2;
  beta[3] = 4.7081E-3; //GaSb
  beta[4] = 5.4101E-3; //GaAsSb
  beta[5] = 0.0; //error

  // As-As: (shouldn't occur)
  param[0].d = 0;
  param[0].alpha1 = 0;
  // GaAs:
  param[1].d = 4.62595;
  param[1].alpha1 = 2.6311E-2;//2.0652E-2;
  // Ga-Ga: (shouldn't occur)
  param[2].d = 0;
  param[2].alpha1 = 0;
  // InAs:
  param[3].d = 4.97625;
  param[3].alpha1 = 2.2454E-2;//1.3921E-2;
  // In-Ga: (shouldn't occur)
  param[4].d = 0;
  param[4].alpha1 = 0;
  // In-In: (shouldn't occur)
  param[5].d = 0;
  param[5].alpha1 = 0;
  // As-Sb: (shouldn't occur)
  param[6].d = 0;
  param[6].alpha1 = 0;
  // GaSb
  param[7].d = 4.9870;
  param[7].alpha1 = 2.0425E-2;
  // In-Sb: (shouldn't occur)
  param[8].d = 0;
  param[8].alpha1 = 0;
  // Sb-Sb: (shouldn't occur)
  param[9].d = 0;
  param[9].alpha1 = 0;
}

double VFF::f2(const double r, const int i, const int j) const {
  const VFFParam& p1 = param[ind(i,j)];
  double y1 = p1.d * p1.d;
  double y2 = r * r - y1;
  return 0.375 * y2 * y2 * p1.alpha1 / y1;
}

double VFF::drF2(const double r, const int i, const int j) const {
  const VFFParam& p1 = param[ind(i,j)];
  double y1 = p1.d * p1.d;
  double y2 = r * r - y1;
  return 1.5 * r * y2 * p1.alpha1 / y1;


}

double VFF::h(const double r1, const double r2, const double costheta,
                   const int i, const int j, const int k ) const {
  int index1 = ind(i,j);
  int index2 = ind(i,k);
  const VFFParam& p1 = param[index1];
  const VFFParam& p2 = param[index2];
  int n;
  switch (index1 * index2){
    case 1: //Ga-As-Ga
    n = 0;
    break;
    case 9: //In-As-In
    n = 1;
    break;
    case 3: //Ga-As-In
    n = 2;
    break;
    case 14: //Ga-Sb-Ga
    n = 3;
    break;
    case 6: //As-Ga-Sb
    n = 4;
    break;
    default: //something else (error)
    n = 5;
    break;
  }
  double y1 = p1.d * p2.d;
  double y2 = r1 * r2 * costheta - y1 * cosphi[n];
  return 0.375 * y2 * y2 * beta[n] / y1;
}

void VFF::dh(const double r1, const double r2, const double costheta,
                  const int i, const int j, const int k,
                  double& dr1H, double& dr2H, double& dthetaH) const {
  int index1 = ind(i,j);
  int index2 = ind(i,k);
  const VFFParam& p1 = param[index1];
  const VFFParam& p2 = param[index2];
  int n;
  switch (index1 * index2){
    case 1: //Ga-As-Ga
    n = 0;
    break;
    case 9: //In-As-In
    n = 1;
    break;
    case 3: //Ga-As-In
    n = 2;
    break;
    case 14: //Ga-Sb-Ga
    n = 3;
    break;
    case 6: //As-Ga-Sb
    n = 4;
    break;
    default: //something else (error)
    n = 5;
    break;
  }
  double y1 = p1.d * p2.d;
  double y2 = r1 * r2 * costheta - y1 * cosphi[n];
  //double y3 = r1 * r1 - p1.d * p1.d;
  dr1H = 0.75 * beta[n] * y2 * r2 * costheta / y1;
  dr2H = 0.75 * r1 * costheta * beta[n] * y2 / y1;
  dthetaH = 0.75 * r1 * r2 * beta[n] * y2 / y1;  //really dcosthetaH
}
