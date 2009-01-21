// $Id: Minimize1D.cc,v 1.3 2006/05/29 15:48:45 jshumwa Exp $
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
#include "Minimize1D.h"
#include <cmath>
#include <cstdlib>
#include <iostream>

Minimize1D::Minimize1D(Function& f) : f(f) {
}

Minimize1D::~Minimize1D() {
}

void Minimize1D::bracketMin(double& ax, double& bx, double& cx,
                            double& fa, double& fb, double& fc) {
  if (fa<fb) { //Switch a and b if b not downhill.
    double temp;
    temp=ax; ax=bx; bx=temp;
    temp=fa; fa=fb; fb=temp;
//std::cout << "REVERSED" << std::endl;
  }
  cx=bx+(1+GOLD)*(bx-ax);
  fc=f(cx);
  while (fb>fc) { //Loop until we bracket.
    //Compute u by parabolic extapolation from a,b,c
    double r=(bx-ax)*(fb-fc);
    double q=(bx-cx)*(fb-fa);
    double u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0*(q-r>0)?1:-1*
		     (fabs(q-r)>1e-20)?fabs(q-r):1e-20);
    double ulim=bx+20*(cx-bx);
    double fu=0;
    if ( (bx-u)*(u-cx) > 0.0) {
      //Parabolic u is between b and c
      fu=f(u);
      if (fu<fc) {
        //Got a minimum between b and c.
        ax=bx; bx=u; fa=fb; fb=fu;
        return;
      } else if (fu>fb) {
        //Got a minimum between a and u.
        cx=u; fc=fu;
        return;
      }
      // Parabola is no use, use default maganification
      u=cx+(1+GOLD)*(cx-bx);
      fu=f(u);
    } else if ( (cx-u)*(u-ulim) > 0.0 ) {
      //Parabolic u is between c and limit
      fu=f(u);
      if (fu<fc) {
        bx=cx; cx=u; u=cx+(1+GOLD)*(cx-bx);
        fb=fc; fc=fu; fu=f(u);
      }
    } else if ( (u-ulim)*(ulim-cx) >= 0.0) {
      //Limit parabolic u to max value
      u=ulim; fu=f(u);
    } else {
      // Reject parabolic u, use default maganification
      u=cx+(1+GOLD)*(cx-bx);
      fu=f(u);
    }
    ax=bx; bx=cx; cx=u; // Eliminate oldest point and continue. 
    fa=fb; fb=fc; fc=fu; 
  }
}

double Minimize1D::goldenSearch(const double ax, const double bx, 
                                const double cx, double& xmin){
  double x0=ax;
  double x3=cx;
  double x1,x2;
  if (fabs(cx-bx) > fabs(bx-ax)) {
    x1=bx; 
    x2=bx+(1.0-GOLD)*(cx-bx);
  } else {
    x2=bx; 
    x1=bx-(1.0-GOLD)*(bx-ax);
  }
  double f1=f(x1);
  double f2=f(x2);
//std::cout << "x1=" << x1 << ", f1=" << f1 << std::endl;
//std::cout << "x2=" << x2 << ", f2=" << f2 << std::endl;
  while(fabs(x3-x0) > 3e-8*(fabs(x1)+fabs(x2))) {
    if (f1 < f2) {
      x3=x2; x2=x1; x1=GOLD*x2+(1-GOLD)*x0;
      f2=f1; f1=f(x1);
    } else {
      x0=x1; x1=x2; x2=GOLD*x1+(1-GOLD)*x3;
      f1=f2; f2=f(x2);
    }
  }
  if (f1 < f2) {
    xmin=x1;
    return f1;
  } else {
    xmin=x2;
    return f2;
  }
}

double Minimize1D::brentSearch(const double ax, const double bx, 
                               const double cx, double& xmin){
  const int ITMAX=100;
  const double ZEPS=1e-10;
  const double tol=1e-2;
  //const double tol=3e-8;
  double e=0.,d=0.;
  double a,b;
  if (ax<cx) {a=ax; b=cx;} else {a=cx; b=ax;}
  double x,w,v,fw,fv,fx;
  x=w=v=bx;
  fx=fw=fv=f(x);
  for (int iter=1; iter<=ITMAX; ++iter) {
    double xmid=0.5*(a+b);
    double tol1=tol*fabs(x)+ZEPS;
    double tol2=2.0*tol1;
    if (fabs(x-xmid) <= (tol2-0.5*(b-a))) {
std::cout << "niter=" << iter << std::endl;
      xmin=x;
      return fx;
    }
    if (fabs(e) > tol1) {
      //Construct a trial parabolic fit.
      double r=(x-w)*(fx-fv);
      double q=(x-v)*(fx-fw);
      double p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if (q>0.0) p=-p;
      q=fabs(q);
      double etemp=e;
      e=d;
      if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x)) {
        // Parabolic fit is unacceptable, use golden section step.
        //std::cout << "didn't use parabolic fit: " << iter << std::endl;
        d=(1.-GOLD)*(e=(x>=xmid ? a-x : b-x));
      } else {
        // Parabolic fit is acceptable.
        //std::cout << "used parabolic fit: " << iter << std::endl;
        d=p/q;
        double u=x+d;
        if (u-a < tol2 || b-u <tol2) d= (xmid>x?tol1:-tol1);
      }
    } else {
        //std::cout << "didn't try parabolic fit: " << iter << std::endl;
        d=(1.-GOLD)*(e=(x>=xmid ? a-x : b-x));
    }
    double u=( fabs(d)>=tol1 ? x+d : x+(d>0?tol1:-tol1) );
    double fu=f(u);
    if (fu <= fx) {
      if (u>=x) a=x; else b=x;
      v=w;w=x;x=u;
      fv=fw;fw=fx;fx=fu;
    } else {
      if (u<x) a=u; else b=u;
      if (fu<=fw || w==x) {
        v=w; w=u;
        fv=fw; fw=fu;
      } else if (fu<=fv || v==x || v==w) {
        v=u; fv=fu;
      }
    }
  }
  std::cout << "Warning: too many iterations in brentSearch" << std::endl;
  xmin=x;
  return fx;
}

const double Minimize1D::GOLD=0.618034;
