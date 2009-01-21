// $Id: Mat3.h,v 1.3 2007/03/14 19:47:50 jshumwa Exp $
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
#ifndef __Mat3_h_
#define __Mat3_h_
#include <valarray>
#include <iostream>
#include <blitz/tinyvec.h>

/** Representation of a double valued 3x3 matrix. 
    @author John Shumway */
class Mat3 {
public:
#ifdef ENABLE_FLOAT
  typedef blitz::TinyVector<float,3> Vec3;
#else
  typedef blitz::TinyVector<double,3> Vec3;
#endif
  /** Constructor. */
  Mat3() : xx(0),xy(0),xz(0),yx(0),yy(0),yz(0),zx(0),zy(0),zz(0) {};
  /** The coordinates. */
  double xx,xy,xz,yx,yy,yz,zx,zy,zz;
   /** Allow printing of a matrix to an output stream. */
  friend std::ostream& operator<<(std::ostream&,const Mat3&);
  /** Add */
  Mat3& operator+=(const Mat3& m) {xx+=m.xx; xy+=m.xy; xz+=m.xz;
                                   yx+=m.yx; yy+=m.yy; yz+=m.yz;
                                   zx+=m.zx; zy+=m.zy; zz+=m.zz; return *this;}
  /** Subtract */
  Mat3& operator-=(const Mat3& m) {xx-=m.xx; xy-=m.xy; xz-=m.xz;
                                   yx-=m.yx; yy-=m.yy; yz-=m.yz;
                                   zx-=m.zx; zy-=m.zy; zz-=m.zz; return *this;}
  /** Scaling by a scalar. */
  Mat3& operator*=(const double s) {xx*=s; xy*=s; xz*=s;
                                    yx*=s; yy*=s; yz*=s;
                                    zx*=s; zy*=s; zz*=s; return *this;}
  /** Trace. */
  double trace() {return xx+yy+zz;}
  /** Tensor (outer) product. */
  Mat3& tensorProduct(const Vec3& u, const Vec3& v) {
    xx=u[0]*v[0]; xy=u[0]*v[1]; xz=u[0]*v[2];
    yx=u[1]*v[0]; yy=u[1]*v[1]; yz=u[1]*v[2];
    zx=u[2]*v[0]; zy=u[2]*v[1]; zz=u[2]*v[2]; return *this;}
};

#endif
