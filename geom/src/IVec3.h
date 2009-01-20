// $Id: IVec3.h,v 1.2 2004/04/22 20:09:44 jshumwa Exp $
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
#ifndef __IVec3_h_
#define __IVec3_h_

/** Representation of a integer valued 3-vector. 
    @author John Shumway */
class IVec3 {
public:
  /** Constructor. */
  IVec3() : x(0),y(0),z(0) {};
  /** Constructor. */
  IVec3(const int x, const int y, const int z) : x(x),y(y),z(z) {};
  /** The coordinates. */
  int x,y,z;
};
#endif
