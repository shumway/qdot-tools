// $Id: Vec3.h,v 1.3 2004/10/08 13:20:22 jshumwa Exp $
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
#ifndef __Vec3_h_
#define __Vec3_h_

/** Representation of a 3-vector. 
    @author John Shumway */
class Vec3 {
public:
  /** Constructor. */
  Vec3() : x(0),y(0),z(0) {};
  /** Constructor. */
  Vec3(const float x, const float y, const float z) : x(x),y(y),z(z) {};
  /** The coordinates. */
  float x,y,z;
};
#endif
