// $Id: Mat3.cc,v 1.2 2004/08/14 01:42:07 jshumwa Exp $
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
#include "Mat3.h"


std::ostream& operator<<(std::ostream& os, const Mat3& m) {
  return os << "(" << m.xx << "," << m.xy << "," << m.xz << ")\n"
            << "(" << m.yx << "," << m.yy << "," << m.yz << ")\n"
            << "(" << m.zx << "," << m.zy << "," << m.zz << ")" << std::endl;
}

