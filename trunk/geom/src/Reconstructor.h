// $Id: Reconstructor.h,v 1.1 2005/01/07 17:47:28 jshumwa Exp $
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
#ifndef __Reconstructor_h_
#define __Reconstructor_h_

#include <string>
#include "LatticeUtil.h"
#include <blitz/array.h>

/** Utility to etch dummy atoms from structure for free-standing or
    surface nanostructures.
    @author John Shumway */
class Reconstructor : public LatticeUtil {
public:
  /** Constructor. */
  Reconstructor(const LatticeUtil&);
private:
  void findNewNeighbors(const blitz::Array<int,3>&, 
                        const int ix, const int iy, const int ithis);
};
#endif
