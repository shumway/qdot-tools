// $Id: GridFactory.h,v 1.1.1.1 2004/05/03 16:49:21 jshumwa Exp $
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
#ifndef __GridFactory_h_
#define __GridFactory_h_

#include "Grid.h"
#include <blitz/tinyvec.h>
/** Class for constructing grids.
@author John shumway
@version $Revision: 1.1.1.1 $
*/
class GridFactory{
public:
  // Typedefs.
  typedef blitz::TinyVector<int,3> IVec;
  typedef blitz::TinyVector<float,3> Vec;
  // Constructor.
  GridFactory(const IVec&, const Vec&);
  // Destructor.
  ~GridFactory() {}
  // Get a pointer to a new Grid.
  template <class T>
  Grid<T>* getNewGrid() const {return new Grid<T>(extent,delta);}
private:
  // Grid size.
  IVec extent;
  // Grid spacing.
  Vec delta;
};
#endif
