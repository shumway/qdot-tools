// $Id: GridFactory.cc,v 1.1.1.1 2004/05/03 16:49:21 jshumwa Exp $
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
#include "GridFactory.h"
#include "Grid.h"

GridFactory::GridFactory(const IVec& extent, const Vec& delta)
  : extent(extent), delta(delta) {
  std::cout << "Will produce grids with extent " << extent << std::endl;
  std::cout << " and grid spacing " << delta << std::endl;
}
