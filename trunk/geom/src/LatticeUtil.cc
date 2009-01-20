// $Id: LatticeUtil.cc,v 1.2 2004/04/22 20:09:44 jshumwa Exp $
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
#include "LatticeUtil.h"
#include "CoordsH5.h"
#include "LatticeIndexingH5.h"
#include "NeighborsH5.h"
#include "SuperCellH5.h"
#include "SpeciesH5.h"

LatticeUtil::LatticeUtil() : cell(0), species(0), indexing(0), coords(0),
  neighbor(0) {
}

LatticeUtil::~LatticeUtil() {
  delete cell;
  delete species;
  delete indexing;
  delete coords;
  delete neighbor;
}
