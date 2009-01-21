// $Id: Checkpoint.cc,v 1.7 2004/08/14 19:45:59 jshumwa Exp $
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
#include "Checkpoint.h"
#include "CoordsH5.h"
#include "SpeciesH5.h"
#include "SuperCellH5.h"
#include "NeighborsH5.h"
#include <iostream>

Checkpoint::Checkpoint(const CoordsH5& coords, const SpeciesH5& species, 
  const NeighborsH5& neighbors, const SuperCellH5& cell, 
  const std::string& fn,const int interval)
  : coords(coords), species(species), neighbors(neighbors), cell(cell),
    filename(fn), interval(interval) {
  std::cout << "Checkpointing to file " << filename
            << " at an interval " << interval << std::endl;
}

void Checkpoint::write(int i) const {
  if(i==0){
    std::string command="cp struct.h5 ";
    command+=filename;
    system(command.c_str());
    return;
  }
  if ((i+1)%interval != 0) return;
  cell.h5Write(filename,cell.REWRITE);
  coords.h5Write(filename,coords.REWRITE);
}
