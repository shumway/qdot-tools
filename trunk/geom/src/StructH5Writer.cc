// $Id: StructH5Writer.cc,v 1.5 2005/01/07 17:47:28 jshumwa Exp $
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
#include "StructH5Writer.h"
#include "SuperCell.h"
#include "CoordsH5.h"
#include "LatticeIndexingH5.h"
#include "Nanostructure.h"
#include "NeighborsH5.h"
#include "Structure.h"
#include "SpeciesH5.h"
#include "SuperCellH5.h"
#include "ZincBlende.h"
#include "EtchedLattice.h"
#include "Reconstructor.h"
#include "Binary.h"
#include "CommentH5.h"
#include "EtchedLattice.h"
#include <fstream>

StructH5Writer::StructH5Writer(const SuperCell& cell, const Nanostructure& ns) :
  cell(cell), structs(ns), lattice(0) {

  // Create the bulk zinc blende lattice.
  int nspecies=ns.getNSpecies();

/*  zb = new ZincBlende(cell.getNX(),cell.getNY(),cell.getNZ(),cell.getA(),
                 ns.getSpeciesName(cell.getBulkMaterial()->getAnionIndex()-1),
                 ns.getSpeciesName(cell.getBulkMaterial()->getCationIndex()-1),
                 nspecies, 
                 cell.getBulkMaterial()->getAnionIndex(),
                 cell.getBulkMaterial()->getCationIndex()); */
  lattice = new ZincBlende(cell,ns);

  // Get the coordinates and species for the bulk zinc blend lattice.
  LatticeIndexingH5& indexingH5(lattice->getLatticeIndexing());
  SpeciesH5& speciesH5(lattice->getSpecies());

  // Place the nanostructures in the bulk lattice.
  for (int i=0; i<ns.getNumber(); ++i) {
    ns.getStructure(i)->build(speciesH5,indexingH5);
  }
} 

StructH5Writer::~StructH5Writer() {
 delete lattice;
}

void StructH5Writer::write(const std::string& filename) {
  lattice->getSuperCell().h5Write(filename,StructH5::NEW);
  lattice->getCoords().h5Write(filename);
  lattice->getSpecies().h5Write(filename);
  lattice->getLatticeIndexing().h5Write(filename);
  lattice->getNeighbors().h5Write(filename);
}

void StructH5Writer::etch(const int etchIndex) {
  LatticeUtil *etched = new EtchedLattice(*lattice,etchIndex);
  delete lattice;
  lattice = etched;
}

void StructH5Writer::reconstruct() {
  LatticeUtil *reconstructed = new Reconstructor(*lattice);
  delete lattice;
  lattice = reconstructed;
}
