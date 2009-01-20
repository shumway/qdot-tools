// $Id: EtchedLattice.cc,v 1.2 2004/06/11 01:09:10 jshumwa Exp $
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
#include "EtchedLattice.h"
#include "CoordsH5.h"
#include "LatticeIndexingH5.h"
#include "NeighborsH5.h"
#include "SuperCellH5.h"
#include "SpeciesH5.h"
#include "Vec3.h"
#include "IVec3.h"
#include <blitz/array.h>

EtchedLattice::EtchedLattice(const LatticeUtil& original,
                             const int etchedAtomID) {
  std::cout << "Etching atomID " << etchedAtomID << std::endl;
  // First count the number of remaining atoms.
  int natoms=original.getCoords().coords.size();
  int nspecies=original.getSpecies().name.size();
  int nnmax=original.getNeighbors().nmax;
std::cout << natoms << "," << nspecies << "," << nnmax << std::endl;
  int nremain=0;
  for (int i=0; i<natoms; ++i) {
    if (original.getSpecies().species[i]!=etchedAtomID) ++nremain;
  }
std::cout << nremain << std::endl;
  // Now allocate new structures.
  cell = new SuperCellH5(original.getSuperCell());
  species = new SpeciesH5(nspecies-1,nremain);
  for (int ispec=0; ispec<nspecies; ++ispec) {
    if (ispec+1<etchedAtomID) {
      species->name[ispec]=original.getSpecies().name[ispec];
    } else {
      species->name[ispec-1]=original.getSpecies().name[ispec];
    }
  }
  indexing = new LatticeIndexingH5(original.getLatticeIndexing());
  indexing->index.resize(nremain);
  coords = new CoordsH5(nremain);
  neighbor = new NeighborsH5(nnmax,nremain);
  // Setup index for remaining atoms.
  blitz::Array<int,1> index(nremain), revIndex(natoms); revIndex=-1;
  int iatom=0;
  for (int i=0; i<natoms; ++i) {
    if (original.getSpecies().species[i]!=etchedAtomID) {
      index(iatom)=i; revIndex(i)=iatom++;
    }
  }
  // Now copy data into new structures.
  for (int iatom=0; iatom<nremain; ++iatom) {
    coords->coords[iatom]=original.getCoords().coords[index(iatom)];
    int ispec=original.getSpecies().species[index(iatom)];
    species->species[iatom]=(ispec<etchedAtomID)?ispec:ispec-1;
    indexing->index[iatom]=original.getLatticeIndexing().index[index(iatom)];
    for (int ii=0; ii<nnmax; ++ii) {
      int id=original.getNeighbors().nn[index(iatom)*nnmax+ii];
      neighbor->nn[iatom*nnmax+ii]=revIndex(id-1)+1;
    }
  }
}
