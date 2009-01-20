// $Id: Reconstructor.cc,v 1.2 2005/01/21 01:01:16 jshumwa Exp $
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
#include "Reconstructor.h"
#include "CoordsH5.h"
#include "LatticeIndexingH5.h"
#include "NeighborsH5.h"
#include "SuperCellH5.h"
#include "SpeciesH5.h"
#include "Vec3.h"
#include "IVec3.h"
#include <blitz/array.h>

Reconstructor::Reconstructor(const LatticeUtil& original) {
std::cout << "reconstructing surface" << std::endl;
  // Now allocate new structures.
  cell = new SuperCellH5(original.getSuperCell());
  species = new SpeciesH5(original.getSpecies());
  indexing = new LatticeIndexingH5(original.getLatticeIndexing());
  coords = new CoordsH5(original.getCoords());
  neighbor = new NeighborsH5(original.getNeighbors());
  const int natom = coords->coords.size();
  const int nmax = neighbor->nmax;
  std::valarray<int>& nn(neighbor->nn);
  std::valarray<LatticeIndexingH5::indices>& index(indexing->index);
  blitz::Array<int,1> dist(nmax+1);
  dist=0;
  // Find new, reconstructed neighbors.
  const int nx=indexing->n1.x, ny=indexing->n2.y;
  blitz::Array<int,3> table(nx,ny,16);
  table=-1;
  for (int iatom=0; iatom<natom; ++iatom) {
    int nnbr=0;
    for (int i=0; i<nmax; ++i) nnbr += (nn[i+4*iatom]>0)?1:0;
    dist(nnbr)++;
    if (nnbr<4) {
      //std::cout << index[iatom].i << ", "
      //          << index[iatom].j << ", "
      //          << index[iatom].k << ", "
      //          << index[iatom].ipos << "    ";
      int ix = 2*(index[iatom].i+index[iatom].k);
      int iy = 2*(index[iatom].i+index[iatom].j);
      int iz = 2*(index[iatom].j+index[iatom].k);
      ix+=index[iatom].ipos; iy+=index[iatom].ipos; iz+=index[iatom].ipos;
      //std::cout << ix << ", " << iy << ", " << iz << std::endl;
      // put atom indices on a table.
      int i=0;
      while (table(ix/4,iy/4,i)!=-1) ++i;
      table(ix/4,iy/4,i) = iatom;
    } 
  }
  std::cout << dist << std::endl;
  for (int ix=0; ix<nx; ++ix) {
    for (int iy=0; iy<ny; ++iy) {
      for (int i=0; i<16; ++i) {
        if (table(ix,iy,i)>-1) {
          findNewNeighbors(table,ix,iy,i);
        }
      }
    }
  }
  dist=0;
  for (int iatom=0; iatom<natom; ++iatom) {
    int nnbr=0;
    for (int i=0; i<nmax; ++i) nnbr += (nn[i+4*iatom]>0)?1:0;
    dist(nnbr)++;
  }
  std::cout << dist << std::endl;
}

void Reconstructor::findNewNeighbors(const blitz::Array<int,3>& table,
    const int ix, const int iy, const int ithis) {
  const int iatom=table(ix,iy,ithis);
  std::valarray<int>& nn(neighbor->nn);
  const int nmax = neighbor->nmax;
  blitz::TinyVector<int,3> idim=table.shape(); 
  const int nx=idim[0],ny=idim[1];
  std::cout << "Finding neighbors for " << iatom << std::endl;
  LatticeIndexingH5::indices thisindex=indexing->index[iatom];
  // Find a reference bond.
  int katom=-1; for (int i=0; katom==-1; ++i) katom=nn[nmax*iatom+i]-1;
  LatticeIndexingH5::indices indexk=indexing->index[katom];
  int kx=2*(thisindex.i-indexk.i+thisindex.k-indexk.k)
          +(thisindex.ipos-indexk.ipos);
  int ky=2*(thisindex.i-indexk.i+thisindex.j-indexk.j)
          +(thisindex.ipos-indexk.ipos);
  int kz=2*(thisindex.j-indexk.j+thisindex.k-indexk.k)
          +(thisindex.ipos-indexk.ipos);
  kx = (kx+nx*6)%(nx*4)-nx*2;
  ky = (ky+ny*6)%(ny*4)-ny*2;
  // Search for neighbors. 
  for (int idx=-1; idx<=1; ++idx) {
    int ix2=(ix+nx+idx)%nx;
    for (int idy=-1; idy<=1; ++idy) {
      int iy2=(iy+ny+idy)%ny;
      for (int i=0; i<idim[2]; ++i) {
        const int jatom=table(ix2,iy2,i);
        if (jatom==-1) break;
        LatticeIndexingH5::indices index=indexing->index[jatom];
        int dx=2*(thisindex.i-index.i+thisindex.k-index.k)
                 +(thisindex.ipos-index.ipos);
        int dy=2*(thisindex.i-index.i+thisindex.j-index.j)
                 +(thisindex.ipos-index.ipos);
        int dz=2*(thisindex.j-index.j+thisindex.k-index.k)
                +(thisindex.ipos-index.ipos);
        dx = (dx+nx*6)%(nx*4)-nx*2;
        dy = (dy+ny*6)%(ny*4)-ny*2;

        if (dz==0 && abs(dx)==2 && abs(dy)==2 && dx*kx+dy*ky==0) {
          std::cout << iatom << ", " << jatom << ":\t";
          std::cout << dx << ", " << dy << ", " << dz  << "   "
                    << kx << ", " << ky << ", " << kz << std::endl;
          for (int ii=0; ii<nmax; ++ii) {
            if (nn[iatom*nmax+ii]==0) {
              nn[iatom*nmax+ii]=jatom+1;
              break;
            }
          }
        }
      }
    }
  }
}

