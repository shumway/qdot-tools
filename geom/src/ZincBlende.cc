// $Id: ZincBlende.cc,v 1.3 2006/08/08 18:02:48 jshumwa Exp $
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
#include "ZincBlende.h"
#include "CoordsH5.h"
#include "LatticeIndexingH5.h"
#include "Nanostructure.h"
#include "NeighborsH5.h"
#include "Material.h"
#include "SuperCell.h"
#include "SuperCellH5.h"
#include "SpeciesH5.h"
#include "Vec3.h"
#include "IVec3.h"

ZincBlende::ZincBlende(/*const int  nx, const int ny, const int nz, 
                       const double a, const std::string& anionName,
                       const std::string& cationName, const int nspecies,
                       const int anionID, const int cationID) { */
const SuperCell& cellInfo, const Nanostructure& ns) {
  int nx=cellInfo.getNX();
  int ny=cellInfo.getNY();
  int nz=cellInfo.getNZ();
  double a=cellInfo.getA();
  double strainxx=cellInfo.getStrain(0);
  double strainyy=cellInfo.getStrain(1);
  double strainzz=cellInfo.getStrain(2);
  int natoms = nx*ny*nz*8;

  // Set up the supercell.
  cell = new SuperCellH5(); 
  cell->a1[0]=nx*a*(1+strainxx);
  cell->a2[1]=ny*a*(1+strainyy);
  cell->a3[2]=nz*a*(1+strainzz);

  // Allocate the species table.
  int nspecies=ns.getNSpecies();
  species = new SpeciesH5(nspecies,natoms);
  // Set the species names.
  for (int ispecies=0; ispecies<nspecies; ++ispecies) {
    species->name[ispecies]=ns.getSpeciesName(ispecies);
  }

  // Allocate the lattice indexing table.
//  indexing = new LatticeIndexingH5(natoms,8,offset);
//  indexing->nCell[0]=nx; indexing->nCell[1]=ny; indexing->nCell[2]=nz;
  Vec3 a1(0.5,0.5,0.),a2(0,0.5,0.5),a3(0.5,0.,0.5);
  Vec3 r0(0.,0.,0.),r1(0.25,0.25,0.25);
  std::valarray<Vec3> offset(2); offset[0]=r0;offset[1]=r1;
  IVec3 n1(nx,-nx,nx),n2(ny,ny,-ny),n3(-nz,nz,nz);
  indexing = new LatticeIndexingH5(a,a1,a2,a3,offset,natoms,n1,n2,n3);

  // Allocate the atomic coordinates.
  coords = new CoordsH5(natoms);

  // Create the lattice.
  const Material& material(*cellInfo.getBulkMaterial());
  int iatom=0;
  Vec3 pt;
  LatticeIndexingH5::indices ind;
  LatticeIndexingH5::indices prim_ind;
  for (ind.k=0; ind.k<nz; ++ind.k) {
    for (ind.j=0; ind.j<ny; ++ind.j) {
      for (ind.i=0; ind.i<nx; ++ind.i) {
        for (ind.ipos=0; ind.ipos<8; ++ind.ipos) {  
          //species->species[iatom] = (ind.ipos & 2)==0 ? anionID : cationID;
          prim_ind.i=ind.i+ind.j-ind.k;
          prim_ind.j=ind.j+ind.k-ind.i;
          prim_ind.k=ind.k+ind.i-ind.j;
          switch (ind.ipos) {
          case (0): {
            prim_ind.ipos=0; break;} 
          case (1): {
            prim_ind.ipos=0; prim_ind.i++; break;} 
          case (2): {
            prim_ind.ipos=1; break;} 
          case (3): {
            prim_ind.ipos=1; prim_ind.i++; break;} 
          case (4): {
            prim_ind.ipos=0; prim_ind.k++; break;} 
          case (5): {
            prim_ind.ipos=0; prim_ind.j++; break;} 
          case (6): {
            prim_ind.ipos=1; prim_ind.k++; break;} 
          case (7): {
            prim_ind.ipos=1; prim_ind.j++; break;} 
          }
          indexing->index[iatom]=prim_ind;
          pt.x = a*(1+strainxx)
                  *(0.5*(prim_ind.i+prim_ind.k)+offset[prim_ind.ipos].x);
          pt.y = a*(1+strainyy)
                  *(0.5*(prim_ind.i+prim_ind.j)+offset[prim_ind.ipos].y);
          pt.z = a*(1+strainzz)
                  *(0.5*(prim_ind.j+prim_ind.k)+offset[prim_ind.ipos].z);
          coords->coords[iatom] = pt;
          species->species[iatom] = (prim_ind.ipos==0)
                                    ? material.getAnionIndex(pt)
                                    : material.getCationIndex(pt);
//cout << prim_ind.i << "," <<  prim_ind.j << "," <<  prim_ind.k << ": " 
//     << pt.x << "," <<  pt.y << "," <<  pt.z  << endl;
          ++iatom;
        }
      }
    }
  }

  // Generate the nearest-neighbor lists directly from the cell indexing.
  neighbor = findNeighbors(*indexing);
}

NeighborsH5* 
ZincBlende::findNeighbors(const LatticeIndexingH5& indexing) const {
  int natoms = indexing.index.size();
  NeighborsH5* nbr = new NeighborsH5(4,natoms);
  int ncell[3];
  ncell[0]=indexing.n1.x;
  ncell[1]=indexing.n2.y;
  ncell[2]=indexing.n3.z;
  // First generate an ordered neighbor lists.
  // (There are 23 atoms that are nearest neighbors to the eight atoms in
  // each base cell.  Construct the ordered 23 atom list for each cell.)
  std::valarray<int> nlist(23*ncell[0]*ncell[1]*ncell[2]);
  for (int iatom=1; iatom<=natoms; ++iatom) {
    // Generate a base index from primative indexing.
    const LatticeIndexingH5::indices& prim_ind(indexing.index[iatom-1]);
    LatticeIndexingH5::indices ind;
    ind.i = (prim_ind.i+prim_ind.k)/2;
    ind.j = (prim_ind.i+prim_ind.j)/2;
    ind.k = (prim_ind.j+prim_ind.k)/2;
    switch ((prim_ind.i&1) + 2*(prim_ind.j&1) + 4*(prim_ind.k&1)) {
    case 0: case 7:
      ind.ipos = (prim_ind.ipos==0) ? 0 : 2;
      break;
    case 1: case 6:
      ind.ipos = (prim_ind.ipos==0) ? 1 : 3;
      break;
    case 2: case 5:
      ind.ipos = (prim_ind.ipos==0) ? 5 : 7;
      break;
    case 3: case 4:
      ind.ipos = (prim_ind.ipos==0) ? 4 : 6;
      break;
    };
    switch (ind.ipos) {
    case 0: {
      int ii = (ind.i+ncell[0]-1)%ncell[0];
      int jj = (ind.j+ncell[1]-1)%ncell[1];
      int kk = (ind.k+ncell[2]-1)%ncell[2];
      nlist[0+23*(ind.i+ncell[0]*(ind.j+ncell[1]*ind.k))]=iatom;
      nlist[15+23*(ii+ncell[0]*(jj+ncell[1]*ind.k))]=iatom;
      nlist[20+23*(ii+ncell[0]*(ind.j+ncell[1]*kk))]=iatom;
      nlist[22+23*(ind.i+ncell[0]*(jj+ncell[1]*kk))]=iatom;
      break;}
    case 1: {
      int kk = (ind.k+ncell[2]-1)%ncell[2];
      nlist[1+23*(ind.i+ncell[0]*(ind.j+ncell[1]*ind.k))]=iatom;
      nlist[21+23*(ind.i+ncell[0]*(ind.j+ncell[1]*kk))]=iatom;
      break;}
    case 2: {
      nlist[2+23*(ind.i+ncell[0]*(ind.j+ncell[1]*ind.k))]=iatom;
      break;}
    case 3: {
      int ii = (ind.i+1)%ncell[0];
      int jj = (ind.j+1)%ncell[1];
      nlist[3+23*(ind.i+ncell[0]*(ind.j+ncell[1]*ind.k))]=iatom;
      nlist[8+23*(ii+ncell[0]*(jj+ncell[1]*ind.k))]=iatom;
      nlist[17+23*(ind.i+ncell[0]*(jj+ncell[1]*ind.k))]=iatom;
      nlist[19+23*(ii+ncell[0]*(ind.j+ncell[1]*ind.k))]=iatom;
      break;}
    case 4: {
      int jj = (ind.j+ncell[1]-1)%ncell[1];
      nlist[4+23*(ind.i+ncell[0]*(ind.j+ncell[1]*ind.k))]=iatom;
      nlist[14+23*(ind.i+ncell[0]*(jj+ncell[1]*ind.k))]=iatom;
      break;}
    case 5: {
      int ii = (ind.i+ncell[0]-1)%ncell[0];
      nlist[5+23*(ind.i+ncell[0]*(ind.j+ncell[1]*ind.k))]=iatom;
      nlist[13+23*(ii+ncell[0]*(ind.j+ncell[1]*ind.k))]=iatom;
      break;}
    case 6: {
      int ii = (ind.i+1)%ncell[0];
      int kk = (ind.k+1)%ncell[2];
      nlist[6+23*(ind.i+ncell[0]*(ind.j+ncell[1]*ind.k))]=iatom;
      nlist[10+23*(ii+ncell[0]*(ind.j+ncell[1]*kk))]=iatom;
      nlist[11+23*(ind.i+ncell[0]*(ind.j+ncell[1]*kk))]=iatom;
      nlist[18+23*(ii+ncell[0]*(ind.j+ncell[1]*ind.k))]=iatom;
      break;}
    case 7: {
      int jj = (ind.j+1)%ncell[1];
      int kk = (ind.k+1)%ncell[2];
      nlist[7+23*(ind.i+ncell[0]*(ind.j+ncell[1]*ind.k))]=iatom;
      nlist[9+23*(ind.i+ncell[0]*(jj+ncell[1]*kk))]=iatom;
      nlist[12+23*(ind.i+ncell[0]*(ind.j+ncell[1]*kk))]=iatom;
      nlist[16+23*(ind.i+ncell[0]*(jj+ncell[1]*ind.k))]=iatom;
      break;}
    }
  }

  // Use ordered neighbor lists to fill in nearest neighbor lists.
  // (For each of the eight atoms in the cell, use the nn table to pick out
  // the 4 neighbors for the atom from list of 23 neighbors for the cell.)
  for (int iatom=0; iatom<natoms; ++iatom) {
    const int nn[] = { 2, 8, 9,10,  3, 2,11,12,  4, 5, 1, 0,  13,14,15, 1,
                       6,16,17, 2,  7,18, 2,19, 20,21,13, 4,  21,22,14, 5};
    // Generate a base index from primative indexing.
    const LatticeIndexingH5::indices& prim_ind(indexing.index[iatom]);
    LatticeIndexingH5::indices ind;
    ind.i = (prim_ind.i+prim_ind.k)/2;
    ind.j = (prim_ind.i+prim_ind.j)/2;
    ind.k = (prim_ind.j+prim_ind.k)/2;
    switch ((prim_ind.i&1) + 2*(prim_ind.j&1) + 4*(prim_ind.k&1)) {
    case 0: case 7: {
      ind.ipos = (prim_ind.ipos==0) ? 0 : 2;
      break;}
    case 1: case 6: {
      ind.ipos = (prim_ind.ipos==0) ? 1 : 3;
      break;}
    case 2: case 5: {
      ind.ipos = (prim_ind.ipos==0) ? 5 : 7;
      break;}
    case 3: case 4: {
      ind.ipos = (prim_ind.ipos==0) ? 4 : 6;
      break;}
    }
    for (int irel=0; irel<4; ++irel) {
      nbr->nn[irel+4*iatom] = nlist[nn[irel+4*ind.ipos] 
                               +23*(ind.i+ncell[0]*(ind.j+ncell[1]*ind.k))];
    }
  }
  return nbr;
}
