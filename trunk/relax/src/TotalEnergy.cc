// $Id: TotalEnergy.cc,v 1.7 2007/03/14 19:47:50 jshumwa Exp $
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
#ifdef _OPENMP
#include <omp.h>
#endif
#include "TotalEnergy.h"
#include <iostream>
#include "CoordsH5.h"
#include "SpeciesH5.h"
#include "NeighborsH5.h"
#include "SuperCellH5.h"
#include "StillWeb.h"
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>

TotalEnergy::TotalEnergy(const CoordsH5& coords, const NeighborsH5& nbr,
    const SpeciesH5& species, const SuperCellH5& cell, 
    const AtomicPotential& pot)
  : coords(coords), nbr(nbr), species(species), cell(cell), pot(pot) {
//#ifdef _OPENMP
//    numThreads(omp_get_num_threads()),
//#else
//    numThreads(1),
//#endif
//    delta(nbr.nmax,numThreads), d(nbr.nmax,numThreads), 
//    spec(nbr.nmax,numThreads), igej(nbr.nmax,numThreads) {
}
 

void TotalEnergy::compute() {
  int natom=coords.coords.size();
  double stretchEnergy=0;
  double bindingEnergy=0;
  double bendEnergy=0;
  #pragma omp parallel
  {
#ifdef ENABLE_FLOAT
  typedef blitz::TinyVector<float,3> Vec3;
#else
  typedef blitz::TinyVector<double,3> Vec3;
#endif
  blitz::Array<Vec3,1> delta(nbr.nmax);
  blitz::Array<double,1> d(nbr.nmax); 
  blitz::Array<int,1> spec(nbr.nmax);
  blitz::Array<bool,1> igej(nbr.nmax);
  // Loop over all atoms.
  #pragma omp for reduction(+:stretchEnergy,bindingEnergy,bendEnergy)
  for (int iatom=0; iatom<natom; ++iatom) {
//    int ithread=0;
//#ifdef _OPENMP
//    ithread=omp_get_thread_num();
//#endif
//std::cout << ithread << " of " << numThreads << std::endl;
    Vec3 pos=coords[iatom];
    int ispec=species.species[iatom]-1;
    // Find neigbors within cutoff.
    int nnn=0;
    for (int i=0; i<nbr.nmax; ++i) {
      int inbr=nbr.getNeighbor(iatom,i)-1;
      if (inbr!=-1) {
        delta(nnn)=coords[inbr];delta(nnn)-=pos;
        cell.pbc(delta(nnn));
        d(nnn)=sqrt(dot(delta(nnn),delta(nnn)));
        spec(nnn)=species.species[inbr]-1;
        igej(nnn)= (iatom>inbr);
        if (d(nnn)<=pot.getCutoff(ispec,spec(nnn))) nnn++;
      }
    }
    // Loop over bonds and accumulate energy.
    for (int j=0; j<nnn; ++j) {
      if (igej(j)) {
        bindingEnergy+=pot.bind(ispec,spec(j));
        stretchEnergy+=pot.f2(d(j),ispec,spec(j));
      }
      for (int k=0; k<j; ++k) {
        double costheta=dot(delta(j),delta(k))/(d(j)*d(k));
        bendEnergy+=pot.h(d(j),d(k),costheta,ispec,spec(j),spec(k));
      }
    }
  }
  }
  this->stretchEnergy=stretchEnergy;
  this->bindingEnergy=bindingEnergy;
  this->bendEnergy=bendEnergy;
  energy=stretchEnergy+bendEnergy;
}

void TotalEnergy::report() {
  std::cout << "Stretch energy = " << stretchEnergy*27.211 << " eV"<< std::endl;
  std::cout << "Bend energy = " << bendEnergy*27.211 << " eV" << std::endl;
  std::cout << "Binding energy = " << bindingEnergy*27.211 << " eV" << std::endl;
  std::cout << "Total energy = " << (bindingEnergy+energy)*27.211
            << " eV" << std::endl;
}

