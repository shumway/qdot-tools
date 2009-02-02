// $Id: Forces.cc,v 1.6 2007/11/20 22:29:21 jshumwa Exp $
/*
    Copyright (C) 2004-7 John B. Shumway, Jr.

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
#include "Forces.h"
#include <iostream>
#include <valarray>
#include <cstring>
#include "CoordsH5.h"
#include "SpeciesH5.h"
#include "NeighborsH5.h"
#include "SuperCellH5.h"
#include "AtomicPotential.h"

Forces::Forces(const CoordsH5& coords, const NeighborsH5& nbr,
    const SpeciesH5& species, const SuperCellH5& cell, 
    const AtomicPotential& pot)
  : coords(coords), nbr(nbr), species(species), cell(cell), pot(pot),
    force(coords.coords.size()) {
};
 

void Forces::compute() {
  int natom=force.size();
  //for (int iatom=0; iatom<natom; ++iatom) force[iatom]=0;
  // Replace with memset.
  memset(&force[0],0,natom*sizeof(Vec3));
  const int NGROUP=100;
  #pragma omp parallel
  {
  std::valarray<Vec3> delta(nbr.nmax);
  std::valarray<double> d(nbr.nmax);
  std::valarray<int> spec(nbr.nmax);
  std::valarray<int> id(nbr.nmax);
#ifdef _OPENMP
  std::valarray<int> otherID(nbr.nmax*(1+2*(nbr.nmax-1))*NGROUP);
  std::valarray<Vec3> otherForce(nbr.nmax*(1+2*(nbr.nmax-1))*NGROUP);
#endif
  Vec3 pos, rjhat, f;
  // Loop over all atoms.
#pragma omp for
  for (int iatom0=0; iatom0<natom; iatom0+=NGROUP) {
#ifdef _OPENMP
  int nother=0;
#endif
  for (int iatom=iatom0; iatom<iatom0+NGROUP && iatom<natom; ++iatom) {
    pos=coords[iatom];
    int ispec=species.species[iatom]-1;
    // Find neigbors within cutoff.
    int nnn=0;
    for (int i=0; i<nbr.nmax; ++i) {
      int inbr=nbr.getNeighbor(iatom,i)-1;
      if (inbr!=-1) {
        id[nnn]=inbr;
        delta[nnn]=coords[inbr];delta[nnn]-=pos;
        cell.pbc(delta[nnn]);
        d[nnn]=sqrt(dot(delta[nnn],delta[nnn]));
        spec[nnn]=species.species[inbr]-1;
        if (d[nnn]<=pot.getCutoff(0,0)) nnn++;
      }
    }
    // Loop over bonds and accumulate force.
    for (int j=0; j<nnn; ++j) {
      rjhat=delta[j]; rjhat/=d[j];
      //Add strech force.
      double drF2=pot.drF2(d[j],ispec,spec[j]);
      f=rjhat; f*=-0.5*drF2;
#ifdef _OPENMP
      otherForce[nother]=f; otherID[nother++]=id[j];
#else
      force[id[j]]+=f;
#endif
      force[iatom]-=f;
      for (int k=0; k<j; ++k) {
        Vec3 rkhat=delta[k]; rkhat/=d[k];
        //Add bend force (includes additional stretch parts).
        double costheta=dot(delta[j],delta[k])/(d[j]*d[k]);
        Vec3 rjkhat=rjhat; rjkhat*=costheta; rjkhat-=rkhat;
        Vec3 rkjhat=rkhat; rkjhat*=costheta; rkjhat-=rjhat;
        double dr1H=0, dr2H=0, dthetaH=0;
        pot.dh(d[j],d[k],costheta,ispec,spec[j],spec[k],dr1H,dr2H,dthetaH);
        //Derivative of bend potential wrt r1.
        f=rjhat; f*=-dr1H;
#ifdef _OPENMP
        otherForce[nother]=f; otherID[nother++]=id[j];
#else
        force[id[j]]+=f;
#endif
        force[iatom]-=f;
        //Derivative of bend potential wrt r2.
        f=rkhat; f*=-dr2H;
#ifdef _OPENMP
        otherForce[nother]=f; otherID[nother++]=id[k];
#else
        force[id[k]]+=f;
#endif
        force[iatom]-=f;
        //Derivative of bend potential wrt costheta;
        f=rjkhat; f*=dthetaH/d[j];
#ifdef _OPENMP
        otherForce[nother]=f; otherID[nother++]=id[j];
#else
        force[id[j]]+=f;
#endif
        force[iatom]-=f;
        f=rkjhat; f*=dthetaH/d[k];
#ifdef _OPENMP
        otherForce[nother]=f; otherID[nother++]=id[k];
#else
        force[id[k]]+=f;
#endif
        force[iatom]-=f;
      }
    }
    }
#ifdef _OPENMP
#pragma omp critical
    for (int j=0; j<nother; ++j) force[otherID[j]]+=otherForce[j];
#endif
  } 
  } 
}

