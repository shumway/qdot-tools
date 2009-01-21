// $Id: Stress.cc,v 1.5 2007/11/20 22:29:22 jshumwa Exp $
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
#include "Stress.h"
#include <iostream>
#include <cstring>
#include "CoordsH5.h"
#include "SpeciesH5.h"
#include "NeighborsH5.h"
#include "SuperCellH5.h"
#include "AtomicPotential.h"

Stress::Stress(const CoordsH5& coords, const NeighborsH5& nbr,
    const SpeciesH5& species, const SuperCellH5& cell, 
    const AtomicPotential& pot)
  : coords(coords), nbr(nbr), species(species), cell(cell), pot(pot),
    stress(coords.coords.size()), volume(coords.coords.size()) {
};
 

void Stress::compute() {
  int natom=stress.size();
  //for (int iatom=0; iatom<natom; ++iatom) stress[iatom]=0;
  // Replace with memset.
  memset(stress.data(),0,natom*sizeof(SymMat));
  // Loop over all atoms.
  for (int iatom=0; iatom<natom; ++iatom) {
    Vec3 pos=coords[iatom];
    int ispec=species.species[iatom]-1;
    // Find neigbors within cutoff.
    std::valarray<Vec3> delta(nbr.nmax);
    std::valarray<double> d(nbr.nmax);
    std::valarray<int> spec(nbr.nmax);
    std::valarray<int> id(nbr.nmax);
    int nnn=0;
    volume(iatom)=1;
    for (int i=0; i<nbr.nmax; ++i) {
      int inbr=nbr.getNeighbor(iatom,i)-1;
      if (inbr!=-1) {
        id[nnn]=inbr;
        delta[nnn]=coords[inbr];delta[nnn]-=pos;
        cell.pbc(delta[nnn]);
        d[nnn]=sqrt(dot(delta[nnn],delta[nnn]));
        if (d[nnn]<=pot.getCutoff(0,0)) {
          spec[nnn]=species.species[inbr]-1;
          volume(iatom)*=d[nnn];
          nnn++;
        }
      }
    }
    volume(iatom)=pow(volume(iatom),3./nnn)*8./sqrt(27);
    // Loop over bonds and accumulate stress.
    for (int j=0; j<nnn; ++j) {
      Vec3 rjhat=delta[j]; rjhat/=d[j];
      //Add strech force.
      double drF2=pot.drF2(d[j],ispec,spec[j]);
      Vec3 f=rjhat; f*=0.25*drF2; // f=-.5*force
      SymMat s; s.tensorProduct(f,delta[j]);
      stress(id[j])+=s;
      stress(iatom)+=s;
      for (int k=0; k<j; ++k) {
        Vec3 rkhat=delta[k]; rkhat/=d[k];
        //Add bend force (includes additional stretch parts).
        double costheta=dot(delta[j],delta[k])/(d[j]*d[k]);
        Vec3 rjkhat=rjhat; rjkhat*=costheta; rjkhat-=rkhat;
        Vec3 rkjhat=rkhat; rkjhat*=costheta; rkjhat-=rjhat;
        double dr1H=0, dr2H=0, dthetaH=0;
        pot.dh(d[j],d[k],costheta,ispec,spec[j],spec[k],dr1H,dr2H,dthetaH);
        //Derivative of bend potential wrt r1.
        f=rjhat; f*=0.5*dr1H;
        s.tensorProduct(f,delta[j]);
        stress(id[j])+=s;
        stress(iatom)+=s;
        //Derivative of bend potential wrt r2.
        f=rkhat; f*=0.5*dr2H;
        s.tensorProduct(f,delta[k]);
        stress(id[k])+=s;
        stress(iatom)+=s;
        //Derivative of bend potential wrt costheta;
        f=rjkhat; f*=-0.5*dthetaH/d[j];
        s.tensorProduct(f,delta[j]);
        stress(id[j])+=s;
        stress(iatom)+=s;
        f=rkjhat; f*=-0.5*dthetaH/d[k];
        s.tensorProduct(f,delta[k]);
        stress(id[k])+=s;
        stress(iatom)+=s;
      }
    }
  } 
  stress/=volume;
}

