// $Id: CellRelaxer.cc,v 1.3 2007/03/15 19:03:00 jshumwa Exp $
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
#include "CellRelaxer.h"
#include <iostream>
#include "CoordsH5.h"
#include "SpeciesH5.h"
#include "NeighborsH5.h"
#include "SuperCellH5.h"
#include "TotalEnergy.h"
#include "AtomicPotential.h"

CellRelaxer::CellRelaxer(CoordsH5& coords, const NeighborsH5& nbr,
    const SpeciesH5& species, SuperCellH5& cell, 
    const AtomicPotential& pot, TotalEnergy& en, const int mode)
  : minimizer(*this),
    coords(coords), nbr(nbr), species(species), 
    cell(cell), pot(pot), en(en), mode(mode) {
};
 

void CellRelaxer::relax() {
  en.compute();
  SymMat stress=computeStress();
  lineMinimize();
  stress=computeStress();
}

double CellRelaxer::operator()(const double x) {
  int natom=coords.coords.size();
  if (x!=xLast) {
    double scale=(1.0+x)/(1.0+xLast);
    if (mode==ZONLY) {
      for (int iatom=0; iatom<natom;++iatom) coords[iatom][2]*=scale;
      cell.a1[2]*=scale; cell.a2[2]*=scale; cell.a3[2]*=scale;
    } else if (mode==VOLUME) {
      for (int iatom=0; iatom<natom;++iatom) coords[iatom]*=scale;
      cell.a1*=scale; cell.a2*=scale; cell.a3*=scale;
    }
    cell.computeRecipricalVectors();
    xLast=x;
    en.compute();
  }
  return en.getEnergy();
}

void CellRelaxer::lineMinimize() {
  // Start with points ax and bx.
  double ax=xLast=0.0;
  double ena=(*this)(ax);
  double bx=0.01;
  double enb=(*this)(bx);
  std::cout << "ax="<<ax<<", bx="<<bx<< std::endl;
  std::cout << "ena="<<27.211*ena<<", enb="<<27.211*enb<< std::endl;
  // Bracket the minimum.
  double cx=0,enc=0;
  std::cout << "Bracketing minimum." << std::endl;
  minimizer.bracketMin(ax,bx,cx,ena,enb,enc);
  std::cout << "ax="<<ax<<", bx="<<bx<<", cx="<<cx << std::endl;
  std::cout << "ena="<<27.211*ena<<", enb="<<27.211*enb<<", enc="<<27.211*enc <<
 std::endl;
  // Find solution using golden section search.
  double xmin, emin;
  std::cout << "Searching for minimum." << std::endl;
  //emin = minimizer.goldenSearch(ax,bx,cx,xmin);
  emin = minimizer.brentSearch(ax,bx,cx,xmin);
  (*this)(xmin);
  en.compute();
  en.report();
}

CellRelaxer::SymMat CellRelaxer::computeStress() {
  std::cout << "Calculating stress on SuperCell." << std::endl;
  int natom=coords.coords.size();
  SymMat stress; stress*=0.0;
  #pragma omp parallel
  {
  SymMat mystress; mystress*=0.0;
  std::valarray<Vec3> delta(nbr.nmax);
  std::valarray<double> d(nbr.nmax);
  std::valarray<int> spec(nbr.nmax);
  std::valarray<int> id(nbr.nmax);
  // Loop over all atoms.
  #pragma omp for
  for (int iatom=0; iatom<natom; ++iatom) {
    Vec3 pos=coords[iatom];
    int ispec=species.species[iatom]-1;
    // Find neigbors within cutoff.
    int nnn=0;
    //volume(iatom)=1;
    for (int i=0; i<nbr.nmax; ++i) {
      int inbr=nbr.getNeighbor(iatom,i)-1;
      if (inbr!=-1) {
        id[nnn]=inbr;
        delta[nnn]=coords[inbr];delta[nnn]-=pos;
        cell.pbc(delta[nnn]);
        d[nnn]=sqrt(dot(delta[nnn],delta[nnn]));
        if (d[nnn]<=pot.getCutoff(0,0)) {
          spec[nnn]=species.species[inbr]-1;
          //volume(iatom)*=d[nnn];
          nnn++;
        }
      }
    }
    //volume(iatom)=pow(volume(iatom),3./nnn)*8./sqrt(27);
    // Loop over bonds and accumulate stress.
    for (int j=0; j<nnn; ++j) {
      Vec3 rjhat=delta[j]; rjhat/=d[j];
      //Add strech force.
      double drF2=pot.drF2(d[j],ispec,spec[j]);
      Vec3 f=rjhat; f*=0.25*drF2; // f=-.5*force
      SymMat s; s.tensorProduct(f,delta[j]);
      mystress+=2*s;
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
        mystress+=2*s;
        //Derivative of bend potential wrt r2.
        f=rkhat; f*=0.5*dr2H;
        s.tensorProduct(f,delta[k]);
        mystress+=2*s;
        //Derivative of bend potential wrt costheta;
        f=rjkhat; f*=-0.5*dthetaH/d[j];
        s.tensorProduct(f,delta[j]);
        mystress+=2*s;
        f=rkjhat; f*=-0.5*dthetaH/d[k];
        s.tensorProduct(f,delta[k]);
        mystress+=2*s;
      }
    }
  } 
  #pragma omp critical
  stress+=mystress;
  }
  //for (int iatom=0; iatom<natom; ++iatom) stress[iatom]*=(1./volume[iatom]);
  double volume=dot(cell.a1,cross(cell.a2,cell.a3));
  stress/=volume;
  std::cout << stress << std::endl;
  return stress;
}

