// $Id: main.cc,v 1.16 2007/12/28 02:19:40 jshumwa Exp $
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
#include <iostream>
#include <iomanip>
#include <valarray>
#include <cstdlib>
#include "CellRelaxer.h"
#include "CoordsH5.h"
#include "SpeciesH5.h"
#include "NeighborsH5.h"
#include "SuperCellH5.h"
#include "StillWeb.h"
#include "Spring.h"
#include "TotalEnergy.h"
#include "Forces.h"
#include "ConjGrad.h"
#include "TestForce.h"
#include "Checkpoint.h"
#include "InputParser.h"
#include "Stress.h"
#include "SurfaceStress.h"
#include "StressH5.h"
/** @mainpage RELAX code
  * @section overview Overview
  * <p>This program relaxes atomic nanostructures using empirical 
  * interatomic potentials.</p>
  * @section outline Outline of code
  * <p>The structure is read in using the StructH5 subclasses SpeciesH5, 
  * CoordsH5, NeighborsH5, and SuperCellH5. </p>
  * <p>An AtomicPotential object is created, namely StillWeb (a hard-coded
  * Si/Ge Stillinger-Webber potential).  The structure info and
  * AtomicPotential are fed into the TotalEnergy and Forces objects.</p>
  * <p>The structure is minimized using with conjugate-gradients
  * using a ConjGrad object.  Periodic checkpointing is accomplished with
  * a Checkpoint object.  The relaxed structure can be found in the
  * checkpoint output is the relaxed structure.</p>
  * <p>See InputParser for a sample XML Input file </p>
  * @todo Add XML input of potential parameters.
  * @todo Add generalized-VFF.
  * @author John Shumway */
int main(int argc, char** argv) {
  #pragma omp parallel
  std::cout.precision(10);
  std::cout << "relax: atomic strain relaxation code" << std::endl;
  int id=0, nthreads=1;
  #pragma omp parallel private(id)
  {
    #ifdef _OPENMP
    id = omp_get_thread_num();
    #endif
    printf("Hello World from thread %d\n", id);
    #pragma omp barrier
    if ( id == 0 ) {
      #ifdef _OPENMP
      nthreads = omp_get_num_threads();
      #endif
      printf("There are %d threads\n",nthreads);
    }
  }

  InputParser parser("relax.xml");
  parser.parse();
  const std::string infile(parser.getStructInFileName());
  std::cout << "Reading structure from " << infile << std::endl;
  SpeciesH5 species(infile);
  CoordsH5 coords(infile);
  SuperCellH5 cell(infile);
  NeighborsH5 nbr(infile);
  Checkpoint* checkpoint=parser.parseCheckpoint(coords, species, nbr, cell);
  if (!checkpoint) std::cout << "No checkpointing request" << std::endl;
  int natom=coords.coords.size();
  AtomicPotential *potential=parser.parsePotential();
  TotalEnergy en(coords,nbr,species,cell,*potential);
  Forces f(coords,nbr,species,cell,*potential);
  // Thermal anneal if requested.
  int annealSweeps=parser.getAnnealSweeps();
  if (annealSweeps>0) {
  
  }
  // Minimize using conjugate gradients.
  int maxiter=parser.getMaxIter();
  if (maxiter>0) {
    double ftol=parser.getFTol();
    CellRelaxer *cellRelaxer
           = parser.parseCellRelaxer(coords,nbr,species,cell,*potential,en);
    //cellRelaxer = new CellRelaxer(coords,nbr,species,cell,*potential,en,mode);
    ConjGrad conjGrad(f,en,coords,checkpoint,ftol,maxiter,cellRelaxer,
                      parser.getCellRelaxInterval(),parser.getScale());
    conjGrad.minimize();
    std::string outfile=parser.getRelaxedFilename();
    if (outfile!="") {
      std::string command="cp struct.h5 ";
      command+=outfile;
      system(command.c_str());
      cell.h5Write(outfile,cell.REWRITE);
      coords.h5Write(outfile,coords.REWRITE);
    }
  }
  // Calculate stresses.
  std::string stressfile=parser.getStressFilename();
  if (stressfile!="") {
    std::cout << "Calculating stresses" << std::endl;
    Stress stress(coords,nbr,species,cell,*potential);
    stress.compute();
    StressH5 stressh5(stress);
    std::string command="cp struct.h5 ";
    command+=stressfile;
    system(command.c_str());
    cell.h5Write(stressfile,cell.REWRITE);
    coords.h5Write(stressfile,coords.REWRITE);
    stressh5.h5Write(stressfile,stressh5.REWRITE);
  }
  // Calculate stresses.
  std::string surfacefile=parser.getSurfaceFilename();
  if (surfacefile!="") {
    std::cout << "Calculating surface stresses" << std::endl;
    SurfaceStress surface(coords,nbr,species,cell,*potential);
  }

  // Test forces, but first randomize positions.
  /*for (int i=0; i<natom; ++i) {
    coords[i].x+=0.2*(0.5*-(double(rand())/RAND_MAX));
    coords[i].y+=0.2*(0.5*-(double(rand())/RAND_MAX));
    coords[i].z+=0.2*(0.5*-(double(rand())/RAND_MAX));
  }
  TestForce test(f,en,coords);
  test.test(); */
  delete checkpoint;
  delete potential;
}

