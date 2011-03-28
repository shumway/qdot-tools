// $Id: main.cc,v 1.12 2007/12/06 23:59:59 jshumwa Exp $
/*
    Copyright (C) 2004,2008,2009 John B. Shumway, Jr.

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
#include <iostream>
#include <cstdlib>
#include <fstream>
#include "SpeciesH5.h"
#include "CoordsH5.h"
#include "NeighborsH5.h"
#include "LatticeIndexingH5.h"
#include "SuperCellH5.h"
#include "CommentH5.h"
#include <hdf5.h>
#include "Grid.h"
#include <blitz/array.h>
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>
#include "GridFactory.h"
#include "CompositionGrids.h"
#include "StressGrid.h"
#include "Stress.h"
#include "StressH5.h"
#include "Strain.h"
#include "StrainFromStress.h"
#include "ZincBlendeEMA.h"
#include "ZincBlendeElasticity.h"

/** @mainpage getema docs
  * @section overview Overview
  * <p>Compute band offsets and effective mass tensors from grided
  * composition and stress data. </p>
  * @section outline Outline of code
  * <p>The structure is read in using the StructH5 subclasses SpeciesH5, 
  * CoordsH5, and SuperCellH5. 
  * The output is an HDF5 file.</p>
  * @section Method
  * <p>The stress and strain tensors are related by the elastic moduli,
  * @f[  \sigma_{ij} = \sum_{kl} C_{ij\;kl} \epsilon_{kl},  @f] 
  * and elastic constants,
  * @f[  \epsilon_{ij} = \sum_{kl} S_{ij\;kl} \sigma_{kl}.  @f]
  * Since the stress and strain tensors are symmetric, we denote the
  * elastic moduli and constants  as a 6 by 6 matrix, with the index mapping
  * @f[ (1,2,3,4,5,6) \leftrightarrow (xx,yy,zz,yz,xz,xy). @f]
  * A cubic system has three unique elastic moduli, 
  * @f$C_{11}, C_{12}, C_{44}.@f$ 
  * <p>The elastic constants for a cubic system are 
  * @f[
  * S_{11} = \frac{C_{11}^2-C_{12}^2}{C_{11}^3+2C_{12}^3-3C_{11}C_{12}^2}
  * @f]
  * @f[
  * S_{12} = \frac{C_{12}(C_{12}-C_{11})}{C_{11}^3+2C_{12}^3-3C_{11}C_{12}^2}
  * @f]
  * @f[ S_{44} = \frac{1}{C_{44}}@f]
  * </p>
  * <p>The band offsets are given by</p>
  * @author John Shumway */
typedef blitz::TinyVector<float,3> Vec;
typedef blitz::TinyVector<int,3> IVec;
int main(int argc, char** argv) {
  double scale=1.0;
  if (argc==2) scale=atof(argv[1]);
  bool isInGaAs=true;
  bool isGaAsSb=true;
  H5open();
  std::cout << "getema: calculate effective mass approximation grids"
          << std::endl;
  // Open the relax.h5 file.
  const std::string infile="relaxed.h5";
  std::cout << "Reading structure from " << infile << std::endl;
  SpeciesH5 species(infile);
  CoordsH5 coords(infile);
  NeighborsH5 neighbors(infile);
  SuperCellH5 cell(infile);
  LatticeIndexingH5 lattice(infile);
  Stress stress;
  StressH5 stressH5(stress);
  stressH5.h5Read(infile);
  CommentH5 comment(infile);
  // Identify the atom types.
  int indSi=-1,indGe=-1,indAs=-1,indGa=-1,indIn=-1,indSb=-1;
  for (int i=0; i<species.getNSpecies(); ++i) {
    std::string name=species.getSpeciesName(i);
    if (name=="Si") indSi=i;
    if (name=="Ge") indGe=i;
    if (name=="As") indAs=i;
    if (name=="Ga") indGa=i;
    if (name=="In") indIn=i;
    if (name=="Sb") indSb=i;
  }
  if (indSi>=0 && indGe>=0) {
    isInGaAs=false;
    isGaAsSb=false;
  } else if (indSb>0) {
    isGaAsSb=true;
    isInGaAs=false;
  } else if (indGa<0 || indIn<0) {
    std::cout << "Can't identify material as Si/Ge, InGaAs or GaAsSb, exiting"
              << std::endl; 
    exit(-1);
  }
  std::cout << "Calculating for " << (isInGaAs?"InGaAs":isGaAsSb?"GaAsSb":"Si/Ge") << std::endl; //!!
  Vec delta=lattice.a;
  IVec extent;
  extent[0]=(int)((cell.a1[0]+0.5*delta[0])/delta[0]);
  extent[1]=(int)((cell.a2[1]+0.5*delta[1])/delta[1]);
  extent[2]=(int)((cell.a3[2]+0.5*delta[2])/delta[2]);
  GridFactory gridFactory(extent,delta);
  CompositionGrids composition(gridFactory,species,coords);
  composition.calculate();
  StressGrid stressGrid(gridFactory,stress,coords);
  stressGrid.calculate();
  ElasticConstants* elasticity(0);
  if (isInGaAs) {
    elasticity = new ZincBlendeElasticity(gridFactory,
      ZincBlendeElasticity::GAAS, ZincBlendeElasticity::INAS,indGa,indIn);
  } else if (isGaAsSb) {
    elasticity = new ZincBlendeElasticity(gridFactory,
      ZincBlendeElasticity::GAAS, ZincBlendeElasticity::GASB,indAs,indSb);    
  } else {
    elasticity = new ZincBlendeElasticity(gridFactory,
      ZincBlendeElasticity::SI, ZincBlendeElasticity::GE,indSi,indGe);
  }
  elasticity->calculate(composition);
  StrainGrid *strain=0;
  if (isInGaAs || isGaAsSb) {
    strain = new Strain(gridFactory,coords,composition,neighbors,cell,isInGaAs);
    //strain = new StrainFromStress(gridFactory,*elasticity,stressGrid);
  } else { 
    strain = new StrainFromStress(gridFactory,*elasticity,stressGrid);
  }
  EMAModel* ema(0);
  if (isInGaAs) {
    ema = new ZincBlendeEMA(gridFactory,
              ZincBlendeEMA::GAAS, ZincBlendeEMA::INAS,indGa,indIn);
  } else if (isGaAsSb) {
    ema = new ZincBlendeEMA(gridFactory,
              ZincBlendeEMA::GAAS, ZincBlendeEMA::GASB,indAs,indSb);
  } else {
    ema = new ZincBlendeEMA(gridFactory,
              ZincBlendeEMA::SI, ZincBlendeEMA::GE,indSi,indGe);
  }
  ema->calculate(composition,*strain);

  //Output grids to HDF5 file.
  std::string outfile="emagrids.h5";
  std::cout << "Writting grids to " << outfile << std::endl;
  hid_t fileID = H5Fcreate(outfile.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  //Write composition grids.
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  hid_t grpID = H5Gcreate(fileID,"composition",0,H5P_DEFAULT,H5P_DEFAULT); 
#else
  hid_t grpID = H5Gcreate(fileID,"composition",0);
#endif
  hsize_t dims3[] = {extent[0],extent[1],extent[2]};
  hid_t dspaceID = H5Screate_simple(3,dims3,0);
  for (int i=0; i<species.getNSpecies(); ++i) {
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
    hid_t dsetID = H5Dcreate2(grpID,species.getSpeciesName(i).c_str(),
                              H5T_NATIVE_FLOAT, dspaceID, H5P_DEFAULT,
                              H5P_DEFAULT, H5P_DEFAULT); 
#else
    hid_t dsetID = H5Dcreate(grpID,species.getSpeciesName(i).c_str(),
                             H5T_NATIVE_FLOAT, dspaceID, H5P_DEFAULT); 
#endif
    H5Dwrite(dsetID, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             composition.getGrid(i).getData().data());
    H5Dclose(dsetID);
  }
  H5Sclose(dspaceID);
  H5Gclose(grpID);
  //Write stress grids.
  hsize_t dims4[] = {extent[0],extent[1],extent[2],6};
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  grpID = H5Gcreate(fileID,"stress",0,H5P_DEFAULT,H5P_DEFAULT); 
#else
  grpID = H5Gcreate(fileID,"stress",0);
#endif
  dspaceID = H5Screate_simple(4,dims4,0);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  hid_t dsetID = H5Dcreate2(grpID, "stress", H5T_NATIVE_FLOAT, dspaceID,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 
#else
  hid_t dsetID = H5Dcreate(grpID, "stress", H5T_NATIVE_FLOAT, dspaceID,
                           H5P_DEFAULT); 
#endif
  H5Dwrite(dsetID, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           stressGrid.getGrid().getData().data());
  H5Dclose(dsetID);
  H5Sclose(dspaceID);
  dspaceID = H5Screate_simple(3,dims3,0);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  dsetID = H5Dcreate(grpID, "trace", H5T_NATIVE_FLOAT, dspaceID, H5P_DEFAULT,
                      H5P_DEFAULT, H5P_DEFAULT); 
#else
  dsetID = H5Dcreate(grpID, "trace", H5T_NATIVE_FLOAT, dspaceID, H5P_DEFAULT); 
#endif
  H5Dwrite(dsetID, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           stressGrid.getTraceGrid().getData().data());
  H5Dclose(dsetID);
  H5Sclose(dspaceID);
  dspaceID = H5Screate_simple(3,dims3,0);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  dsetID = H5Dcreate(grpID,"biaxial",H5T_NATIVE_FLOAT,dspaceID,H5P_DEFAULT,
                      H5P_DEFAULT, H5P_DEFAULT); 
#else
  dsetID = H5Dcreate(grpID,"biaxial",H5T_NATIVE_FLOAT,dspaceID,H5P_DEFAULT); 
#endif
  H5Dwrite(dsetID, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           stressGrid.getBiaxialGrid().getData().data());
  H5Dclose(dsetID);
  H5Sclose(dspaceID);
  H5Gclose(grpID);
  //Write strain grids.
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  grpID = H5Gcreate(fileID,"strain",0,H5P_DEFAULT,H5P_DEFAULT); 
#else
  grpID = H5Gcreate(fileID,"strain",0);
#endif
  dspaceID = H5Screate_simple(4,dims4,0);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  dsetID = H5Dcreate(grpID, "strain", H5T_NATIVE_FLOAT, dspaceID, H5P_DEFAULT,
                      H5P_DEFAULT, H5P_DEFAULT); 
#else
  dsetID = H5Dcreate(grpID, "strain", H5T_NATIVE_FLOAT, dspaceID, H5P_DEFAULT); 
#endif
  H5Dwrite(dsetID, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           strain->getGrid().getData().data());
  H5Dclose(dsetID);
  H5Sclose(dspaceID);
  dspaceID = H5Screate_simple(3,dims3,0);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  dsetID = H5Dcreate(grpID, "trace", H5T_NATIVE_FLOAT, dspaceID, H5P_DEFAULT,
                      H5P_DEFAULT, H5P_DEFAULT); 
#else
  dsetID = H5Dcreate(grpID, "trace", H5T_NATIVE_FLOAT, dspaceID, H5P_DEFAULT); 
#endif
  H5Dwrite(dsetID, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           strain->getTraceGrid().getData().data());
  H5Dclose(dsetID);
  H5Sclose(dspaceID);
  dspaceID = H5Screate_simple(3,dims3,0);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  dsetID = H5Dcreate(grpID,"biaxial",H5T_NATIVE_FLOAT,dspaceID,H5P_DEFAULT,
                      H5P_DEFAULT, H5P_DEFAULT); 
#else
  dsetID = H5Dcreate(grpID,"biaxial",H5T_NATIVE_FLOAT,dspaceID,H5P_DEFAULT); 
#endif
  H5Dwrite(dsetID, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           strain->getBiaxialGrid().getData().data());
  H5Dclose(dsetID);
  H5Sclose(dspaceID);
  H5Gclose(grpID);
  //Write band offset grids.
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  grpID = H5Gcreate(fileID,"boffset",0,H5P_DEFAULT,H5P_DEFAULT); 
#else
  grpID = H5Gcreate(fileID,"boffset",0);
#endif
  dspaceID = H5Screate_simple(3,dims3,0);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  dsetID = H5Dcreate2(grpID, "ve", H5T_NATIVE_FLOAT, dspaceID, H5P_DEFAULT,
                      H5P_DEFAULT, H5P_DEFAULT); 
#else
  dsetID = H5Dcreate(grpID, "ve", H5T_NATIVE_FLOAT, dspaceID, H5P_DEFAULT); 
#endif
  H5Dwrite(dsetID, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           ema->getVeGrid().getData().data());
  H5Dclose(dsetID);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  dsetID = H5Dcreate2(grpID, "vh", H5T_NATIVE_FLOAT, dspaceID, H5P_DEFAULT,
                      H5P_DEFAULT, H5P_DEFAULT); 
#else
  dsetID = H5Dcreate(grpID, "vh", H5T_NATIVE_FLOAT, dspaceID, H5P_DEFAULT); 
#endif
  H5Dwrite(dsetID, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           ema->getVhGrid().getData().data());
  H5Dclose(dsetID);
  delta*=scale;
  hsize_t dims1[] = {1};
  dspaceID = H5Screate_simple(1,dims1,0);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  dsetID = H5Dcreate2(grpID, "a", H5T_NATIVE_FLOAT, dspaceID, H5P_DEFAULT,
                      H5P_DEFAULT, H5P_DEFAULT); 
#else
  dsetID = H5Dcreate(grpID, "a", H5T_NATIVE_FLOAT, dspaceID, H5P_DEFAULT); 
#endif
  H5Dwrite(dsetID, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           delta.data());
  H5Dclose(dsetID);
  H5Sclose(dspaceID);
  H5Gclose(grpID);
  //Write the descriptive comment.
  comment.h5Write(outfile);
  H5close();
}

