// $Id: LatticeIndexingH5.cc,v 1.5 2007/03/14 19:47:50 jshumwa Exp $
/*
    Copyright (C) 2004, 2009 John B. Shumway, Jr.

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
#include "LatticeIndexingH5.h"
#include <hdf5.h>
#include <iostream>

LatticeIndexingH5::LatticeIndexingH5(const double a, const Vec3& a1, 
  const Vec3& a2, const Vec3& a3, const std::valarray<Vec3>& offset,
  const int natoms, const IVec3 n1, const IVec3 n2, const IVec3 n3)
  : index(natoms), 
  a(a),a1(a1),a2(a2),a3(a3),nAtomPerPrimCell(offset.size()),offset(offset),
  n1(n1),n2(n2),n3(n3) {
}

LatticeIndexingH5::LatticeIndexingH5(const std::string& filename) {
  h5Read(filename);
}

void LatticeIndexingH5::h5Read(const std::string& filename) {
  std::cout << "Reading lattice indexing from " << filename << std::endl;
  H5open();
  hid_t fileID = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  hid_t grpID = H5Gopen2(fileID,"latticeIndexing",H5P_DEFAULT);
#else
  hid_t grpID = H5Gopen(fileID,"latticeIndexing");
#endif
  { // Read the lattice constant in atomic units.
    hid_t attrID = H5Aopen_name(grpID,"a");
    H5Aread(attrID,H5T_NATIVE_DOUBLE,&a);
    H5Aclose(attrID);
  }
  { // Read the lattice index array.
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
    hid_t dsetID = H5Dopen2(grpID,"index",H5P_DEFAULT);
#else
    hid_t dsetID = H5Dopen(grpID,"index");
#endif
    hid_t dspaceID = H5Dget_space(dsetID);
    hsize_t dims[2];
    H5Sget_simple_extent_dims(dspaceID,dims,NULL);
    H5Sclose(dspaceID);
    index.resize(dims[0]);
    H5Dread(dsetID,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,&index[0]);
    H5Dclose(dsetID);
  }
  H5Gclose(grpID);
  H5Fclose(fileID);
  H5close();
}


int LatticeIndexingH5::getNAtomPerPrimCell() const {
  return nAtomPerPrimCell;
}

const LatticeIndexingH5::Vec3&
LatticeIndexingH5::getPrimAtomOffset(const int i) const {
  return offset[i];
}

const LatticeIndexingH5::Vec3
 LatticeIndexingH5::getIdealPosition(const int i) const {
  indices ind=index[i];
  Vec3 off(offset[ind.ipos]);
  Vec3 p;
  p[0] = (ind.i*a1[0]+ind.j*a2[0]+ind.k*a3[0]+off[0])*a;
  p[1] = (ind.i*a1[1]+ind.j*a2[1]+ind.k*a3[1]+off[1])*a;
  p[2] = (ind.i*a1[2]+ind.j*a2[2]+ind.k*a3[1]+off[2])*a;
  return p;
}

const LatticeIndexingH5::Vec3 
LatticeIndexingH5::getIdealPositionOverA(const int i) const {
  indices ind=index[i];
  Vec3 off(offset[ind.ipos]);
  Vec3 p;
  p[0] = (ind.i*a1[0]+ind.j*a2[0]+ind.k*a3[0]+off[0]);
  p[1] = (ind.i*a1[1]+ind.j*a2[1]+ind.k*a3[1]+off[1]);
  p[2] = (ind.i*a1[2]+ind.j*a2[2]+ind.k*a3[1]+off[2]);
  return p;
}

int LatticeIndexingH5::getIPos(const int i) const {
  return index[i].ipos;
}

void LatticeIndexingH5::h5Write(const std::string& filename, 
                                const int mode) const {
  std::cout << "Writting lattice indexing info to " << filename << std::endl;
  H5open();
  hid_t fileID=0;
  if (mode==NEW) {
    fileID = H5Fcreate(filename.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  } else {
    fileID = H5Fopen(filename.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
  }
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  hid_t grpID = H5Gcreate2(fileID,"latticeIndexing",0,H5P_DEFAULT,H5P_DEFAULT);
#else
  hid_t grpID = H5Gcreate(fileID,"latticeIndexing",0);
#endif
  // Write out the lattice constant in atomic units.
  {hid_t aspaceID = H5Screate(H5S_SCALAR);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  hid_t attrID = H5Acreate2(grpID,"a",H5T_NATIVE_FLOAT,aspaceID,H5P_DEFAULT,
                            H5P_DEFAULT);
#else
  hid_t attrID = H5Acreate(grpID,"a",H5T_NATIVE_FLOAT,aspaceID,H5P_DEFAULT);
#endif
  H5Awrite(attrID,H5T_NATIVE_DOUBLE,(double*)&a);
  H5Aclose(attrID);
  H5Sclose(aspaceID);}
  // Write out the primative cell lattice vectors in units of a.
  {const hsize_t dims[] = {3};
  hid_t aspaceID = H5Screate_simple(1,dims,0);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  hid_t attrID = H5Acreate2(grpID,"a1",H5T_NATIVE_FLOAT,aspaceID,H5P_DEFAULT,
                            H5P_DEFAULT);
#else
  hid_t attrID = H5Acreate(grpID,"a1",H5T_NATIVE_FLOAT,aspaceID,H5P_DEFAULT);
#endif
  H5Awrite(attrID,H5T_NATIVE_FLOAT,(Vec3*)&a1);
  H5Aclose(attrID);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  attrID = H5Acreate2(grpID,"a2",H5T_NATIVE_FLOAT,aspaceID,H5P_DEFAULT,
                      H5P_DEFAULT);
#else
  attrID = H5Acreate(grpID,"a2",H5T_NATIVE_FLOAT,aspaceID,H5P_DEFAULT);
#endif
  H5Awrite(attrID,H5T_NATIVE_FLOAT,(Vec3*)&a2);
  H5Aclose(attrID);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  attrID = H5Acreate2(grpID,"a3",H5T_NATIVE_FLOAT,aspaceID,H5P_DEFAULT,
                      H5P_DEFAULT);
#else
  attrID = H5Acreate(grpID,"a3",H5T_NATIVE_FLOAT,aspaceID,H5P_DEFAULT);
#endif
  H5Awrite(attrID,H5T_NATIVE_FLOAT,(Vec3*)&a3);
  H5Aclose(attrID);
  H5Sclose(aspaceID);}
  // Write out the number of atoms in the primative cell.
  {hid_t aspaceID = H5Screate(H5S_SCALAR);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  hid_t attrID = H5Acreate2(grpID,"nAtomsPerPrimCell",H5T_NATIVE_INT,
                            aspaceID,H5P_DEFAULT,H5P_DEFAULT);
#else
  hid_t attrID = H5Acreate(grpID,"nAtomsPerPrimCell",H5T_NATIVE_INT,
                           aspaceID,H5P_DEFAULT);
#endif
  H5Awrite(attrID,H5T_NATIVE_INT,(int*)&nAtomPerPrimCell);
  H5Aclose(attrID);
  H5Sclose(aspaceID);}
  // Write out the positions of the atoms in the primative cell.
  { hsize_t dims[] = {nAtomPerPrimCell,3};
  hid_t aspaceID = H5Screate_simple(2,dims,0);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  hid_t attrID = H5Acreate2(grpID,"primativeAtomOffsets",
                            H5T_NATIVE_FLOAT,aspaceID,H5P_DEFAULT,H5P_DEFAULT);
#else
  hid_t attrID = H5Acreate(grpID,"primativeAtomOffsets",
                           H5T_NATIVE_FLOAT,aspaceID,H5P_DEFAULT);
#endif
  H5Awrite(attrID,H5T_NATIVE_FLOAT,(Vec3*)&offset[0]);
  H5Aclose(attrID);
  H5Sclose(aspaceID);}
  // Write out the type index array.
  {int natoms = index.size();
  hsize_t dims[] = {natoms,4};
  hid_t dspaceID = H5Screate_simple(2,dims,0);
  hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
  dims[0] = (natoms<50000) ? natoms : 50000;
  H5Pset_chunk(plist,2,dims);
  H5Pset_deflate(plist,1);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  hid_t dsetID = H5Dcreate2(grpID,"index",H5T_NATIVE_INT,dspaceID,plist,
                            H5P_DEFAULT,H5P_DEFAULT);
#else
  hid_t dsetID = H5Dcreate(grpID,"index",H5T_NATIVE_INT,dspaceID,plist);
#endif
  H5Pclose(plist);
  plist = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_buffer(plist,10000000,0,0);
  H5Dwrite(dsetID,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,plist,&index[0]);
  H5Dclose(dsetID);
  H5Pclose(plist);
  H5Sclose(dspaceID);}
  // Write supercell vectors.
  {const hsize_t dims[] = {3};
  hid_t aspaceID = H5Screate_simple(1,dims,0);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  hid_t attrID = H5Acreate2(grpID,"n1",H5T_NATIVE_INT,aspaceID,H5P_DEFAULT,
                            H5P_DEFAULT);
#else
  hid_t attrID = H5Acreate(grpID,"n1",H5T_NATIVE_INT,aspaceID,H5P_DEFAULT);
#endif
  H5Awrite(attrID,H5T_NATIVE_INT,(IVec3*)&n1);
  H5Aclose(attrID);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  attrID = H5Acreate2(grpID,"n2",H5T_NATIVE_INT,aspaceID,H5P_DEFAULT,
                      H5P_DEFAULT);
#else
  attrID = H5Acreate(grpID,"n2",H5T_NATIVE_INT,aspaceID,H5P_DEFAULT);
#endif
  H5Awrite(attrID,H5T_NATIVE_INT,(IVec3*)&n2);
  H5Aclose(attrID);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  attrID = H5Acreate2(grpID,"n3",H5T_NATIVE_INT,aspaceID,H5P_DEFAULT,
                      H5P_DEFAULT);
#else
  attrID = H5Acreate(grpID,"n3",H5T_NATIVE_INT,aspaceID,H5P_DEFAULT);
#endif
  H5Awrite(attrID,H5T_NATIVE_INT,(IVec3*)&n3);
  H5Aclose(attrID);
  H5Sclose(aspaceID);}
  // Close the group.
  H5Gclose(grpID);
  H5Fclose(fileID);
  H5close();
}

