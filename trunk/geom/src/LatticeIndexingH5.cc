// $Id: LatticeIndexingH5.cc,v 1.5 2006/07/01 19:18:52 jshumwa Exp $
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
#include "LatticeIndexingH5.h"
#include "hdf5.h"
#include <iostream>

LatticeIndexingH5::LatticeIndexingH5(const double a, const Vec3& a1, 
  const Vec3& a2, const Vec3& a3, const std::valarray<Vec3>& offset,
  const int natoms, const IVec3 n1, const IVec3 n2, const IVec3 n3)
  : a(a),a1(a1),a2(a2),a3(a3),nAtomPerPrimCell(offset.size()),offset(offset),
  index(natoms),n1(n1),n2(n2),n3(n3) {
}

int LatticeIndexingH5::getNAtomPerPrimCell() const {
  return nAtomPerPrimCell;
}

const Vec3& LatticeIndexingH5::getPrimAtomOffset(const int i) const {
  return offset[i];
}

const Vec3 LatticeIndexingH5::getIdealPosition(const int i) const {
  indices ind=index[i];
  Vec3 off(offset[ind.ipos]);
  Vec3 p;
  p.x = (ind.i*a1.x+ind.j*a2.x+ind.k*a3.x+off.x)*a;
  p.y = (ind.i*a1.y+ind.j*a2.y+ind.k*a3.y+off.y)*a;
  p.z = (ind.i*a1.z+ind.j*a2.z+ind.k*a3.z+off.z)*a;
  return p;
}

const Vec3 LatticeIndexingH5::getIdealPositionOverA(const int i) const {
  indices ind=index[i];
  Vec3 off(offset[ind.ipos]);
  Vec3 p;
  p.x = (ind.i*a1.x+ind.j*a2.x+ind.k*a3.x+off.x);
  p.y = (ind.i*a1.y+ind.j*a2.y+ind.k*a3.y+off.y);
  p.z = (ind.i*a1.z+ind.j*a2.z+ind.k*a3.z+off.z);
  return p;
}

int LatticeIndexingH5::getIPos(const int i) const {
  return index[i].ipos;
}

void LatticeIndexingH5::h5Write(const std::string& filename, 
                                const int mode) const {
  std::cout << "Writing lattice indexing info to " << filename << std::endl;
  H5open();
  hid_t fileID;
  switch(mode) {
  case (NEW):
    fileID = H5Fcreate(filename.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
    break;
  case (APPEND):
    fileID = H5Fopen(filename.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
    break;
  }
  hid_t grpID = H5Gcreate(fileID,"latticeIndexing",0);
  // Write out the lattice constant in atomic units.
  {hid_t aspaceID = H5Screate(H5S_SCALAR);
  hid_t attrID = H5Acreate(grpID,"a",H5T_NATIVE_FLOAT,aspaceID,H5P_DEFAULT);
  H5Awrite(attrID,H5T_NATIVE_DOUBLE,(double*)&a);
  H5Aclose(attrID);
  H5Sclose(aspaceID);}
  // Write out the primative cell lattice vectors in units of a.
  {const hsize_t dims[] = {3};
  hid_t aspaceID = H5Screate_simple(1,dims,0);
  hid_t attrID = H5Acreate(grpID,"a1",H5T_NATIVE_FLOAT,aspaceID,H5P_DEFAULT);
  H5Awrite(attrID,H5T_NATIVE_FLOAT,(Vec3*)&a1);
  H5Aclose(attrID);
  attrID = H5Acreate(grpID,"a2",H5T_NATIVE_FLOAT,aspaceID,H5P_DEFAULT);
  H5Awrite(attrID,H5T_NATIVE_FLOAT,(Vec3*)&a2);
  H5Aclose(attrID);
  attrID = H5Acreate(grpID,"a3",H5T_NATIVE_FLOAT,aspaceID,H5P_DEFAULT);
  H5Awrite(attrID,H5T_NATIVE_FLOAT,(Vec3*)&a3);
  H5Aclose(attrID);
  H5Sclose(aspaceID);}
  // Write out the number of atoms in the primative cell.
  {hid_t aspaceID = H5Screate(H5S_SCALAR);
  hid_t attrID = H5Acreate(grpID,"nAtomsPerPrimCell",H5T_NATIVE_INT,
                           aspaceID,H5P_DEFAULT);
  H5Awrite(attrID,H5T_NATIVE_INT,(int*)&nAtomPerPrimCell);
  H5Aclose(attrID);
  H5Sclose(aspaceID);}
  // Write out the positions of the atoms in the primative cell.
  { hsize_t dims[] = {nAtomPerPrimCell,3};
  hid_t aspaceID = H5Screate_simple(2,dims,0);
  hid_t attrID = H5Acreate(grpID,"primativeAtomOffsets",
                           H5T_NATIVE_FLOAT,aspaceID,H5P_DEFAULT);
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
  hid_t dsetID = H5Dcreate(grpID,"index",H5T_NATIVE_INT,dspaceID,plist);
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
  hid_t attrID = H5Acreate(grpID,"n1",H5T_NATIVE_INT,aspaceID,H5P_DEFAULT);
  H5Awrite(attrID,H5T_NATIVE_INT,(IVec3*)&n1);
  H5Aclose(attrID);
  attrID = H5Acreate(grpID,"n2",H5T_NATIVE_INT,aspaceID,H5P_DEFAULT);
  H5Awrite(attrID,H5T_NATIVE_INT,(IVec3*)&n2);
  H5Aclose(attrID);
  attrID = H5Acreate(grpID,"n3",H5T_NATIVE_INT,aspaceID,H5P_DEFAULT);
  H5Awrite(attrID,H5T_NATIVE_INT,(IVec3*)&n3);
  H5Aclose(attrID);
  H5Sclose(aspaceID);}
  // Close the group.
  H5Gclose(grpID);
  H5Fclose(fileID);
  H5close();
}
