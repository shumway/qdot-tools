// $Id: SuperCellH5.cc,v 1.5 2007/03/14 19:47:50 jshumwa Exp $
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
#include "SuperCellH5.h"
#include <hdf5.h>
#include <iostream>
#include <blitz/tinyvec-et.h>

SuperCellH5::SuperCellH5() {
  a1[0]=0; a1[1]=0; a1[2]=0;
  a2[0]=0; a2[1]=0; a2[2]=0;
  a3[0]=0; a3[1]=0; a3[2]=0;
}

SuperCellH5::SuperCellH5(const std::string& filename) {
  h5Read(filename);
}

void SuperCellH5::h5Read(const std::string& filename) {
  std::cout << "Reading supercell info from " << filename << std::endl;
  H5open();
  hid_t fileID = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  hid_t grpID = H5Gopen(fileID,"superCell");
  hid_t attrID = H5Aopen_name(grpID,"a1");
#ifdef ENABLE_FLOAT
  H5Aread(attrID,H5T_NATIVE_FLOAT,&a1);
#else
  H5Aread(attrID,H5T_NATIVE_DOUBLE,&a1);
#endif
  H5Aclose(attrID);
  attrID = H5Aopen_name(grpID,"a2");
#ifdef ENABLE_FLOAT
  H5Aread(attrID,H5T_NATIVE_FLOAT,&a2);
#else
  H5Aread(attrID,H5T_NATIVE_DOUBLE,&a2);
#endif
  H5Aclose(attrID);
  attrID = H5Aopen_name(grpID,"a3");
#ifdef ENABLE_FLOAT
  H5Aread(attrID,H5T_NATIVE_FLOAT,&a3);
#else
  H5Aread(attrID,H5T_NATIVE_DOUBLE,&a3);
#endif
  H5Aclose(attrID);
  H5Gclose(grpID);
  H5Fclose(fileID);
  H5close();
  computeRecipricalVectors();
}

void SuperCellH5::computeRecipricalVectors() {
  double invol=1/(a1[0]*(a2[1]*a3[2]-a2[2]*a3[1]) +
                  a1[2]*(a2[2]*a3[0]-a2[0]*a3[2]) +
                  a1[3]*(a2[0]*a3[1]-a2[1]*a3[0]));
  b1[0]=invol*(a2[1]*a3[2]-a2[2]*a3[1]);
  b1[1]=invol*(a2[2]*a3[0]-a2[0]*a3[2]);
  b1[2]=invol*(a2[0]*a3[1]-a2[1]*a3[0]);
  b2[0]=invol*(a3[1]*a1[2]-a3[2]*a1[1]);
  b2[1]=invol*(a3[2]*a1[0]-a3[0]*a1[2]);
  b2[2]=invol*(a3[0]*a1[1]-a3[1]*a1[0]);
  b3[0]=invol*(a1[1]*a2[2]-a1[2]*a2[1]);
  b3[1]=invol*(a1[2]*a2[0]-a1[0]*a2[2]);
  b3[2]=invol*(a1[0]*a2[1]-a1[1]*a2[0]);
}

void SuperCellH5::pbc(Vec3& v) const {
  const double rcut2=25.;
  if (dot(v,v)>rcut2) {
    v-=a1*floor(dot(b1,v)+0.5);
    v-=a2*floor(dot(b2,v)+0.5);
    v-=a3*floor(dot(b3,v)+0.5);
  }
}

void SuperCellH5::h5Write(const std::string& filename, const int mode) const {
  std::cout << "Writing supercell info to " << filename << std::endl;
  H5open();
  hid_t fileID;
  if (mode==NEW) {
    fileID = H5Fcreate(filename.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  } else {
    fileID = H5Fopen(filename.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
  }
  if (mode!=REWRITE) {
  hid_t grpID = H5Gcreate(fileID,"superCell",0);
  const hsize_t dims[] = {3};
  hid_t aspaceID = H5Screate_simple(1,dims,0);
  hid_t attrID = H5Acreate(grpID,"a1",H5T_NATIVE_FLOAT,aspaceID,H5P_DEFAULT);
#ifdef ENABLE_FLOAT
  H5Awrite(attrID,H5T_NATIVE_FLOAT,&a1);
#else
  H5Awrite(attrID,H5T_NATIVE_DOUBLE,&a1);
#endif
  H5Aclose(attrID);
  attrID = H5Acreate(grpID,"a2",H5T_NATIVE_FLOAT,aspaceID,H5P_DEFAULT);
#ifdef ENABLE_FLOAT
  H5Awrite(attrID,H5T_NATIVE_FLOAT,&a2);
#else
  H5Awrite(attrID,H5T_NATIVE_DOUBLE,&a2);
#endif
  H5Aclose(attrID);
  attrID = H5Acreate(grpID,"a3",H5T_NATIVE_FLOAT,aspaceID,H5P_DEFAULT);
#ifdef ENABLE_FLOAT
  H5Awrite(attrID,H5T_NATIVE_FLOAT,&a3);
#else
  H5Awrite(attrID,H5T_NATIVE_DOUBLE,&a3);
#endif
  H5Aclose(attrID);
  H5Sclose(aspaceID);
  H5Gclose(grpID);
  } else {
    hid_t grpID = H5Gopen(fileID,"superCell");
    hid_t attrID = H5Aopen_name(grpID,"a1");
#ifdef ENABLE_FLOAT
    H5Awrite(attrID,H5T_NATIVE_FLOAT,&a1);
#else
    H5Awrite(attrID,H5T_NATIVE_DOUBLE,&a1);
#endif
    attrID = H5Aopen_name(grpID,"a2");
#ifdef ENABLE_FLOAT
    H5Awrite(attrID,H5T_NATIVE_FLOAT,&a2);
#else
    H5Awrite(attrID,H5T_NATIVE_DOUBLE,&a2);
#endif
    attrID = H5Aopen_name(grpID,"a3");
#ifdef ENABLE_FLOAT
    H5Awrite(attrID,H5T_NATIVE_FLOAT,&a3);
#else
    H5Awrite(attrID,H5T_NATIVE_DOUBLE,&a3);
#endif
    H5Gclose(grpID);
  }
  H5Fclose(fileID);
  H5close();
}
