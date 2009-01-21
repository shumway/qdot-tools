// $Id: StressH5.cc,v 1.6 2007/03/14 19:47:50 jshumwa Exp $
/*
    Copyright (C) 2004,2008 John B. Shumway, Jr.

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
#include "StressH5.h"
#include "Stress.h"
#include <hdf5.h>
#include <iostream>

void StressH5::h5Read(const std::string& filename) {
  std::cout << "Reading coordinates from " << filename << std::endl;
  H5open();
  hid_t fileID = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  //Read the stress tensors.
  hid_t dsetID = H5Dopen(fileID,"stress/stress");
  hid_t dspaceID = H5Dget_space(dsetID);
  hsize_t dims[2];
  H5Sget_simple_extent_dims(dspaceID,dims,NULL);
  H5Sclose(dspaceID);
  stress.stress.resize(dims[0]);
#ifdef ENABLE_FLOAT
  H5Dread(dsetID,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,
    stress.stress.data());
#else
  H5Dread(dsetID,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,
    stress.stress.data());
#endif
  H5Fclose(fileID);
  H5close();
}

void StressH5::h5Write(const std::string& filename, const int mode) const {
  //H5::Exception::dontPrint();
  std::cout << "Writing stresses to " << filename << std::endl;
  H5open();
  hid_t fileID;
  if (mode==NEW) {
    fileID = H5Fcreate(filename.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  } else {
    fileID = H5Fopen(filename.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
  }
  hid_t dsetID; 
  hid_t grpID; 
  hid_t plist; 
  const int natoms = stress.stress.size();
  if (mode!=REWRITE) {
    grpID = H5Gcreate(fileID,"stress",0);
    // Write out the stress.
    hsize_t dims[] = {natoms,6};
    hid_t dspaceID = H5Screate_simple(2,dims,0);
    hid_t plist = H5Pcreate(H5P_DATASET_CREATE); 
    dims[0] = (natoms<50000) ? natoms : 50000;
    H5Pset_chunk(plist,2,dims);
    H5Pset_deflate(plist,1);
    hid_t dsetID = H5Dcreate(grpID,"stress",H5T_NATIVE_FLOAT,dspaceID,plist);
    H5Pclose(plist);
    plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_buffer(plist,100000000,0,0); 
#ifdef ENABLE_FLOAT
    H5Dwrite(dsetID,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,plist,
             stress.stress.data());
#else
    H5Dwrite(dsetID,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,plist,
             stress.stress.data());
#endif
    H5Dclose(dsetID);
    H5Pclose(plist);
    H5Sclose(dspaceID);
  } else {
    dsetID = H5Dopen(fileID,"stress/stress");
#ifdef ENABLE_FLOAT
    H5Dwrite(dsetID,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,0,stress.stress.data());
#else
    H5Dwrite(dsetID,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,0,stress.stress.data());
#endif
    H5Dclose(dsetID);
  }
  if (mode!=REWRITE) {
    hid_t aspaceID = H5Screate(H5S_SCALAR);
    hid_t attrID = H5Acreate(grpID,"nAtom",H5T_NATIVE_INT,aspaceID,H5P_DEFAULT);
    H5Awrite(attrID,H5T_NATIVE_INT,&natoms);
    H5Aclose(attrID);
    H5Sclose(aspaceID);
    H5Gclose(grpID);
  }
  // Close the file.
  H5Fclose(fileID);
  H5close();
}
