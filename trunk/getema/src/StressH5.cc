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
  std::cout << "Reading stress tensors from " << filename << std::endl;
  H5open();
  hid_t fileID = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  //Read the stress tensors.
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  hid_t dsetID = H5Dopen2(fileID,"stress/stress",H5P_DEFAULT);
#else
  hid_t dsetID = H5Dopen(fileID,"stress/stress");
#endif
  hid_t dspaceID = H5Dget_space(dsetID);
  hsize_t dims[2];
  H5Sget_simple_extent_dims(dspaceID,dims,NULL);
  H5Sclose(dspaceID);
  stress.stress.resize(dims[0]);
  H5Dread(dsetID,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,
    stress.stress.data());
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
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
    grpID = H5Gcreate2(fileID,"stress",0,H5P_DEFAULT,H5P_DEFAULT); 
#else
    grpID = H5Gcreate(fileID,"stress",0);
#endif
    // Write out the stress.
    hsize_t dims[] = {natoms,6};
    hid_t dspaceID = H5Screate_simple(2,dims,0);
    hid_t plist = H5Pcreate(H5P_DATASET_CREATE); 
    dims[0] = (natoms<50000) ? natoms : 50000;
    H5Pset_chunk(plist,2,dims);
    H5Pset_deflate(plist,1);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
    hid_t dsetID = H5Dcreate2(grpID,"stress",H5T_NATIVE_FLOAT,dspaceID,plist,
                              H5P_DEFAULT,H5P_DEFAULT);
#else
    hid_t dsetID = H5Dcreate(grpID,"stress",H5T_NATIVE_FLOAT,dspaceID,plist);
#endif
    H5Pclose(plist);
    plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_buffer(plist,100000000,0,0); 
    H5Dwrite(dsetID,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,plist,
             stress.stress.data());
    H5Dclose(dsetID);
    H5Pclose(plist);
    H5Sclose(dspaceID);
  } else {
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
    dsetID = H5Dopen2(fileID,"stress/stress",H5P_DEFAULT);
#else
    dsetID = H5Dopen(fileID,"stress/stress");
#endif
    H5Dwrite(dsetID,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,0,stress.stress.data());
    H5Dclose(dsetID);
  }
  if (mode!=REWRITE) {
    hid_t aspaceID = H5Screate(H5S_SCALAR);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
    hid_t attrID = H5Acreate2(grpID,"nAtom",H5T_NATIVE_INT,aspaceID,H5P_DEFAULT,
                              H5P_DEFAULT);
#else
    hid_t attrID = H5Acreate(grpID,"nAtom",H5T_NATIVE_INT,aspaceID,H5P_DEFAULT);
#endif
    H5Awrite(attrID,H5T_NATIVE_INT,&natoms);
    H5Aclose(attrID);
    H5Sclose(aspaceID);
    H5Gclose(grpID);
  }
  // Close the file.
  H5Fclose(fileID);
  H5close();
}
