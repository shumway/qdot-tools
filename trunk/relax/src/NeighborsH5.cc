// $Id: NeighborsH5.cc,v 1.2 2004/08/14 01:42:07 jshumwa Exp $
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
#include "NeighborsH5.h"
#include <hdf5.h>
#include <iostream>

NeighborsH5::NeighborsH5 (const int nmax, const int natoms) 
  : nmax(nmax), nn(nmax*natoms) {
}

NeighborsH5::NeighborsH5 (const std::string& filename) {
  h5Read(filename);
}

void NeighborsH5::h5Read(const std::string& filename) {
  std::cout << "Reading neighbor info from " << filename << std::endl;
  H5open();
  hid_t fileID = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  //Read the neighbor list.
  {
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
    hid_t dsetID = H5Dopen2(fileID,"neighbors/neighbor",H5P_DEFAULT);
#else
    hid_t dsetID = H5Dopen(fileID,"neighbors/neighbor");
#endif
    hid_t dspaceID = H5Dget_space(dsetID);
    hsize_t dims[2];
    H5Sget_simple_extent_dims(dspaceID,dims,NULL);
    int natom=dims[0]; nmax=dims[1];
    std::cout << "natom " << natom << std::endl;
    std::cout << "nmax " << nmax << std::endl;
    nn.resize(natom*nmax);
    H5Dread(dsetID,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,&nn[0]);
    H5Sclose(dspaceID);
    H5Dclose(dsetID);
  }
  H5Fclose(fileID);
  H5close();
}

void NeighborsH5::h5Write(const std::string& filename, const int mode) const {
  std::cout << "Writing neighbor info to " << filename << std::endl;
  H5open();
  hid_t fileID=0;
  if (mode==NEW) {
    fileID = H5Fcreate(filename.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  } else {
    fileID = H5Fopen(filename.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
  }
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  hid_t grpID = H5Gcreate2(fileID,"neighbors",0,H5P_DEFAULT,H5P_DEFAULT);
#else
  hid_t grpID = H5Gcreate(fileID,"neighbors",0);
#endif
  // Write out the type index array.
  int natoms = nn.size()/nmax;
  hsize_t dims[] = {natoms,nmax};
  hid_t dspaceID = H5Screate_simple(2,dims,0);
  hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
  dims[0] = (natoms<50000) ? natoms : 50000;
  H5Pset_chunk(plist,2,dims);
  H5Pset_deflate(plist,1);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  hid_t dsetID = H5Dcreate2(grpID,"neighbor",H5T_NATIVE_INT,dspaceID,plist,
                            H5P_DEFAULT,H5P_DEFAULT);
#else
  hid_t dsetID = H5Dcreate(grpID,"neighbor",H5T_NATIVE_INT,dspaceID,plist);
#endif
  H5Pclose(plist);
  plist = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_buffer(plist,10000000,0,0);
  H5Dwrite(dsetID,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,plist,&nn[0]);
  H5Dclose(dsetID);
  H5Pclose(plist);
  H5Sclose(dspaceID);
  // Close the file
  H5Gclose(grpID);
  H5Fclose(fileID);
  H5close();
}
