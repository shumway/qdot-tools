// $Id: NeighborsH5.cc,v 1.3 2004/05/18 22:33:59 matth2 Exp $
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
#include "NeighborsH5.h"
#include "hdf5.h"
#include <iostream>

NeighborsH5::NeighborsH5 (const int nmax, const int natoms) 
  : nmax(nmax), nn(nmax*natoms) {
}

void NeighborsH5::h5Write(const std::string& filename, const int mode) const {
  std::cout << "Writing neighbor info to " << filename << std::endl;
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
  hid_t grpID = H5Gcreate(fileID,"neighbors",0);
  // Write out the type index array.
  int natoms = nn.size()/nmax;
  hsize_t dims[] = {natoms,nmax};
  hid_t dspaceID = H5Screate_simple(2,dims,0);
  hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
  dims[0] = (natoms<50000) ? natoms : 50000;
  H5Pset_chunk(plist,2,dims);
  H5Pset_deflate(plist,1);
  hid_t dsetID = H5Dcreate(grpID,"neighbor",H5T_NATIVE_INT,dspaceID,plist);
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
