// $Id: CoordsH5.cc,v 1.4 2004/10/08 13:20:22 jshumwa Exp $
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
#include "CoordsH5.h"
#include "hdf5.h"
#include <iostream>

CoordsH5::CoordsH5 (const int natom) 
  : coords(natom) {
}

void CoordsH5::h5Write(const std::string& filename, const int mode) const {
  std::cout << "Writing coordinates to " << filename << std::endl;
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
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  hid_t grpID = H5Gcreate2(fileID, "coords", H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);
#else
  hid_t grpID = H5Gcreate(fileID, "coords", 0);
#endif
  // Write out the coordinates.
  int natoms = coords.size();
  hsize_t dims[] = {natoms,3};
  hid_t dspaceID = H5Screate_simple(2,dims,0);
  hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
  dims[0] = (natoms<50000) ? natoms : 50000;
  H5Pset_chunk(plist,2,dims);
  H5Pset_deflate(plist,1);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  hid_t dsetID = H5Dcreate2(grpID, "coord",
                         H5T_NATIVE_FLOAT, dspaceID, H5P_DEFAULT,
                         plist, H5P_DEFAULT);
#else
  hid_t dsetID = H5Dcreate(grpID,"coord",H5T_NATIVE_FLOAT,dspaceID,plist);
#endif
  H5Pclose(plist);
  plist = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_buffer(plist,100000000,0,0);
  H5Dwrite(dsetID,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,plist,&coords[0]);
  H5Dclose(dsetID);
  H5Pclose(plist);
  H5Sclose(dspaceID);
  // Write the natoms attribute.
  hid_t aspaceID = H5Screate(H5S_SCALAR);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  hid_t attrID = H5Acreate(grpID,"nAtom",H5T_NATIVE_INT,aspaceID,H5P_DEFAULT,
                           H5P_DEFAULT);
#else
  hid_t attrID = H5Acreate(grpID,"nAtom",H5T_NATIVE_INT,aspaceID,H5P_DEFAULT);
#endif
  H5Awrite(attrID,H5T_NATIVE_INT,&natoms);
  H5Aclose(attrID);
  H5Sclose(aspaceID);
  // Close the file
  H5Gclose(grpID);
  H5Fclose(fileID);
  H5close();
}
