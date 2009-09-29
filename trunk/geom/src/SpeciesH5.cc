// $Id: SpeciesH5.cc,v 1.3 2004/05/18 22:33:59 matth2 Exp $
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
#include "SpeciesH5.h"
#include "hdf5.h"
#include <iostream>

SpeciesH5::SpeciesH5 (const int nspecies, const int natoms) 
  : name(nspecies), species(natoms) {
}

void SpeciesH5::h5Write(const std::string& filename, const int mode) const {
  std::cout << "Writing species info to " << filename << std::endl;
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
  hid_t grpID = H5Gcreate2(fileID, "species", H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);
#else
  hid_t grpID = H5Gcreate(fileID,"species",0);
#endif
  // Write out the type index array.
  int natoms = species.size();
  hsize_t dims[] = {natoms};
  hid_t dspaceID = H5Screate_simple(1,dims,0);
  hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
  dims[0] = (natoms<50000) ? natoms : 50000;
  H5Pset_chunk(plist,1,dims);
  H5Pset_deflate(plist,9);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  hid_t dsetID = H5Dcreate2(grpID, "species",
                         H5T_NATIVE_INT, dspaceID, H5P_DEFAULT,
                         plist, H5P_DEFAULT);
#else
  hid_t dsetID = H5Dcreate(grpID,"species",H5T_NATIVE_INT,dspaceID,plist);
#endif
  H5Pclose(plist);
  plist = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_buffer(plist,10000000,0,0);
  H5Dwrite(dsetID,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,plist,&species[0]);
  H5Dclose(dsetID);
  H5Pclose(plist);
  H5Sclose(dspaceID);
  // Write out the type name array (need to create an char[8] type).
  int nspecies = name.size();
  hid_t char8 = H5Tcopy(H5T_C_S1); H5Tset_size(char8,8);
  dims[0] = nspecies;
  hid_t aspaceID = H5Screate_simple(1,dims,0);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  hid_t attrID = H5Acreate(grpID,"name",char8,aspaceID,H5P_DEFAULT,
                           H5P_DEFAULT);
#else
  hid_t attrID = H5Acreate(grpID,"name",char8,aspaceID,H5P_DEFAULT);
#endif
  std::valarray<char> buffer('\000',8*nspecies);
  for (int i=0; i<nspecies; ++i) name[i].copy(&buffer[i*8],8,0);
  H5Awrite(attrID,char8,&buffer[0]);
  H5Aclose(attrID);
  H5Sclose(aspaceID);
  H5Tclose(char8);
  // Write the nspecies attribute.
  aspaceID = H5Screate(H5S_SCALAR);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  attrID = H5Acreate(grpID,"nSpecies",H5T_NATIVE_INT,aspaceID,H5P_DEFAULT,
                           H5P_DEFAULT);
#else
  attrID = H5Acreate(grpID,"nSpecies",H5T_NATIVE_INT,aspaceID,H5P_DEFAULT);
#endif
  H5Awrite(attrID,H5T_NATIVE_INT,&nspecies);
  H5Aclose(attrID);
  H5Sclose(aspaceID);
  // Close the file
  H5Gclose(grpID);
  H5Fclose(fileID);
  H5close();
}
