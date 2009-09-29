// $Id: SuperCellH5.cc,v 1.3 2004/05/18 22:33:59 matth2 Exp $
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
#include "SuperCellH5.h"
#include "hdf5.h"
#include <iostream>

SuperCellH5::SuperCellH5() {
  a1[0]=0; a1[1]=0; a1[2]=0;
  a2[0]=0; a2[1]=0; a2[2]=0;
  a3[0]=0; a3[1]=0; a3[2]=0;
}

void SuperCellH5::h5Write(const std::string& filename, const int mode) const {
  std::cout << "Writing supercell info to " << filename << std::endl;
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
  hid_t grpID = H5Gcreate2(fileID, "superCell", H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);
#else
  hid_t grpID = H5Gcreate(fileID, "superCell", 0);
#endif
  const hsize_t dims[] = {3};
  hid_t aspaceID = H5Screate_simple(1,dims,0);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  hid_t attrID = H5Acreate(grpID,"a1",H5T_NATIVE_FLOAT,aspaceID,H5P_DEFAULT,
                           H5P_DEFAULT);
#else
  hid_t attrID = H5Acreate(grpID,"a1",H5T_NATIVE_FLOAT,aspaceID,H5P_DEFAULT);
#endif
  H5Awrite(attrID,H5T_NATIVE_DOUBLE,a1);
  H5Aclose(attrID);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  attrID = H5Acreate(grpID,"a2",H5T_NATIVE_FLOAT,aspaceID,H5P_DEFAULT,
                           H5P_DEFAULT);
#else
  attrID = H5Acreate(grpID,"a2",H5T_NATIVE_FLOAT,aspaceID,H5P_DEFAULT);
#endif
  H5Awrite(attrID,H5T_NATIVE_DOUBLE,a2);
  H5Aclose(attrID);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  attrID = H5Acreate(grpID,"a3",H5T_NATIVE_FLOAT,aspaceID,H5P_DEFAULT,
                           H5P_DEFAULT);
#else
  attrID = H5Acreate(grpID,"a3",H5T_NATIVE_FLOAT,aspaceID,H5P_DEFAULT);
#endif
  H5Awrite(attrID,H5T_NATIVE_DOUBLE,a3);
  H5Aclose(attrID);
  H5Sclose(aspaceID);
  H5Gclose(grpID);
  H5Fclose(fileID);
  H5close();
}
