// $Id:
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
#include "CommentH5.h"
#include "hdf5.h"
#include <iostream>
#include <fstream>

CommentH5::CommentH5() {
  comment="";
}

CommentH5::CommentH5(const std::string& filename){
  std::string buffer;
  std::ifstream in(filename.c_str());
  while(!std::getline(in, buffer).eof())
    comment+=buffer+="\n";
//std::cout << comment << std::endl;
}

void CommentH5::h5Write(const std::string& filename, const int mode) const {
  std::cout << "Writing comment to " << filename << std::endl;
  H5open();
  hid_t fileID;
  switch(mode) {
  case (NEW):
    fileID = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    break;
  case (APPEND):
    fileID = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    break;
  }
  hid_t dspaceID = H5Screate(H5S_SCALAR);
  hid_t dtypeID = H5Tcopy(H5T_C_S1);
  size_t size = comment.size();
  H5Tset_size(dtypeID, size);
  hid_t dsetID = H5Dcreate(fileID, "comment", dtypeID, dspaceID, H5P_DEFAULT);
  H5Dwrite(dsetID, dtypeID, H5S_ALL, H5S_ALL, H5P_DEFAULT, comment.data());
  H5Dclose(dsetID);
  H5Tclose(dtypeID);
  H5Sclose(dspaceID);
  H5Fclose(fileID);
  H5close();
}
