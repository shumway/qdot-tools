#include "SuperCellH5.h"
#include <hdf5.h>
#include <iostream>
#include <blitz/tinyvec-et.h>

SuperCellH5::SuperCellH5() : a1(0.0), a2(0.0), a3(0.0) {
}

SuperCellH5::SuperCellH5(const std::string& filename) {
  h5Read(filename);
}

void SuperCellH5::h5Read(const std::string& filename) {
  std::cout << "Reading supercell info from " << filename << std::endl;
  H5open();
  hid_t fileID = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  hid_t grpID = H5Gopen2(fileID,"superCell",H5P_DEFAULT);
#else
  hid_t grpID = H5Gopen(fileID,"superCell");
#endif
  hid_t attrID = H5Aopen_name(grpID,"a1");
  H5Aread(attrID,H5T_NATIVE_FLOAT,&a1);
  H5Aclose(attrID);
  attrID = H5Aopen_name(grpID,"a2");
  H5Aread(attrID,H5T_NATIVE_FLOAT,&a2);
  H5Aclose(attrID);
  attrID = H5Aopen_name(grpID,"a3");
  H5Aread(attrID,H5T_NATIVE_FLOAT,&a3);
  H5Aclose(attrID);
  H5Gclose(grpID);
  H5Fclose(fileID);
  H5close();
  computeRecipricalVectors();
}

void SuperCellH5::computeRecipricalVectors() {
  double invol=1/(a1[0]*(a2[1]*a3[2]-a2[2]*a3[1]) +
                  a1[1]*(a2[2]*a3[0]-a2[0]*a3[2]) +
                  a1[2]*(a2[0]*a3[1]-a2[1]*a3[0]));
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

void SuperCellH5::pbc(Vec& v) const {
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
  switch(mode) {
  case (NEW):
    fileID = H5Fcreate(filename.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
    break;
  case (APPEND):
    fileID = H5Fopen(filename.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
    break;
  }
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  hid_t grpID = H5Gcreate2(fileID,"superCell",0,H5P_DEFAULT,H5P_DEFAULT);
#else
  hid_t grpID = H5Gcreate(fileID,"superCell",0);
#endif
  const hsize_t dims[] = {3};
  hid_t aspaceID = H5Screate_simple(1,dims,0);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  hid_t attrID = H5Acreate2(grpID,"a1",H5T_NATIVE_FLOAT,aspaceID,H5P_DEFAULT,
                            H5P_DEFAULT);
#else
  hid_t attrID = H5Acreate(grpID,"a1",H5T_NATIVE_FLOAT,aspaceID,H5P_DEFAULT);
#endif
  H5Awrite(attrID,H5T_NATIVE_DOUBLE,&a1);
  H5Aclose(attrID);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  attrID = H5Acreate2(grpID,"a2",H5T_NATIVE_FLOAT,aspaceID,H5P_DEFAULT,
                      H5P_DEFAULT);
#else
  attrID = H5Acreate(grpID,"a2",H5T_NATIVE_FLOAT,aspaceID,H5P_DEFAULT);
#endif
  H5Awrite(attrID,H5T_NATIVE_DOUBLE,&a2);
  H5Aclose(attrID);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
  attrID = H5Acreate2(grpID,"a3",H5T_NATIVE_FLOAT,aspaceID,H5P_DEFAULT,
                      H5P_DEFAULT);
#else
  attrID = H5Acreate(grpID,"a3",H5T_NATIVE_FLOAT,aspaceID,H5P_DEFAULT);
#endif
  H5Awrite(attrID,H5T_NATIVE_DOUBLE,&a3);
  H5Aclose(attrID);
  H5Sclose(aspaceID);
  H5Gclose(grpID);
  H5Fclose(fileID);
  H5close();
}
