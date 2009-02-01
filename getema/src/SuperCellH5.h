#ifndef __SuperCellH5_h_
#define __SuperCellH5_h_

#include <string>
#include "StructH5.h"
#include <blitz/tinyvec.h>

/** The supercell.

    Layout in the H5 file:
    <ul> <li> superCell <ul>
      <li> a1 (attribute, real[3]) - 1st supercell lattice vector
      <li> a2 (attribute, real[3]) - 2nd supercell lattice vector
      <li> a3 (attribute, real[3]) - 3rd supercell lattice vector
    </ul></ul>
    @author John Shumway */
class SuperCellH5 : public StructH5 {
public:
  /// Typedefs.
  typedef blitz::TinyVector<float,3> Vec;
  /// Constructor 
  SuperCellH5();
  /// Construct from file.
  SuperCellH5(const std::string& filename);
  /// Supercell lattice vector @f$a_1@f$.
  Vec a1;
  /// Supercell lattice vector @f$a_2@f$.
  Vec a2;
  /// Supercell lattice vector @f$a_3@f$.
  Vec a3;
  /// Read the supercell from a struct.h5 file.
  void h5Read(const std::string& filename);
  /// Write the supercell to a struct.h5 file.
  void h5Write(const std::string& filename, const int mode=APPEND) const;
  /// Compute reciprical lattice vectors.
  void computeRecipricalVectors();
  /** Set to the smallest displacement with PBC.
   * @bug Only projects to primative unit cell, not-neccesarily
   * the smallest vector. */
  void pbc(Vec&) const;
  /// Reciprical supercell lattice vector @f$b_1,b_2,b_3@f$.
  Vec b1,b2,b3;
};
#endif
