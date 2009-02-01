#ifndef __SpeciesH5_h_
#define __SpeciesH5_h_

#include <vector>
#include <string>
#include <blitz/array.h>
#include "StructH5.h"

/** The species (atom types).

    Layout in the H5 file:
    \begin{itemize}
    \item species \begin{itemize}
      \item name (attribute, string[nSpecies]) - Species names.
      \item nSpecies (attribute, integer) - Number of species.
      \item species (attribute, integer[nPart]) - Species type of each atom
              (enumerated as integers starting from 1).
    \end{itemize}\end{itemize}
    @author John Shumway
    @version $Revision: 1.1.1.1 $ */
class SpeciesH5 : public StructH5 {
public:
  /// Typedefs.
  typedef blitz::Array<int,1> IArray;
  /// Constructor.
  SpeciesH5(const int nspecies, const int natoms);
  /// Construct from a file.
  SpeciesH5(const std::string& filename);
  /// Read the supercell from a struct.h5 file.
  void h5Read(const std::string& filename);
  /// Write the supercell to a struct.h5 file.
  void h5Write(const std::string& filename, const int mode=APPEND) const;
  /// Get the number of species.
  int getNSpecies() const {return name.size();}
  /// Get species name by index.
  std::string& getSpeciesName(const int i) {return name[i];}
  /// Get the number of particles.
  int getNPart() const {return species.size();}
  /// Get species type by particle index.
  int operator()(const int i) const {return species(i);}
private:
  /// The species names.
  std::vector<std::string> name;
  /// The species indexing.
  mutable IArray species;
};
#endif
