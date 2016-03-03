This is a collection of nanostructure simulation and modeling programs for self assembled quantum dots and other nanostructures.

## Atomic Nanostructure Layout Code: geom ##
This program, `geom`, is a simple utility for turning nanostructure descriptions into atomic layout. These atomic positions are common input for several types of quantum dot models, such as tight-binding or empirical pseudopotentials. The user specifies the size, shape, and composition of a coherent III-V or group IV semiconductor heterostructure. Shapes include lens, pyramidal, and conical dots, as well as quantum wells, and may be filled with different alloy or pure materials. Output is an [HDF5 file](http://www.hdfgroup.org/HDF5/), `struct.h5`, containing atom positions, chemical identities, lattice indexing, and neighbor tables. The object-oriented design makes it easy for users to add new shapes and materials to this utility.
  * Libraries used: HDF5, libxml2

## Strain Relaxation Code: relax ##
This program, `relax`, uses a conjugate-gradient algorithm to relax atomic positions in coherently strained semiconductor nanostructures. A typical simulation of millions of atoms in and around an embeded quantum dot runs in a few hours on a typical PC. Input is an HDF5 file, struct.h5, containing the atomic positions, identities, and neighbor tables. The program output, `relaxed.h5`, is the relaxed atomic positions in the same [HDF5](http://www.hdfgroup.org/HDF5/) format as the input file, suitable for input into an atomistic electronic structure model. Also calculates stress for each atom.
  * Libraries used: HDF5, libxml2, blitz

## Effective Mass Model Generator: getema ##
This program, `getema`, extracts an effective mass model `boffset.h5` from a relaxed nanostructure `relaxed,h5`. The strained band offsets in `boffset.h5` can be used as input to our _**[pi](http://code.google.com/p/pi-qmc)**_ path integral quantum Monte Carlo code.
  * Libraries used: HDF5, libxml2, blitz
