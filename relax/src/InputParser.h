// $Id: InputParser.h,v 1.10 2007/12/28 02:29:59 jshumwa Exp $
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

#include "XMLParser.h"
class AtomicPotential;
class Checkpoint;
class SpeciesH5;
class CoordsH5;
class SpeciesH5;
class NeighborsH5;
class SuperCellH5;
class TotalEnergy;
class CellRelaxer;

/** Input parser for relax code. Reads potential information and simulation
  paramters from relax.xml input file.
 
  Sample input file:
  \verbatim
  <AtomisticStructure>
    <Structure infile="struct.h5"/>
    <Relax maxiter="200" ftol="1e-6" outfile="relaxed.h5">
      <CellRelax interval="10" mode="zonly"/>
      <Checkpoint interval="10" file="checkpoint.h5"/>
    </Relax>
    <CalculateStresses outfile="relaxed.h5"/>
    <!-- To test consistency of forces and energies, use the following
    <TestForces/> -->
    <ForceModel>
      <VFF/>
    </ForceModel>
  </AtomisticStructure>
  \endverbatim
  @author John Shumway
  @version $Revision: 1.10 $
 */
class InputParser : public XMLParser {
public:
  /// Constructor.
  InputParser(const std::string& filename);
  /// Virtual destructor.
  ~InputParser();
  /// Parse potential, returns a new AtomicPotential.
  AtomicPotential* parsePotential();
  /// Parse method.
  void parse(const xmlXPathContextPtr&);
  /// Parse method.
  void parse() {parse(context);}
  /// Get the name of the input struct.h5 file.
  const std::string getStructInFileName() {return structInFileName;}
  /// Construct a Checkpoint object if tag is present.
  Checkpoint* parseCheckpoint(CoordsH5&, SpeciesH5&, NeighborsH5&,
                              SuperCellH5&);
  /// Construct a Checkpoint object if tag is present.
  CellRelaxer* parseCellRelaxer(CoordsH5&, NeighborsH5&, SpeciesH5&,
                                SuperCellH5&, AtomicPotential&, TotalEnergy&);
  /// Get the maximum number of ConjGradient iterations.
  int getMaxIter() {return maxiter;}
  /// Get the force tolerance for ConjGradient convergence.
  double getFTol() {return ftol;}
  /// Get the name of the file to write the relaxed structure to.
  std::string getRelaxedFilename() {return relaxedFilename;}
  /// Get the name of the file to write the stress and structure to.
  std::string getStressFilename() {return stressFilename;}
  /// Get the name of the file to write the surface stress and structure to.
  std::string getSurfaceFilename() {return surfaceFilename;}
  /// Get the interval between cell relaxations.
  int getCellRelaxInterval() {return cellRelaxInterval;}
  /// Get the maximum number of ConjGradient iterations.
  int getAnnealSweeps() {return annealSweeps;}
  /// Get the scale factor for line minimization.
  double getScale() {return scale;}
private:
  /// The xml document to parse from.
  xmlDocPtr doc;
  /// The context to parse from.
  xmlXPathContextPtr context;
  /// The name of the structure input file.
  std::string structInFileName;
  /// The maximum number of conjugant gradient iterations.
  int maxiter;
  /// The force tolerance for conjugate gradients convergence.
  double ftol;
  /// The filename to write the relaxed structure to.
  std::string relaxedFilename;
  /// The filename to write the stress and structure to.
  std::string stressFilename;
  /// The filename to write the surface stress and structure to.
  std::string surfaceFilename;
  /// The interval between cell relaxations.
  int cellRelaxInterval;
  /// The number of sweeps of Monte Carlo annealing.
  int annealSweeps;
  /// Scale for line minimization (set smaller than 1.0 if unstable).
  double scale;
};
