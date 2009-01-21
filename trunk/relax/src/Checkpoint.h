// $Id: Checkpoint.h,v 1.3 2004/04/27 21:02:47 jshumwa Exp $
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
#ifndef _Checkpoint_h_
#define _Checkpoint_h_
class CoordsH5;
class SpeciesH5;
class NeighborsH5;
class SuperCellH5;
#include <string>

/** Checkpoint a stress run. */
class Checkpoint {
public:
  /// Constructor.
  Checkpoint(const CoordsH5&, const SpeciesH5&, const NeighborsH5&, 
             const SuperCellH5&, const std::string& filename,
             const int interval);
  /// Write to file.
  void write(const int) const;
  /// Get the name of the checkpoint file.
  std::string getFileName() {return filename;}
private:
  const CoordsH5& coords;
  const SpeciesH5& species;
  const NeighborsH5& neighbors;
  const SuperCellH5& cell;
  const std::string filename;
  const int interval;
};
#endif
