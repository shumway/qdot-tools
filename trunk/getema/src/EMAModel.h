// $Id: EMAModel.h,v 1.2 2006/07/02 06:27:58 jshumwa Exp $
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
#ifndef __EMAModel_h_
#define __EMAModel_h_

/** Base class for calculating and storing position-dependent effective mass
model parameters.
@author John Shumway
@version $Revision: 1.2 $ */
class EMAModel {
public:
  /// Get the electron offset grid.
  virtual const Grid<float>& getVeGrid() const=0;
  /// Get the hole offset grid.
  virtual const Grid<float>& getVhGrid() const=0;
  /// Calculate the grids.
  virtual void calculate(const CompositionGrids&, const StrainGrid&)=0;
};
#endif
