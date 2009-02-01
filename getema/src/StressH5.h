// $Id: StressH5.h,v 1.1 2004/04/28 01:54:59 jshumwa Exp $
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
#ifndef __StressH5_h_
#define __StressH5_h_
#include "StructH5.h"
class Stress;

/** Class for HDF5 I/O of stress tensor.
 * @author John Shumway */
class StressH5 : public StructH5 {
public:
  /// Constructor.
  StressH5(Stress& stress) : stress(stress) {};
  /// Write the obejct to a struct.h5 file.
  virtual void h5Read(const std::string& filename);
  /// Write the obejct to a struct.h5 file.
  virtual void h5Write(const std::string& filename, const int mode=APPEND)const;
private:
  /// Reference to stress object.
  Stress& stress;
};
#endif
