// $Id: StructH5.h,v 1.3 2004/04/28 01:54:59 jshumwa Exp $
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
#ifndef __StructH5_h_
#define __StructH5_h_

/** Base class of representations of struct.h5 elements.
    See subclasses for details of the Struct.h5 file format.

    For information about the HDF5 file format and libraries, see
    \URL{http://hdf.ncsa.uiuc.edu/HDF5}.
    @author John Shumway */

#include <string>

class StructH5 {
public:
  /** Write the object to a struct.h5 file. */
  virtual void h5Write(const std::string& filename, 
                       const int mode=APPEND) const=0;
  static const int NEW=0;
  static const int APPEND=1;
  static const int REWRITE=2;
};

#endif
