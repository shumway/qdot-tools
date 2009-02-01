// $Id: Grid.h,v 1.5 2006/07/02 06:27:58 jshumwa Exp $
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
#ifndef __Grid_h_
#define __Grid_h_
#include <blitz/array.h>
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>
/** Class for grids.
 @author John Shumway
 @version $Revision: 1.5 $
 */
template<class T>
class Grid {
public:
  /// Typedefs.
  typedef blitz::TinyVector<float,3> Vec;
  typedef blitz::TinyVector<int,3> IVec;
  typedef blitz::Array<T,3> Array;
  /// Constructor.
  Grid(const IVec& extent, const Vec& delta)
  : data(extent), extent(extent), delta(delta), deltaInv(1./delta) {
    data=0.0;
  }
  /// Destructor.
  ~Grid() {}
  /// Get data.
  Array& getData() {return data;}
  /// Get data.
  const Array& getData() const {return data;}
  /// Add data value to grid point in lower left corner of cell containging pt.
  void addData(const Vec &pt, const T &value) {
    IVec index = floor(pt*deltaInv);
    (index += extent) %= extent; 
    data(index)+=value;
  }
  /// Get the grid extents.
  IVec getExtent() const {return extent;}
  /// Get grid value.
  T& operator()(const int i, const int j, const int k) {return data(i,j,k);}
  /// Get grid value.
  const T operator()(const int i, const int j, const int k) const
                    {return data(i,j,k);}
private:
  /// The data.
  Array data;
  /// The grid size. 
  const IVec extent;
  /// The grid spacing.
  const Vec delta, deltaInv;
};
#endif
