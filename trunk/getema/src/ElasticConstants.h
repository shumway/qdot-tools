// $Id: ElasticConstants.h,v 1.1.1.1 2004/05/03 16:49:21 jshumwa Exp $
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
#ifndef __ElasticConstants_h_
#define __ElasticConstants_h_
class CompositionGrids;

/** Base class for calculating and storing position-dependent elastic constants.
@author John Shumway
@version $Revision: 1.1.1.1 $ */
class ElasticConstants {
public:
  /// Typedef.
  typedef blitz::TinyMatrix<float,6,6> Mat;
  /// Calculate the elastic contants for a given composition.
  virtual void calculate(const CompositionGrids&)=0;
  /// Get the matrix of elastic constants for a grid point.
  virtual const Mat& operator() (const int i, const int j, const int k) const=0;
};
#endif
