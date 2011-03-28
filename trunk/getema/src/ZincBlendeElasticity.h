// $Id: ZincBlendeElasticity.h,v 1.2 2006/07/02 06:27:58 jshumwa Exp $
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
#ifndef __ZincBlendeElasticity_h_
#define __ZincBlendeElasticity_h_

class GridFactory;
class CompositionGrids;
#include "Grid.h"
#include "ElasticConstants.h"

/// Class for storing and calculating elasticity for Cubic Alloys.
/// @author John Shumway
/// @version $Revision: 1.2 $
/// @bug Hard coded for common anion, common cation, or binary alloys.
class ZincBlendeElasticity : public ElasticConstants {
public:
  /// Storage for elastic coefficients.
  struct coef {
    coef(double c11, double c22, double c44) : c11(c11), c12(c12), c44(c44) {}
    double c11, c12, c44;
  };
  /// Constructor
  ZincBlendeElasticity(const GridFactory&, const coef& mat1, const coef& mat2,
    const int ind1, const int ind2);
  /// Virtual destructor
  virtual ~ZincBlendeElasticity();
  /// Calculate the grids.
  void calculate(const CompositionGrids&);
  /// Get the matrix of elastic constants for a grid point.
  virtual const Mat& operator() (const int i, const int j, const int k) const;
  /// Elastic constants for some common materials.
  static const coef GAAS, INAS, GASB, SI, GE;
private:
  /// Elastic constants  for the two materials.
  const coef mat1, mat2;
  /// Factory for creating new Grid objects.
  const GridFactory &factory;
  /// Indicies for the components of alloy materials.
  const int ind1, ind2;
  /// Cubic elasticity constants on a grid.
  Grid<float> *c11,*c12,*c44;
  /// Elasticity matrix.
  mutable Mat mat;
  /// Conversion for tenMbar to atomic units.
  static const double tenMbar;
};
#endif
