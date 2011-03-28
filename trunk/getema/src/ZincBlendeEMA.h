// $Id: ZincBlendeEMA.h,v 1.4 2006/07/03 00:05:58 jshumwa Exp $
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
#ifndef __ZincBlendeEMA_h_
#define __ZincBlendeEMA_h_

class GridFactory;
class CompositionGrids;
class StrainGrid;
#include "Grid.h"
#include "EMAModel.h"
#include <blitz/tinyvec.h>
#include <complex>

/// Class for storing and calculating elasticity for ZincBlend Alloys.
/// Uses parameters from Chris G. Van de Walle, 
/// "Band lineups and deformation potentials in the model-solid theory," 
/// <a href="http://link.aps.org/abstract/PRB/v39/p1871">
/// PRB 39, 1871-1883 (1989)</a>.
/// @author John Shumway
/// @version $Revision: 1.4 $
/// @bug Hard coded for common anion, common cation, or binary alloys.
class ZincBlendeEMA : public EMAModel {
public:
  /// Typedef.
  typedef blitz::TinyVector<int,3> IVec;
  typedef blitz::Array<std::complex<double>,1> CArray1;
  typedef blitz::Array<std::complex<double>,2> CArray2;
  typedef blitz::Array<double,1> Array1;
  /// Structure for band parameters.
  struct param {
    param(double ec, double ev, double ac, double av, double b, double d,
          double so, double xi) : ec(ec*eVtoHa), ev(ev*eVtoHa), ac(ac*eVtoHa),
      av(av*eVtoHa), b(b*eVtoHa), d(d*eVtoHa), so(so*eVtoHa), xi(xi*eVtoHa) {}
    double ec,ev,ac,av,b,d,so,xi;
  };
  /// Constructor
  ZincBlendeEMA(const GridFactory&, const param &mat1, const param &mat2,
                const int ind1, const int ind2);
  /// Virtual destructor
  virtual ~ZincBlendeEMA();
  /// Calculate the grids.
  virtual void calculate(const CompositionGrids&, const StrainGrid& strain);
  /// Get the electron offset grid.
  virtual const Grid<float>& getVeGrid() const {return ve;}
  /// Get the hole offset grid.
  virtual const Grid<float>& getVhGrid() const {return vh;}
  /// Parameters for some materials.
  static const param GAAS, INAS, GASB, SI, GE;
private:
  /// Factory for creating new Grid objects.
  const GridFactory &factory;
  /// Indicies for the three components of common anion ZB materials.
  const int ind1, ind2;
  /// Parameters for the two materials.
  const param mat1, mat2;
  /// Grids for electron and hole band energies.
  Grid<float> &ve,&vh;
  /// Conversion eVtoHa.
  static const double eVtoHa;
};
#endif
