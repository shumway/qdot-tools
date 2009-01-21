// $Id: ConjGrad.h,v 1.7 2007/12/28 02:19:40 jshumwa Exp $
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
#ifndef __ConjGrad_h_
#define __ConjGrad_h_

#include "Minimize1D.h"
#include <blitz/tinyvec.h>
#include <valarray>
class Forces;
class TotalEnergy;
class CoordsH5;
class Checkpoint;
class CellRelaxer;

/** Conjugate gradients algorthm for minimizing energy. 
 * @author John Shumway */
class ConjGrad : public Minimize1D::Function {
public:
#ifdef ENABLE_FLOAT
  typedef blitz::TinyVector<float,3> Vec3;
#else
  typedef blitz::TinyVector<double,3> Vec3;
#endif
  /// Constructor.
  ConjGrad(Forces& f, TotalEnergy& en, CoordsH5& coords,
    const Checkpoint* checkpt, const double ftol=2e-6, const int itermax=1000,
    CellRelaxer* cellRelaxer=0, const int cellRelaxSteps=1,
    const double scale=1.0);
  /// Destructor.
  ~ConjGrad();
  /// Minimize the energy using conjugate gradients algorithm.
  void minimize();
  /// Energy function for 1D minimizer.
  double operator() (double const);
private:
  /// One dimensional funciton minimizer.
  Minimize1D minimizer;
  /// Minimize the energy for coords along direction.
  void lineMinimize(); 
  /// Computation and storage for the forces.
  Forces& f;
  /// Computation of total energy.
  TotalEnergy& en;
  /// Direction for line minimization.
  std::valarray<Vec3> direction;
  /// The coordinates to be optimized.
  CoordsH5& coords;
  /// Last x value for operator() (used for maintaining accurate coordinates).
  double xLast;
  /// Checkpointer.
  const Checkpoint* checkpt;
  /// Maximum force tolerance for convergence.
  const double ftol;
  /// Maximum number of conjugant gradients iterations.
  const int itermax;
  /// Cell relaxation algorithm.
  CellRelaxer* cellRelaxer;
  /// Steps to go before relaxing the cell.
  const int cellRelaxSteps;
  /// Scale factor for line minimization (make smaller for better stability).
  const double scale;
};
#endif
