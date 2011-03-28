// $Id: VFF.h,v 1.2 2004/10/15 00:46:26 jshumwa Exp $
/*
    Copyright (C) 2004 Matt Harowitz and John Shumway

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
#ifndef __VFF_h_
#define __VFF_h_
#include <valarray>
#include "AtomicPotential.h"

/** Generalized Valence Force Field potential.
 *
 * The pairwise interaction is a function of separation,
 * @f[
 * f_2\left(r_{ij}\right)=\frac{3}{8}\frac{\alpha_{ij}
 * \left({r_{ij}}^2-{d_{ij}}^2\right)}{d_{ij}}.
 * @f]  *
 * There is interaction between angular bonds,
 * @f[
 * h\left(r_{ij},r_{ik}\right)=\frac{3}{8}\frac{\beta_{ijk}
 * \left(r_{ij}r_{ik}\cos\left(\theta_{ijk}\right)
 * -d_{ij}d_{ik}\cos\left(\phi_{ijk}\right)\right)}
 * {d_{ij}d_{ik}}.
 * @f]
 *
 * @todo add higher order alpha, sigma.
 * @todo get alpha, beta, etc from elasticity coefficients
 * @author Matthew Harowitz, John Shumway */

 class VFF : public AtomicPotential {
 public:
  /// Constructor for InGaAs.
  VFF();
  /// Destructor.
  virtual ~VFF() {};
  /// Evaluate the pair binding energy.
  virtual double bind(const int i, const int j) const {return 0;}
  /// Evaluate the pair potential.
  virtual double f2(const double r, const int i, const int j) const;
  /// Evaluate the derivative of the pair potential.
  virtual double drF2(const double r, const int i, const int j) const;
  /// Evaluate the angle potential.
  virtual double h(const double r1, const double r2, const double costheta,
           const int i, const int j, const int k) const;
  /// Evaluate the derivatives of the angle potential.
  virtual void dh(const double r1, const double r2, const double costheta,
          const int i, const int j, const int k,
	  double& dr1H, double& dr2H, double& dthetaH) const;
  /// Return the cutoff value.  Not needed for VFF.  Included to match StillWeb
  virtual double getCutoff(const int i, const int j) const {
    return 10.;
  }

  /// Subclass for storage of species dependent parameters.
class VFFParam {
  public:
    double d;
    double alpha1;
 };
    double cosphi[6];
    double beta[6];

private:
  std::valarray<VFFParam> param;
};
#endif
