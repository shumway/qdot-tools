// $Id: AtomicPotential.h,v 1.2 2004/10/15 00:46:26 jshumwa Exp $
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
#ifndef __AtomicPotential_h_
#define __AtomicPotential_h_
#include <valarray>

/** Empirical inter-atomic potential.
 *
 * The pairwise interaction is a function of separation,
 * @f$ f_2(r) @f$.
 * There is interaction between angular bonds,
 * @f$ * h(r_{ij},r_{ik},\theta_{jik}) * @f$ .
 *
 * @author John Shumway */
class AtomicPotential {
public:
  /// Constructor.
  AtomicPotential(){};
  /// Destructor.
  virtual ~AtomicPotential(){};
  /// Evaluate the pair binding energy.
  virtual double bind(const int i, const int j) const=0;
  /// Evaluate the pair potential.
  virtual double f2(const double r, const int i, const int j) const=0;
  /// Evaluate the derivative of the pair potential.
  virtual double drF2(const double r, const int i, const int j) const=0;
  /// Evaluate the angle potential.
  virtual double h(const double r1, const double r2, const double costheta,
                   const int i, const int j, const int k) const=0;
  /// Evaluate the derivatives of the angle potential.
  virtual void dh(const double r1, const double r2, const double costheta,
          const int i, const int j, const int k,
	  double& dr1H, double& dr2H, double& dthetaH) const=0;
  /// Return the cutoff value.
  virtual double getCutoff(const int i, const int j) const=0;
protected:
  int nspec;
  static int ind(const int i, const int j) {
    return (i<j)?i+j*(j+1)/2:j+i*(i+1)/2;
  }
  
};
#endif
