// $Id: Spring.h,v 1.3 2006/05/29 15:48:45 jshumwa Exp $
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
#ifndef __Spring_h_
#define __Spring_h_
#include <valarray>
#include "AtomicPotential.h"

/** Spring potential.
 *
 * The pairwise interaction is a function of separation,
 * @f[
 * f_2(r)= \frac{1}{2}k(r-r_0)^2
 * @f]
 * and zero for larger separations.
 *
 * There is no interaction between angular bonds,
 * @f[
 * h(r_{ij},r_{ik},\theta_{jik}) = 0
 * @f] 
 *
 * @author John Shumway */
class Spring : public AtomicPotential {
public:
  /// Constructor for Si/Ge.
  Spring();
  /// Destructor.
  virtual ~Spring() {};
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
  /// Return the cutoff value.
  virtual double getCutoff(const int i, const int j) const {
    return rcut;
  }

private:
  /// Cutoff (rcut=5).
  static const double rcut;
  /// Spring constant (k=1).
  static const double k;
  /// Bond distance (r0=4.4366).
  static const double r0;
};
#endif
