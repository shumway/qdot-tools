// $Id: StillWeb.h,v 1.3 2004/10/15 02:43:55 jshumwa Exp $
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
#ifndef __StillWeb_h_
#define __StillWeb_h_
#include <valarray>
#include "AtomicPotential.h"

/** Stillinger-Weber potential.
 *
 * Taken from Stillinger and Weber, PRB 31, 5262-5271 (1985),
 * online at http://link.aps.org/abstract/PRB/v31/p5262.
 * The pairwise interaction is a function of separation,
 * @f[
 * f_2(r)= A(Br^{-p}-r^{-q})\exp[(r-a)^{-1}], r<a
 * @f]
 * and zero for larger separations.
 *
 * There is interaction between angular bonds,
 * @f[
 * h(r_{ij},r_{ik},\theta_{jik}) = \lambda\exp[\gamma(r_{ij}-a)^{-1} +
 * \gamma(r_{ik}-a)^{-1}](\cos\theta_{jik}+\frac{1}{3})^2.
 * @f] 
 *
 * Si-Si lattice paramters is 10.246 Bohr radii (T=0K). 
 *
 * @todo Make a table of sqrt values for faster evaluation of h.
 * @author John Shumway */
class StillWeb : public AtomicPotential {
public:
  /// Constructor for Si/Ge.
  StillWeb();
  /// Destructor.
  virtual ~StillWeb() {};
  /// Evaluate the pair binding energy.
  virtual double bind(const int i, const int j) const;
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
    return a/param[ind(i,j)].sigmainv;
  }

  /// Subclass for storage of species dependent parameters.
  class SWParam {
  public:
    double lambda;
    double sigmainv;
    double epsilon;
  };

  double A;
  double B;
  double p;
  double q;
  double a;
  double gamma;
  
private:
  std::valarray<SWParam> param;
};
#endif
