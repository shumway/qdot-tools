// $Id: Minimize1D.h,v 1.2 2006/05/29 15:48:45 jshumwa Exp $
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
#ifndef __Minimize1D_h_
#define __Minimize1D_h_

/**
 *
 * @author John Shumway */
class Minimize1D {
public:
  /// Utility base class for 1D function to minimize.
  class Function {
  public:
    /// Function to minimize.
    virtual double operator()(const double)=0;
  };
  /// Constructor.
  Minimize1D(Function& f);
  /// Destructor.
  ~Minimize1D();
  /// Bracket a minimum.  Input is xa,xb,fa,fb; output is xa,xb,xc,fa,fb,fc.
  void bracketMin(double& ax, double& bx, double& cx,
                  double& fa, double& fb, double& fc);
  /// Perform a golden section search.
  double goldenSearch(const double ax, const double bx, const double cx, 
                      double& xmin);
  /// Perform a Brent's method search.
  double brentSearch(const double ax, const double bx, const double cx, 
                      double& xmin);
private:
  /// Function to minimize.
  Function& f;
  /// Golden ratio, 0.618034.
  static const double GOLD;
};
#endif
