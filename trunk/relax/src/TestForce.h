// $Id: TestForce.h,v 1.3 2007/03/14 19:47:50 jshumwa Exp $
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
#ifndef __TestForce_h_
#define __TestForce_h_

#include <blitz/tinyvec.h>
class Forces;
class TotalEnergy;
class CoordsH5;

/** Test the consistancy of forces and potentials.
 * @author John Shumway */
class TestForce {
public:
#ifdef ENABLE_FLOAT
  typedef blitz::TinyVector<float,3> Vec3;
#else
  typedef blitz::TinyVector<double,3> Vec3;
#endif
  /// Constructor.
  TestForce(Forces& f, TotalEnergy& en, CoordsH5& coords);
  /// Destructor.
  ~TestForce();
  /// Test the consistency of forces and potential.
  void test();
private:
  /// Computation and storage for the forces.
  Forces &f;
  /// Computation of total energy.
  TotalEnergy& en;
  /// The coordinates to be optimized.
  CoordsH5& coords;
};
#endif
