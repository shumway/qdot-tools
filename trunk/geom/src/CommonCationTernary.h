// $Id: CommonCationTernary.h,v 1.3 2004/05/10 21:56:26 jshumwa Exp $
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
#ifndef __CommonCationTernary_h_
#define __CommonCationTernary_h_

#include "Material.h"
class RandomNR;

/** Material class for a common cation ternary alloy.  Construct the material, 
    then set the integer IDs before calling getCationIndex or getCationIndex.  
    Composition is randomly generated, and repeated calls to the same anion
    site will generate random results. 
    @author John Shumway */
class CommonCationTernary : public Material {
public:
  /** Constructor.  The  argument x1 = (1-x2) is the fraction of 
      anion1 atoms.  The last argument seed is an integer seed for the
      random number generator rand()*/
  CommonCationTernary(const std::string& name, const std::string& anion1Name, 
                  const std::string& anion2Name, const std::string& cationName,
                  const double x1, const int seed, RandomNR& rand);
  /** Return the index for an anion at point p. */
  virtual int getAnionIndex(const Vec3& p) const;
  /** Return the index for a cation at point p.  Not reproducable,
      and repeated calls will generate different results. */
  virtual int getCationIndex(const Vec3& p) const;
protected:
  /** The percentage of anion1 atoms. */
  double x1;
  /** The random number generator. */
  RandomNR& rand;
};
#endif
