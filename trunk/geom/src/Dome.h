// $Id: Dome.h,v 1.4 2004/06/14 23:39:57 jshumwa Exp $
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
#ifndef __Dome_h_
#define __Dome_h_
#include "Structure.h"
#include <blitz/array.h>
#include <vector>

class Material;

/** Representation for a dome shape quantum dot structure. 

For a detailed experimental descripton of Ge/Si and GeSi/Si dome shapes
under different growth conditions,
see Chaparro, Zhang, Drucker, and Smith, APL <b>87</b>, 2245 (2000),
http://dx.doi.org/10.1063/1.372168.

Facets are defined by three indicies [ijk] (i>=j) and a height at which
the etended facet would cross the center of the dot.  We take the
growth direction to be [001], so you might want to reorder indicies
in the xml input as [kij] if you are comparing to experimental 
descriptions of growth in the [100] direction:
@verbatim
<!-- 30 nm by 5.4 nm dot -->
<Dome x="64" y="64" z="45.2" material="Ge">
  <Facet k="5" i="1" j="1" height="10.0"/>
  <Facet k="3" i="1" j="1" height="13.0"/>
  <Facet k="2" i="1" j="0" height="13.8"/>
</Dome>
@endverbatim
The height of an [ijk] facet with center height h is
@f[ z = h - \frac{i}{k} x - \frac{j}{k}y. @f]
@version $Revision: 1.4 $
@author John Shumway */
class Dome : public Structure {
public:
  /// Internal class for representing a facet on a dome.
  class Facet {
    public:
    /// Construct Facet by providing indicies and center height.
    Facet(int i, int j, int k, double h) : i(i), j(j), k(k), height(h) {} 
    /// Facet indexing.
    int i,j,k;
    /// (Extrapolated) height of facet at center of dome.
    double height;
  };
  /// Typedefs.
  typedef std::vector<Facet> FArray;
  typedef blitz::TinyVector<double,3> Plane;
  typedef blitz::Array<Plane,1> PArray;
  /** Constructor.  Provide the x,y,z coordinates for the center
      of the base, the height, facets, and material. Assumes four-fold 
      symmetry for facets.  */
  Dome(const double x, const double y, const double z, 
       const FArray&, const Material* material);
  /** Is a given point inside the pyramid? */
  bool isPointInStruct(const Vec3& pt) const;
protected:
  /** The x-coordinate of the center of the base, measured in eight 
      atom lattice units. */
  const double x;
  /** The y-coordinate of the center of the base, measured in eight 
      atom lattice units. */
  const double y;
  /** The z-coordinate of the center of the base, measured in eight 
      atom lattice units. */
  const double z;
  /// The height of the dome.
  double height; 
  /// The number of facets (divided by four-fold symmetry).
  const int nfacet; 
  /// The planes bounding the dome.
  PArray boundingPlane;
};
#endif
