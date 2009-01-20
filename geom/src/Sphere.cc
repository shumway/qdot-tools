// $Id: Sphere.cc,v 1.2 2004/04/22 20:09:44 jshumwa Exp $
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
#include "Sphere.h"
#include "Vec3.h"

Sphere::Sphere(const double x, const double y, const double z, 
               const double radius, const Material* material) :
  Structure(material), x(x), y(y), z(z), radius(radius), r2(radius*radius) {
}

bool Sphere::isPointInStruct(const Vec3& pt) const {
  return (pt.x-x)*(pt.x-x)+(pt.y-y)*(pt.y-y)+(pt.z-z)*(pt.z-z) < r2;
}
/*
void Sphere::build(SpeciesH5& species, LatticeIndexingH5& indexing) const {
  int natoms = species.species.size();
  double r2=radius*radius;
  for (int iatom=0; iatom<natoms; ++iatom) {
    LatticeIndexingH5::indices index = indexing.index[iatom];
    double atomx = index.i + indexing.getOffset(index.ipos)[0];
    double atomy = index.j + indexing.getOffset(index.ipos)[1];
    double atomz = index.k + indexing.getOffset(index.ipos)[2];
    if ( (atomx-x)*(atomx-x)+(atomy-y)*(atomy-y)+(atomz-z)*(atomz-z) < r2 ) {
      if ((index.ipos&2) == 0)  {
        species.species[iatom]=material->getAnionIndex(atomx,atomy,atomz);
      } else { 
        species.species[iatom]=material->getCationIndex(atomx,atomy,atomz);
      }
    }
  }
}
*/
