// $Id: Structure.cc,v 1.2 2004/04/22 20:09:44 jshumwa Exp $
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
#include "Structure.h"
#include "LatticeIndexingH5.h"
#include "Material.h"
#include "SpeciesH5.h"
#include "Vec3.h"

Structure::Structure(const Material* material) : material(material) {
}

void Structure::build(SpeciesH5& species,
                      const LatticeIndexingH5& indexing) const {
  int natoms = species.species.size();
  for (int iatom=0; iatom<natoms; ++iatom) {
    Vec3 pt(indexing.getIdealPositionOverA(iatom));
    if ( isPointInStruct(pt) ) {
      species.species[iatom] =
          ((indexing.getIPos(iatom)&1) == 0)
             ? material->getAnionIndex(pt)
             : material->getCationIndex(pt);
    }
  }
}
