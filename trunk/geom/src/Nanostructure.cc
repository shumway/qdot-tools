// $Id: Nanostructure.cc,v 1.3 2004/06/11 00:48:52 jshumwa Exp $
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
#include "Nanostructure.h"
#include "Structure.h"
#include "Material.h"

Nanostructure::~Nanostructure() {
  for (std::vector<const Structure*>::iterator s=structs.begin();
       s!=structs.end(); ++s) {
    delete *s;
  }
}

void Nanostructure::addStructure(const Structure* s) {
  structs.push_back(s);
}

int Nanostructure::getNumber() const {
  return structs.size();
}

const Structure* Nanostructure::getStructure(const int i) const {
  return structs[i];
}

void Nanostructure::addMaterial(const Material* mat) {
  material.push_back(mat);
  // Read the species names for the material, and set the species IDs.
  for (int ispecies=0; ispecies < mat->speciesName.size(); ++ispecies) {
    std::string name(mat->speciesName[ispecies]);
    int id=0;
    for (id=0; id<speciesName.size() && speciesName[id]!=name; ++id);
    if (id==speciesName.size()) speciesName.push_back(name);
    mat->speciesID[ispecies]=id+1;
  }
}

const Material* Nanostructure::getMaterial(const std::string& name) const {
  for (int imat=0; imat<material.size(); ++imat) {
    if (material[imat]->getName()==name) return material[imat];
  }
  return 0;
}

int Nanostructure::getSpeciesIndex(const std::string& name) const {
  for (unsigned int ispec=0; ispec<speciesName.size(); ++ispec) {
    if (speciesName[ispec]==name) return ispec+1;
  }
  return 0;
}
