// $Id: Nanostructure.h,v 1.3 2004/06/11 00:48:52 jshumwa Exp $
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
#ifndef __Nanostructure_h_
#define __Nanostructure_h_

#include <vector>
#include <string>
class Structure;
class Material;

/** Representation of a nanostrucure by a collection of Structure objects.
    @author John Shumway */
class Nanostructure {
public:
  /** Destrucutor.  Deletes all strucutal objects. */
  ~Nanostructure();

  /** Add a strutural object to the nanostructure. */
  void addStructure(const Structure*);

  /** Add a material. */
  void addMaterial(const Material*);

  /** Get a reference to a material by name. */
  const Material* getMaterial(const std::string& name) const;

  /** Get the number of structural objects in the nanostructure. */ 
  int getNumber() const;

  /** Get a pointer a structure using an index. */
  const Structure* getStructure(const int index) const;

  /** Get the number of species. */
  int getNSpecies() const {return speciesName.size();}

  /** Get species name. */
  const std::string& getSpeciesName(int i) const {return speciesName[i];}

  /// Get species index from a name.
  int getSpeciesIndex(const std::string&) const;

protected:
  /** The list of structures. */
  std::vector<const Structure*> structs;
  /** The list of materials. */
  std::vector<const Material*> material;
  /** The list of species names. */
  std::vector<std::string> speciesName;
};
#endif
