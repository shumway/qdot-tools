// $Id: BinaryAlloy.cc,v 1.2 2004/04/22 20:09:44 jshumwa Exp $
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
#include "BinaryAlloy.h"
#include <cmath>

BinaryAlloy::BinaryAlloy(const std::string& name, const std::string& atom1,
                         const std::string& atom2, const double x1)
  : Material(name,2), x1(x1) {
  speciesName[0]=atom1;
  speciesName[1]=atom2;
  speciesID[0]=1;
  speciesID[1]=2;
}

int BinaryAlloy::getAnionIndex(const Vec3& p) const {
  return (rand()<x1*RAND_MAX) ? speciesID[0] : speciesID[1];
}

int BinaryAlloy::getCationIndex(const Vec3& p) const {
  return (rand()<x1*RAND_MAX) ? speciesID[0] : speciesID[1];
}
