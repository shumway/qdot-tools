// $Id: Binary.cc,v 1.2 2004/04/22 20:09:44 jshumwa Exp $
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
#include "Binary.h"

Binary::Binary(const std::string& name, const std::string& anionName,
                                        const std::string& cationName)
  : Material(name,2) {
  speciesName[0]=anionName;
  speciesName[1]=cationName;
  speciesID[0]=1;
  speciesID[1]=2;
}

int Binary::getAnionIndex(const Vec3& p) const {
  return speciesID[0];
}

int Binary::getCationIndex(const Vec3& p) const {
  return speciesID[1];
}
