// $Id: SuperCell.cc,v 1.3 2006/08/08 18:02:48 jshumwa Exp $
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
#include "SuperCell.h"

SuperCell::SuperCell(const int nx, const int ny, const int nz,
                     const double a, const double strainxx, 
                     const double strainyy, const double strainzz, 
                     const Material* material) :
  nx(nx), ny(ny), nz(nz), a(a), strainxx(strainxx), strainyy(strainyy),
  strainzz(strainzz), bulkMaterial(material) {
}
