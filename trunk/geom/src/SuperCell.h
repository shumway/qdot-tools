// $Id: SuperCell.h,v 1.3 2006/08/08 18:02:48 jshumwa Exp $
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
#ifndef __SuperCell_h_
#define __SuperCell_h_

class Material;

/** Supercell for the nanostructure.
    @author John Shumway */
class SuperCell {
public:
  /** Constructor. */
  SuperCell(const int nx, const int ny, const int nz, const double a,
            const double strainxx, const double strainyy, const double strainzz,
            const Material* material);
  /** Get X-Extent of the supercell in eight-atom primatives. */
  int getNX() const {return nx;}
  /** Get Y-Extent of the supercell in eight-atom primatives. */
  int getNY() const {return ny;}
  /** Get Z-Extent of the supercell in eight-atom primatives. */
  int getNZ() const {return nz;}
  /** Get lattice constant. */
  double getA() const {return a;}
  /** Get the diagonal strain elements. */
  double getStrain(const int idim) const {return (idim==0)?strainxx:(idim==1)
                                                 ?strainyy:strainzz;}
  /** Get pointer to bulk material. */
  const Material* getBulkMaterial() const {return bulkMaterial;}
protected:
  /** X-Extent of the supercell in eight-atom primatives. */
  int nx;
  /** Y-Extent of the supercell in eight-atom primatives. */
  int ny;
  /** Z-Extent of the supercell in eight-atom primatives. */
  int nz;
  /** Lattice constant. */
  double a;
  /** Diagonal strain elements. */
  double strainxx, strainyy, strainzz;
  /** Bulk material. */
  const Material* bulkMaterial;
};
#endif
