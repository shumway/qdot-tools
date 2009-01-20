// $Id: EtchedLattice.h,v 1.1 2004/06/11 00:48:52 jshumwa Exp $
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
#ifndef __EtchedLattice_h_
#define __EtchedLattice_h_

#include <string>
#include "LatticeUtil.h"
class SuperCell;
class Nanostructure;
class LatticeIndexingH5;

/** Utility to etch dummy atoms from structure for free-standing or
    surface nanostructures.
    @author John Shumway */
class EtchedLattice : public LatticeUtil {
public:
  /** Constructor. */
  EtchedLattice(const LatticeUtil&, const int etchedAtomID);
private:
};
#endif
