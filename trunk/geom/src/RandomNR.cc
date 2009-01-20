// $Id: RandomNR.cc,v 1.3 2006/07/05 17:54:37 jshumwa Exp $
/*
    Copyright (C) 2004 Gabriel Bester and John Shumway`

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

#include "RandomNR.h"
#include <iostream>

RandomNR::RandomNR(int iseed) {
  sran1(iseed);
}

void RandomNR::sran1(int iseed){
  if (iseed < 0 ) iseed = -iseed;
  std::cout << "initialising with seed = " <<  iseed << std::endl;
  irand = iseed ? iseed : 1; /* avoid seed = 0*/
  for (int i = 32+7; i >=0 ; i--) {
     /* load shuffle table-- after 8 warm-ups */
     int is = irand/IQ;
     int it = irand-is*IQ;
     irand = IA*it - is*IR;
     if (irand < 0 ) irand += IM;
     shuffle[i%32]=irand;
  }
  ishuffle = shuffle[0];
}

double RandomNR::ran1() {
  /// compute irand - mod(IA*irand,IM) without overflow 
  int is = irand/IQ;
  int it = irand-is*IQ;
  irand = IA*it - IR*is;
  if (irand < 0 ) irand += IM;
  int jshuffle = ishuffle/NDIV;
  ishuffle = shuffle[jshuffle];
  shuffle[jshuffle] = irand;
  double randnew = AM*ishuffle ;
  return (randnew <  RNMX)? randnew : RNMX;
}

const double RandomNR::RNMX=(1.-1.e-14);
const double RandomNR::AM=(1./IM);
