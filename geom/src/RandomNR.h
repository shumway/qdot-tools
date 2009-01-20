// $Id: RandomNR.h,v 1.2 2006/07/05 17:54:37 jshumwa Exp $
/*
    Copyright (C) 2004 Gabriel Bester and John B. Shumway, Jr.

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
#ifndef __RandomNR_h_
#define __RandomNR_h_

/* 
Routines to reproduce random number generator behavior 
described in Numerical Recipies routine ran1.
"Minimal" random number generator of Park and Miller with Bays-Durham
shuffle and added safeguards.
@author Gabriel Bester and John Shumway
@version $Revision: 1.2 $
*/
class RandomNR{
public:
  /// Construct by seeding the random number generator.
  RandomNR(int seed=0);
  /// Seed the random number generator.
  void sran1(int seed);
  /// Get a uniform random deviate between 0.0 and 1.0 (excluding endpoints).
  double ran1();
private:
  /// RNMX should approximate  the largest floating value less than 1.
  static const double RNMX;
  /// 7^5
  static const int IA=16807;
  /// 2*31-1
  static const int IM=2147483647; 
  static const double AM;
  static const int IQ=IM/IA;
  static const int IR=IM-IA*IQ;
  /// Number of divsions of integers into 32 shuffling bins.
  static const int NDIV=(1+(IM-1)/32);
  /// Array for shuffling.
  int shuffle[32];
  /// Pointer for shuffled array.
  int ishuffle;
  /// Stored integer for random number generation.
  int irand;
};
#endif
