// $Id: geom.cc,v 1.4 2004/06/10 19:41:32 jshumwa Exp $
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
#include <iostream>
#include "InputParser.h"
#include "CommentH5.h"

/** @mainpage Atomic Nanostructure Layout Program
@section Overview
      
<p>This program, geom, is a simple utility for turning nanostructure 
descriptions into atomic layout. These atomic positions are common input 
for several types of quantum dot models, such as tight-binding or empirical 
pseudopotentials. The user specifies the size, shape, and composition of a 
coherent III-V or group IV semiconductor heterostructure. Shapes include 
lens, pyramidal, and conical dots, as well as quantum wells, and may be 
filled with different alloy or pure materials. Output is an HDF5 file,
struct.h5, containing atom positions, chemical identities, lattice indexing,
and neighbor tables. The object-oriented design makes it easy for users to 
add new shapes and materials to this utility.</p>

The Nanostructure is described by the types of data:
<ol>
<li>Simulation data, such as SuperCell dimensions, lattice constants,
and embedding material.</li>
<li>Material descriptons, such as alloys, pure materials, or even vacuum
that make up the nanostructure.</li>
<li>Structure shape descriptions, such as dots and quantum wells, 
that can be combined to make composite nanostructures.</li>
</ol>

@section sec1 Getting started
See the InputParser documentation for a sample input file.
@section sec2 Extension capabilities
@author John Shumway
*/
int main(int argc, char** argv) {
  const std::string filename = (argc==2) ? argv[1] : "struct.xml";
  std::cout << std::endl;
  std::cout << "Nanostructure geometry program: " << filename << std::endl;
  std::cout << std::endl;
  CommentH5 comment(filename);
  // Parser xml input file.
  InputParser parser(filename.c_str());
  parser.parse();
  comment.h5Write("struct.h5");
}

