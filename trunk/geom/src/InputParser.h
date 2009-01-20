// $Id: InputParser.h,v 1.3 2004/06/10 19:41:32 jshumwa Exp $
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
#ifndef __InputParser_h_
#define __InputParser_h_

#include <libxml/tree.h>
#include <libxml/parser.h>

/** XML based input parser.  Uses gnu library libxml (http://xmlsoft.org)
  for the parsing.

    This class is not neccesary to use the objects in the geom
    program.  You can hard-code a geometry program that constructs
    the structure, supercell, and material objects directly, with
    no need to use XML (look at InputParser source code).

    This code needs exception handling!  Most errors in the XML
    input currently result in core dumps rather than useful error
    messages.
    
    Sample XML input file:
    \verbatim
    <Nanostructure>
      <Materials>
        <Binary name="GaAs" anion="As" cation="Ga"/>
        <Binary name="InAs" anion="As" cation="In"/>
      </Materials>
      <SuperCell nx="80" ny="80" nz="40" a="10.68" material="GaAs"/>
      <Structures>
        <Well z="19.75" thickness="0.5" material="InAs"/>
        <!-- 252A x 35A lens shaped dot. -->
        <Lens x="40" y="40" z="20" height="6.25" 
              diameter="44.6" material="InAs"/>
      </Structures>
    </Nanostructure>
    \endverbatim

    The syntax for many of the XML elements mattches the class constructors
    (``see also'' classes, listed below).
    @see SuperCell
    @see Binary, CommonAnionTernary, CommonCationTernary
    @see Well, Sphere, Lens
    @todo Add exeption handling.
    @author John Shumway */
class InputParser{
public:
  /** Open an xml file for parsing. */
  InputParser(const char* filename);
  /** Free the xml file and parser resources. */
  ~InputParser();
  /** Parse the xml file. */
  void parse();
protected:
  /** Pointer the xml document. */
  xmlDocPtr doc;
};
#endif
