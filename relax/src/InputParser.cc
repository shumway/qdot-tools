// $Id: InputParser.cc,v 1.9 2007/12/28 02:19:40 jshumwa Exp $
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "InputParser.h"
#include "AtomicPotential.h"
#include "Checkpoint.h"
#include "StillWeb.h"
#include "VFF.h"
#include "CellRelaxer.h"

InputParser::InputParser(const std::string& filename)
  : context (0), annealSweeps(0) {
  doc = xmlParseFile((char*)filename.c_str());
  context = xmlXPathNewContext(doc);
}

InputParser::~InputParser() {
  xmlXPathFreeContext(context);
  xmlFreeDoc(doc);
}

void InputParser::parse(const xmlXPathContextPtr& ctxt) {
  // Read the input filename.
  xmlXPathObjectPtr obj = xmlXPathEval(BAD_CAST"//Structure",ctxt);
  structInFileName=getStringAttribute(obj->nodesetval->nodeTab[0],"infile");
  // Check if relaxation requested.
  obj = xmlXPathEval(BAD_CAST"//Relax",context);
  if (obj->nodesetval->nodeNr>0) {
    xmlNodePtr node=obj->nodesetval->nodeTab[0];
    maxiter=getIntAttribute(node,"maxiter");
    ftol=getDoubleAttribute(node,"ftol");
    scale=getDoubleAttribute(node,"scale");
    if (scale==0.0) scale=1.0;
    relaxedFilename=getStringAttribute(node,"outfile");
  } else {
    maxiter=0;
  }
  // Check if stress calculation is requested.
  obj = xmlXPathEval(BAD_CAST"//CalculateStresses",context);
  if (obj->nodesetval->nodeNr>0) {
    xmlNodePtr node=obj->nodesetval->nodeTab[0];
    stressFilename=getStringAttribute(node,"outfile");
  }
  // Check if stress calculation is requested.
  obj = xmlXPathEval(BAD_CAST"//CalculateSurfaceStress",context);
  if (obj->nodesetval->nodeNr>0) {
    xmlNodePtr node=obj->nodesetval->nodeTab[0];
    surfaceFilename=getStringAttribute(node,"outfile");
  }
}

Checkpoint* InputParser::parseCheckpoint(CoordsH5& coords, 
    SpeciesH5& species, NeighborsH5& neighbors, SuperCellH5& cell) {
  // Read the input filename.
  xmlXPathObjectPtr obj = xmlXPathEval(BAD_CAST"//Checkpoint",context);
  if (obj->nodesetval->nodeNr==0) return 0;
  xmlNodePtr node=obj->nodesetval->nodeTab[0];
  std::string filename=getStringAttribute(node,"file");
  int interval=getIntAttribute(node,"interval");
  return new Checkpoint(coords,species,neighbors,cell,filename,interval);
}

AtomicPotential* InputParser::parsePotential() {
  AtomicPotential* pot=0;
  xmlXPathObjectPtr obj = xmlXPathEval(BAD_CAST"//ForceModel/*",context);
  if (obj->nodesetval->nodeNr==0) return 0;
  xmlNodePtr node=obj->nodesetval->nodeTab[0];
  std::string name=getName(node);
  if (name=="StillingerWeber") {
    pot = new StillWeb();
  } else  if (name=="VFF") {
    pot = new VFF();
  }
  return pot;
}

CellRelaxer* InputParser::parseCellRelaxer(CoordsH5 &coords, NeighborsH5 &nbr,
     SpeciesH5 &species, SuperCellH5 &cell, AtomicPotential &pot, 
     TotalEnergy &en) {
  CellRelaxer* relaxer=0;
  xmlXPathObjectPtr obj = xmlXPathEval(BAD_CAST"//Relax/CellRelax",context);
  if (obj->nodesetval->nodeNr==0) return 0;
  xmlNodePtr node=obj->nodesetval->nodeTab[0];
  cellRelaxInterval=getIntAttribute(node,"interval");
  std::string modeName=getStringAttribute(node,"mode");
  const int mode = (modeName=="volume") ? CellRelaxer::VOLUME
                 : (modeName=="zonly") ? CellRelaxer::ZONLY : 0;
  relaxer=new CellRelaxer(coords,nbr,species,cell,pot,en,mode);
  return relaxer;
}
