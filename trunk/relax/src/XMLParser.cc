// $Id: XMLParser.cc,v 1.3 2007/03/14 19:47:50 jshumwa Exp $
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
#include "XMLParser.h"

double XMLParser::getDoubleAttribute(const xmlNodePtr& node,
                                     const std::string& attName) {
  char* temp = (char*)xmlGetProp(node,BAD_CAST attName.c_str());
  double value = temp?atof((char*)temp):0;
  free(temp);
  return value;
}

int XMLParser::getIntAttribute(const xmlNodePtr& node,
                               const std::string& attName) {
  char* temp = (char*)xmlGetProp(node,BAD_CAST attName.c_str());
  int value = temp?atoi((char*)temp):0;
  free(temp);
  return value;
}

std::string XMLParser::getName(const xmlNodePtr& node) {
  std::string value((char*)node->name);
  return value;
}

bool XMLParser::getBoolAttribute(const xmlNodePtr& node,
                                          const std::string& attName) {
  char* temp = (char*)xmlGetProp(node,BAD_CAST attName.c_str());
  std::string value(temp?(char*)temp:"");
  free(temp);
  return value!="" && (value[0]=='t' || value[0]=='T');
}

std::string XMLParser::getStringAttribute(const xmlNodePtr& node,
                                          const std::string& attName) {
  char* temp = (char*)xmlGetProp(node,BAD_CAST attName.c_str());
  std::string value(temp?(char*)temp:"");
  free(temp);
  return value;
}

#ifdef ENABLE_FLOAT
const blitz::TinyVector<float,3> 
#else
const blitz::TinyVector<double,3> 
#endif
  XMLParser::getVecAttribute(const xmlNodePtr& node) {
#ifdef ENABLE_FLOAT
  blitz::TinyVector<float,3> v;
#else
  blitz::TinyVector<double,3> v;
#endif
  v[0]=getDoubleAttribute(node,"x");
  v[1]=getDoubleAttribute(node,"y");
  v[2]=getDoubleAttribute(node,"z");
  return v;
}

std::string XMLParser::getLinkBase(const xmlNodePtr& node,
                                   const std::string& attName) {
  char* temp = (char*)xmlGetProp(node,BAD_CAST attName.c_str());
  std::string href(temp?(char*)temp:"");
  free(temp);
  return href.substr(0,href.find('|'));
}

std::string XMLParser::getLinkPath(const xmlNodePtr& node,
                                   const std::string& attName) {
  char* temp = (char*)xmlGetProp(node,BAD_CAST attName.c_str());
  std::string href(temp?(char*)temp:"");
  free(temp);
  return href.substr(href.find('|')+1);
}
