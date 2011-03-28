// $Id: InputParser.cc,v 1.9 2007/01/23 17:20:35 jshumwa Exp $
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
#include "InputParser.h"
#include <iostream>
#include <fstream>
#include "SuperCell.h"
#include "Nanostructure.h"
#include "Well.h"
#include "Dome.h"
#include "Hut.h"
#include "Sphere.h"
#include "Lens.h"
#include "Cone.h"
#include "Ellipse.h"
#include "Ellipsoid.h"
#include "Cigar.h"
#include "Cylinder.h"
#include "LongLens.h"
#include "Ring.h"
#include "Pyramid.h"
#include "StructH5Writer.h"
#include "Material.h"
#include "Binary.h"
#include "Pure.h"
#include "BinaryAlloy.h"
#include "CommonAnionTernary.h"
#include "CommonAnionLinear.h"
#include "CommonCationTernary.h"
#include "RadialAlloy.h"
#include "RandomNR.h"

InputParser::InputParser(const char* filename){
  doc = xmlParseFile(filename);
}

InputParser::~InputParser() {
  xmlFreeDoc(doc);
}

void InputParser::parse() {

  SuperCell* cell=0;
  Nanostructure ns;
  int etchIndex=0;
  bool reconstruct=false;
  char * property_ptr;
  RandomNR rand;

  for (xmlNodePtr node = xmlDocGetRootElement(doc)->children; 
       node != NULL; node = node->next) {
    std::string name((char*)node->name);
    std::cout << name << std::endl;
    if (name=="Materials") {
      for (xmlNodePtr matNode=node->children; matNode != NULL;
           matNode = matNode->next) {
        std::string name((char*)matNode->name);
        std::cout << "  " << name << std::endl;
        if (name=="Binary") {
          std::string name((char*)xmlGetProp(matNode,(xmlChar*)"name"));
          std::string anion((char*)xmlGetProp(matNode,(xmlChar*)"anion"));
          std::string cation((char*)xmlGetProp(matNode,(xmlChar*)"cation"));
          ns.addMaterial(new Binary(name,anion,cation));
        } else
        if (name=="Pure") {
          std::string name((char*)xmlGetProp(matNode,(xmlChar*)"name"));
          ns.addMaterial(new Pure(name));
        } else
        if (name=="BinaryAlloy") {
          std::string name((char*)xmlGetProp(matNode,(xmlChar*)"name"));
          std::string atom1((char*)xmlGetProp(matNode,(xmlChar*)"atom1"));
          std::string atom2((char*)xmlGetProp(matNode,(xmlChar*)"atom2"));
          double x1 = atof((char*)xmlGetProp(matNode,(xmlChar*)"x1"));
          ns.addMaterial(new BinaryAlloy(name,atom1,atom2,x1));
        } else
        if (name=="CommonAnionTernary") {
          std::string name((char*)xmlGetProp(matNode,(xmlChar*)"name"));
          std::string anion((char*)xmlGetProp(matNode,(xmlChar*)"anion"));
          std::string cation1((char*)xmlGetProp(matNode,(xmlChar*)"cation1"));
          std::string cation2((char*)xmlGetProp(matNode,(xmlChar*)"cation2"));
          double x1 = atof((char*)xmlGetProp(matNode,(xmlChar*)"x1"));
          int seed = 0;
          property_ptr = (char *) xmlGetProp(matNode,(xmlChar*)"seed");
          if(property_ptr) seed =  atoi(property_ptr);
          ns.addMaterial(new CommonAnionTernary(name,anion,cation1,cation2,
                                                x1,seed,rand));
        } else
        if (name=="CommonCationTernary") {
          std::string name((char*)xmlGetProp(matNode,(xmlChar*)"name"));
          std::string anion1((char*)xmlGetProp(matNode,(xmlChar*)"anion1"));
          std::string anion2((char*)xmlGetProp(matNode,(xmlChar*)"anion2"));
          std::string cation((char*)xmlGetProp(matNode,(xmlChar*)"cation"));
          double x1 = atof((char*)xmlGetProp(matNode,(xmlChar*)"x1"));
          int seed = 0;
          property_ptr = (char *) xmlGetProp(matNode,(xmlChar*)"seed");
          if(property_ptr) seed =  atoi(property_ptr);
          ns.addMaterial(new CommonCationTernary(name,anion1,anion2,
                                                 cation,x1,seed,rand));
        } else
        if (name=="CommonAnionLinear") {
          std::string name((char*)xmlGetProp(matNode,(xmlChar*)"name"));
          std::string anion((char*)xmlGetProp(matNode,(xmlChar*)"anion"));
          std::string cation1((char*)xmlGetProp(matNode,(xmlChar*)"cation1"));
          std::string cation2((char*)xmlGetProp(matNode,(xmlChar*)"cation2"));
          double x1Bottom=atof((char*)xmlGetProp(matNode,(xmlChar*)"x1Bottom"));
          double x1Top = atof((char*)xmlGetProp(matNode,(xmlChar*)"x1Top"));
          double zBottom = atof((char*)xmlGetProp(matNode,(xmlChar*)"zBottom"));
          double zTop = atof((char*)xmlGetProp(matNode,(xmlChar*)"zTop"));
          int seed = 0;
          property_ptr = (char *) xmlGetProp(matNode,(xmlChar*)"seed");
          if(property_ptr) seed =  atoi(property_ptr);
          ns.addMaterial(new CommonAnionLinear(name,anion,cation1,cation2,
                               x1Bottom,x1Top,zBottom,zTop,seed,rand));
        } else
        if (name=="RadialAlloy") {
          std::string name((char*)xmlGetProp(matNode,(xmlChar*)"name"));
          std::string anion((char*)xmlGetProp(matNode,(xmlChar*)"anion"));
          std::string cation1((char*)xmlGetProp(matNode,(xmlChar*)"cation1"));
          std::string cation2((char*)xmlGetProp(matNode,(xmlChar*)"cation2"));
          double x1_0 = atof((char*)xmlGetProp(matNode,(xmlChar*)"x1_0"));
          double x1_r = atof((char*)xmlGetProp(matNode,(xmlChar*)"x1_r"));
          double r = atof((char*)xmlGetProp(matNode,(xmlChar*)"r"));
          int seed = 0;
          property_ptr = (char *) xmlGetProp(matNode,(xmlChar*)"seed");
          if(property_ptr) seed =  atoi(property_ptr);
          ns.addMaterial(new RadialAlloy(name,anion,cation1,cation2,
                                         x1_0,x1_r,r,seed,rand));
        }
      }
    } else
    if (name=="SuperCell") {
      int nx = atoi((char*)xmlGetProp(node,(xmlChar*)"nx"));
      int ny = atoi((char*)xmlGetProp(node,(xmlChar*)"ny"));
      int nz = atoi((char*)xmlGetProp(node,(xmlChar*)"nz"));
      double a = atof((char*)xmlGetProp(node,(xmlChar*)"a"));
      double strainxx=0, strainyy=0, strainzz=0;
      char* value=(char*)xmlGetProp(node,(xmlChar*)"strainxx");
      if (value) strainxx=atof(value);
      value=(char*)xmlGetProp(node,(xmlChar*)"strainyy");
      if (value) strainyy=atof(value);
      value=(char*)xmlGetProp(node,(xmlChar*)"strainzz");
      if (value) strainzz=atof(value);
      std::string matName((char*)xmlGetProp(node,(xmlChar*)"material"));
      const Material* material = ns.getMaterial(matName);
      cell = new SuperCell(nx,ny,nz,a,strainxx,strainyy,strainzz,material);
    } else
    if (name=="Etch") {
      std::string specName((char*)xmlGetProp(node,(xmlChar*)"species"));
      etchIndex=ns.getSpeciesIndex(specName);
    } else
    if (name=="Reconstruct") {
      reconstruct=true;
    } else 
    if (name=="Structures") {
      for (xmlNodePtr structNode=node->children; structNode != NULL;
           structNode = structNode->next) {
        std::string name((char*)structNode->name);
        std::cout << "  " << name << std::endl;
        if (name=="Well") {
          double z = atof((char*)xmlGetProp(structNode,(xmlChar*)"z"));
          double thickness =
                   atof((char*)xmlGetProp(structNode,(xmlChar*)"thickness"));
          std::string matName((char*)xmlGetProp(structNode,
                                                (xmlChar*)"material"));
          const Material* material = ns.getMaterial(matName);
          ns.addStructure(new Well(z,thickness,material));
        } else
        if (name=="Sphere") {
          double x = atof((char*)xmlGetProp(structNode,(xmlChar*)"x"));
          double y = atof((char*)xmlGetProp(structNode,(xmlChar*)"y"));
          double z = atof((char*)xmlGetProp(structNode,(xmlChar*)"z"));
          double r = atof((char*)xmlGetProp(structNode,(xmlChar*)"radius"));
          std::string matName((char*)xmlGetProp(structNode,
                                                (xmlChar*)"material"));
          const Material* material = ns.getMaterial(matName);
          ns.addStructure(new Sphere(x,y,z,r,material));
        } else
        if (name=="Cone") {
          double x = atof((char*)xmlGetProp(structNode,(xmlChar*)"x"));
          double y = atof((char*)xmlGetProp(structNode,(xmlChar*)"y"));
          double z = atof((char*)xmlGetProp(structNode,(xmlChar*)"z"));
          double height=atof((char*)xmlGetProp(structNode,(xmlChar*)"height"));
          double base = atof((char*)xmlGetProp(structNode,(xmlChar*)"base"));
          double top = atof((char*)xmlGetProp(structNode,(xmlChar*)"top"));
          std::string matName((char*)xmlGetProp(structNode,
                                                (xmlChar*)"material"));
          const Material* material = ns.getMaterial(matName);
          ns.addStructure(new Cone(x,y,z,height,base,top,material));
        } else
        if (name=="Ellipsoid") {
          double x = atof((char*)xmlGetProp(structNode,(xmlChar*)"x"));
          double y = atof((char*)xmlGetProp(structNode,(xmlChar*)"y"));
          double z = atof((char*)xmlGetProp(structNode,(xmlChar*)"z"));
          double xr=atof((char*)xmlGetProp(structNode,(xmlChar*)"xr"));
          double yr = atof((char*)xmlGetProp(structNode,(xmlChar*)"yr"));
          double zr = atof((char*)xmlGetProp(structNode,(xmlChar*)"zr"));
          std::string matName((char*)xmlGetProp(structNode,
                                                (xmlChar*)"material"));
          const Material* material = ns.getMaterial(matName);
          ns.addStructure(new Ellipsoid(x,y,z,xr,yr,zr,material));
        } else
        if (name=="Cigar") {
          double x = atof((char*)xmlGetProp(structNode,(xmlChar*)"x"));
          double y = atof((char*)xmlGetProp(structNode,(xmlChar*)"y"));
          double z = atof((char*)xmlGetProp(structNode,(xmlChar*)"z"));
          double radius=atof((char*)xmlGetProp(structNode,(xmlChar*)"radius"));
          double length = atof((char*)xmlGetProp(structNode,
                                                 (xmlChar*)"length"));
          double vx = atof((char*)xmlGetProp(structNode,(xmlChar*)"vx"));
          double vy = atof((char*)xmlGetProp(structNode,(xmlChar*)"vy"));
          double vz = atof((char*)xmlGetProp(structNode,(xmlChar*)"vz"));
          std::string matName((char*)xmlGetProp(structNode,
                                                (xmlChar*)"material"));
          const Material* material = ns.getMaterial(matName);
          ns.addStructure(new Cigar(x,y,z,radius,length,vx,vy,vz,material));
        } else
        if (name=="Cylinder") {
          double x = atof((char*)xmlGetProp(structNode,(xmlChar*)"x"));
          double y = atof((char*)xmlGetProp(structNode,(xmlChar*)"y"));
          double r = atof((char*)xmlGetProp(structNode,(xmlChar*)"radius"));
          std::string matName((char*)xmlGetProp(structNode,
                                                (xmlChar*)"material"));
          const Material* material = ns.getMaterial(matName);
          ns.addStructure(new Cylinder(x,y,r,material));
        } else
        if (name=="Lens") {
          double x = atof((char*)xmlGetProp(structNode,(xmlChar*)"x"));
          double y = atof((char*)xmlGetProp(structNode,(xmlChar*)"y"));
          double z = atof((char*)xmlGetProp(structNode,(xmlChar*)"z"));
          double d = atof((char*)xmlGetProp(structNode,(xmlChar*)"diameter"));
          double h = atof((char*)xmlGetProp(structNode,(xmlChar*)"height"));
          std::string matName((char*)xmlGetProp(structNode,
                                                (xmlChar*)"material"));
          const Material* material = ns.getMaterial(matName);
          ns.addStructure(new Lens(x,y,z,h,d,material));
        }
        if (name=="LongLens") {
          double x = atof((char*)xmlGetProp(structNode,(xmlChar*)"x"));
          double y = atof((char*)xmlGetProp(structNode,(xmlChar*)"y"));
          double z = atof((char*)xmlGetProp(structNode,(xmlChar*)"z"));
          double d1 = atof((char*)xmlGetProp(structNode,(xmlChar*)"majoraxis"));
          double d2 = atof((char*)xmlGetProp(structNode,(xmlChar*)"minoraxis"));
          double h = atof((char*)xmlGetProp(structNode,(xmlChar*)"height"));
          std::string matName((char*)xmlGetProp(structNode,
                                                (xmlChar*)"material"));
          const Material* material = ns.getMaterial(matName);
          ns.addStructure(new LongLens(x,y,z,h,d1,d2,material));
        } else
        if (name=="Pyramid") {
          double x = atof((char*)xmlGetProp(structNode,(xmlChar*)"x"));
          double y = atof((char*)xmlGetProp(structNode,(xmlChar*)"y"));
          double z = atof((char*)xmlGetProp(structNode,(xmlChar*)"z"));
          double h = atof((char*)xmlGetProp(structNode,(xmlChar*)"height"));
          double b = atof((char*)xmlGetProp(structNode,(xmlChar*)"base"));
          double t = atof((char*)xmlGetProp(structNode,(xmlChar*)"top"));
          std::string matName((char*)xmlGetProp(structNode,
                                                (xmlChar*)"material"));
          const Material* material = ns.getMaterial(matName);
          ns.addStructure(new Pyramid(x,y,z,h,b,t,material));
        }
        if (name=="Hut") {
          double x = atof((char*)xmlGetProp(structNode,(xmlChar*)"x"));
          double y = atof((char*)xmlGetProp(structNode,(xmlChar*)"y"));
          double z = atof((char*)xmlGetProp(structNode,(xmlChar*)"z"));
          double h = atof((char*)xmlGetProp(structNode,(xmlChar*)"height"));
          double xb = atof((char*)xmlGetProp(structNode,(xmlChar*)"xbase"));
          double yb = atof((char*)xmlGetProp(structNode,(xmlChar*)"ybase"));
          double xt = atof((char*)xmlGetProp(structNode,(xmlChar*)"xtop"));
          double yt = atof((char*)xmlGetProp(structNode,(xmlChar*)"ytop"));
          std::string matName((char*)xmlGetProp(structNode,
                                                (xmlChar*)"material"));
          const Material* material = ns.getMaterial(matName);
          ns.addStructure(new Hut(x,y,z,h,xb,yb,xt,yt,material));
        }
        if (name=="Dome") {
          double x = atof((char*)xmlGetProp(structNode,(xmlChar*)"x"));
          double y = atof((char*)xmlGetProp(structNode,(xmlChar*)"y"));
          double z = atof((char*)xmlGetProp(structNode,(xmlChar*)"z"));
          std::string matName((char*)xmlGetProp(structNode,
                                                (xmlChar*)"material"));
          const Material* material = ns.getMaterial(matName);
          Dome::FArray facets;
          for (xmlNodePtr facetNode=structNode->children; facetNode != NULL;
                          facetNode = facetNode->next) {
            std::string name((char*)facetNode->name);
            if (name=="Facet") {
              int i = atoi((char*)xmlGetProp(facetNode,(xmlChar*)"i"));
              int j = atoi((char*)xmlGetProp(facetNode,(xmlChar*)"j"));
              int k = atoi((char*)xmlGetProp(facetNode,(xmlChar*)"k"));
              double h = atof((char*)xmlGetProp(facetNode,(xmlChar*)"height"));
              facets.push_back(Dome::Facet(i,j,k,h));
            }
          }
          ns.addStructure(new Dome(x,y,z,facets,material));
        }
        if (name=="Ring" || name=="HalfRing") {
          double x = atof((char*)xmlGetProp(structNode,(xmlChar*)"x"));
          double y = atof((char*)xmlGetProp(structNode,(xmlChar*)"y"));
          double z = atof((char*)xmlGetProp(structNode,(xmlChar*)"z"));
          double h = atof((char*)xmlGetProp(structNode,(xmlChar*)"height"));
          double rout = atof((char*)xmlGetProp(structNode,(xmlChar*)"radiusin"));
          double rin = atof((char*)xmlGetProp(structNode,(xmlChar*)"radiusout"));
          double s11 = atof((char*)xmlGetProp(structNode,(xmlChar*)"stretch"));
          if (s11==0) s11=1;
          std::string matName((char*)xmlGetProp(structNode,
                                                (xmlChar*)"material"));
          const Material* material = ns.getMaterial(matName);
          ns.addStructure(new Ring(x,y,z,h,rin,rout,s11,material));
        }
        if (name=="Ellipse") {
          double x = atof((char*)xmlGetProp(structNode,(xmlChar*)"x"));
          double y = atof((char*)xmlGetProp(structNode,(xmlChar*)"y"));
          double z = atof((char*)xmlGetProp(structNode,(xmlChar*)"z"));
          double h = atof((char*)xmlGetProp(structNode,(xmlChar*)"height"));
          double a = atof((char*)xmlGetProp(structNode,(xmlChar*)"r_x"));
          double a2 = atof((char*)xmlGetProp(structNode,(xmlChar*)"r_y"));
          double r2 = atof((char*)xmlGetProp(structNode,(xmlChar*)"width_x"));
          double r = atof((char*)xmlGetProp(structNode,(xmlChar*)"width_y"));
          double phi = atof((char*)xmlGetProp(structNode,(xmlChar*)"phi"));
          std::string matName((char*)xmlGetProp(structNode,
                                                (xmlChar*)"material"));
          const Material* material = ns.getMaterial(matName);
          ns.addStructure(new Ellipse(x,y,z,phi,h,a,a2,r2,r,material));
        }
      }
    }
  }

  if (cell) {
    StructH5Writer writer(*cell,ns);
    if (etchIndex) writer.etch(etchIndex);
    if (reconstruct) writer.reconstruct();
    writer.write("struct.h5");
  } else {
    std::cout << "Supercell is undefined.  Nothing written." << std::endl;
  }
  delete cell;
}
