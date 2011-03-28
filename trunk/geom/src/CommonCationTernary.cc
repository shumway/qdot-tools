#include "CommonCationTernary.h"
#include <cmath>
#include <stdlib.h>
#include "RandomNR.h"
CommonCationTernary::CommonCationTernary(const  std::string& name, 
          const  std::string& anion1Name, const  std::string& anion2Name,
          const  std::string& cationName, const double x1, const int seed,
	  RandomNR& rand)
  : Material(name,3), x1(x1), rand(rand) {
  speciesName[0]=anion1Name;
  speciesName[1]=anion2Name;
  speciesName[2]=cationName;
  speciesID[0]=1;
  speciesID[1]=2;
  speciesID[2]=3;
  //if(seed) srand48(seed); // Initialize the random number generator
  if(seed) rand.sran1(seed); // Initialize the random number generator
}

int CommonCationTernary::getAnionIndex(const Vec3& p) const {
  //return (drand48()<x1) ? speciesID[0] : speciesID[1];
  return (rand.ran1()<=x1) ? speciesID[0] : speciesID[1];
}

int CommonCationTernary::getCationIndex(const Vec3& p) const {
  return speciesID[2];
}
