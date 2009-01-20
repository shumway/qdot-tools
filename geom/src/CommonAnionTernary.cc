#include "CommonAnionTernary.h"
#include <cmath>
#include <stdlib.h>
#include "RandomNR.h"
CommonAnionTernary::CommonAnionTernary(const std:: string& name, 
          const std::string& anionName, const std::string& cation1Name,
          const std::string& cation2Name, const double x1, const int seed,
          RandomNR& rand)
  : Material(name,3), x1(x1), rand(rand) {
  speciesName[0]=anionName;
  speciesName[1]=cation1Name;
  speciesName[2]=cation2Name;
  speciesID[0]=1;
  speciesID[1]=2;
  speciesID[2]=3;
  //if(seed) srand48(seed);        // Initialize the random number generator
  if(seed) rand.sran1(seed);
}

int CommonAnionTernary::getAnionIndex(const Vec3& p) const {
  return speciesID[0];
}

int CommonAnionTernary::getCationIndex(const Vec3& p) const {
  //return (drand48()<x1) ? speciesID[1] : speciesID[2];
  return (rand.ran1()<=x1) ? speciesID[1] : speciesID[2];
}
