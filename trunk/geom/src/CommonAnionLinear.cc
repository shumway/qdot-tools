#include "CommonAnionLinear.h"
#include "Vec3.h"
#include <cmath>
#include <stdlib.h>
#include "RandomNR.h"
CommonAnionLinear::CommonAnionLinear(const  std::string& name, 
          const  std::string& anionName, const  std::string& cation1Name,
          const  std::string& cation2Name, const double x1Bottom,
          const double x1Top, const double zBottom, const double zTop, 
          const int seed, RandomNR& rand)
  : Material(name,3), 
    x1Bottom(x1Bottom), x1Top(x1Top), zBottom(zBottom), zTop(zTop), rand(rand) {
  speciesName[0]=anionName;
  speciesName[1]=cation1Name;
  speciesName[2]=cation2Name;
  speciesID[0]=1;
  speciesID[1]=2;
  speciesID[2]=3;
  //if(seed) srand48(seed);        // Initialize the random number generator
  if(seed) rand.sran1(seed);
}

int CommonAnionLinear::getAnionIndex(const Vec3& p) const {
  return speciesID[0];
}

int CommonAnionLinear::getCationIndex(const Vec3& p) const {
  double x1 = (p.z<=zBottom)?x1Bottom:(p.z>=zTop)?zTop:
              x1Bottom+(p.z-zBottom)/(zTop-zBottom)*(x1Top-x1Bottom);
  //return (drand48()<x1) ? speciesID[1] : speciesID[2];
  return (rand.ran1()<=x1) ? speciesID[1] : speciesID[2];
}
