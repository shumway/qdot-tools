#include "RadialAlloy.h"
#include "Vec3.h"
#include <cmath>
#include <stdlib.h>
#include "RandomNR.h"
RadialAlloy::RadialAlloy(const  std::string& name, 
          const  std::string& anionName, const  std::string& cation1Name,
          const  std::string& cation2Name, const double x1_0, 
          const double x1_r, const double r, const int seed, RandomNR& rand)
  : CommonAnionTernary(name,anionName,cation1Name,cation2Name,0,seed, rand),
    x1_0(x1_0), x1_r(x1_r), r(r), rand(rand) {
}

int RadialAlloy::getAnionIndex(const Vec3& p) const {
  return speciesID[0];
}

int RadialAlloy::getCationIndex(const Vec3& p) const {
  double rho=sqrt(p.x*p.x+p.y*p.y);
  double x1=(rho>r)?x1_r:x1_0+rho*(x1_r-x1_0)/r;
  //return (drand48()<x1) ? speciesID[1] : speciesID[2];
  return (rand.ran1()<=x1) ? speciesID[1] : speciesID[2];
}
