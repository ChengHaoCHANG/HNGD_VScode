#ifndef Nucleation_hpp
#define Nucleation_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>
#include "Precipitation.hpp"
#include "PhysicsConstants.h"

using namespace std ;

class Nucleation : public Precipitation
{
public:
    Nucleation(Sample* sample, double Kn0);
    void computeKinetics()  ;
    void computeDrivForce() ;

private:
    const double _Kn0 ;
    
    const double _gamma ;
    const double _v ;
    double volume_energy(double T) ;
    
    
};

#endif /* Nucleation_hpp */
