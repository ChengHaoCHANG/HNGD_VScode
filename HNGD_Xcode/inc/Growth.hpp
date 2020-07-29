#ifndef Growth_hpp
#define Growth_hpp

#include <stdio.h>
#include "Precipitation.hpp"
#include "PhysicsConstants.h"

class Growth : public Precipitation
{
public:
    Growth(Sample* sample, double Kmob0, double Kth0, double Eg, double p);
    void computeKinetics()  ;
    void computeDrivForce() ;
    
private:
    const double _Kth0 ;
    const double _Kmob0 ;
    const double _Eg ;
    const double _p ;
    const vector<double> * _temperature ;
    const vector<double> * _tssd ;
    const vector<double> * _Cprec ;
    const vector<double> * _Ctot ;
};


#endif /* Growth_hpp */
