#include <iostream>

#include "Dissolution.hpp"

Dissolution :: Dissolution(Sample* sample, double Kd0, double Ed) :
    Mechanism   (sample),

    _Kd0        (Kd0),
    _Ed         (Ed),

    _temperature(& (sample->returnTemperature())),
    _tssd       (& (sample->returnTSSd())),
    _Css        (& (sample->returnSolutionContent())),
    _Ctot       (& (sample->returnTotalContent()))
{}

void Dissolution :: computeKinetics()
{
    // Dissolution kinetics follow an Arrhenius law
    for(int k=0; k<_nbCells; k++)
        _kinetic_factor[k] = (_Kd0 * exp(-_Ed / (kb * (*_temperature)[k]))) ;
}

void Dissolution :: computeDrivForce()
{
    // Dissolution driving force cannot be negative | hydride will decrease from Cpp to 0 (Ctot -> Css); solution will increased from Css to TSSd
    for(int k=0; k<_nbCells; k++)
        _driving_force[k] = max(min((*_tssd)[k]-(*_Css)[k], (*_Ctot)[k]-(*_Css)[k]), 0.) ;
}
