/**
    This class manages the input reading and output writing for the HNGD model
 */
#ifndef InOut_hpp
#define InOut_hpp

#include <stdio.h>

#include "HNGD.hpp"

using namespace std;

class InOut
{
    public:
    // Read the setting file (and then write Settings In Check)
        static void getSettings(int nb, double* settings, string path_exec, string file_name, string check_name, vector<double>& customDist);
    
    // Read the physical parameters file (and then write PhysicsPara. In Check)
        static void getPhysics(int nb, double* physics, string path_exec, string file_name, string check_name);
    
    // Read the temperature history (and then write ThermalTrea. information In Check)
        static vector<vector<double>> getThermalTreatment(string path_exec, string file_name, string check_name);
    
    // Read the initial hydrogen profile
        static vector<vector<double>> getICHydrogen(string path_exec, string file_name);
    
    // Write the input information in a file for verification purpose
        static void writeSettingsInCheck(double * settings, string path_exec, string check_name);
        static void writePhysicsInCheck(double * physicalParameters, string path_exec, string check_name);

    // Write the output
        static void writeInitialOutput(HNGD hngd, string path_exec, string output_name, int nbNodes, int nbOutput, int nbPosPrint, int* listPosPrint, int geometry, vector<double> customDist) ; 
        static void writeOuput(HNGD hngd, string path_exec, string output_name, int nbNodes, int nbOutput, double t, double temp, int nbPosPrint, int* listPosPrint);            
};

#endif