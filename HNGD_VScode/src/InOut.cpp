#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <string>

#include "InOut.hpp"

using namespace std;

// -------------------------------- Input reading --------------------------------
void InOut::getSettings(int nb, double* settings, string path_exec, string file_name, string check_name, vector<double>& customDist)
{
    ifstream inpset ;
    inpset.open(path_exec + file_name);
    
    if(inpset.fail())
    {
        cout << "Setting file opening failed" << endl ;
        cout << path_exec + file_name;
        exit(1) ;
    }
    
    for(int k=0; k < nb; k++)
    {
      inpset >> settings[k] ;
      inpset.ignore(999,'\n');
      inpset.ignore(999,'\n');
    }
    // get customDist (if exists)
    double pos = 0.0 ;
    while(inpset >> pos) {
      customDist.push_back(pos) ;
    }
    
    InOut::writeSettingsInCheck(settings, path_exec, check_name);
    
    inpset.close();
}

void InOut::getPhysics(int nb, double* physics, string path_exec, string file_name, string check_name)
{
    ifstream inpphys ;
    inpphys.open(path_exec + file_name);
    
    if(inpphys.fail())
    {
        cout << "Physics file opening failed" << endl ;
        exit(1) ;
    }
    
    for(int k=0; k<nb; k++)
    {
      inpphys >> physics[k] ;
      inpphys.ignore(999,'\n');
    }
    
    InOut::writePhysicsInCheck(physics, path_exec, check_name);
    
    inpphys.close();
}

vector<vector<double>> InOut::getThermalTreatment(string path_exec, string file_name, string check_name)
{
    ifstream inptemp ;
    inptemp.open (path_exec + file_name);
    if(inptemp.fail())
    {
        cout << "Temperature file opening failed" << endl ;
        exit(1) ;
    }
    
    ofstream chkinp ;
    chkinp.open(path_exec + check_name, std::ios_base::app);
    chkinp << endl << endl << "time (s), temperature (K), gradient (K/cm)";
    
    vector<double> pos_temp(0) ;
    vector<double> time_temp(0) ;
    vector<vector<double>> temp_inp(0) ;
    
    // Positions
    string line_pos ;
    getline(inptemp, line_pos);
    istringstream stream_pos(line_pos) ;
    double pos ;
    while(stream_pos >> pos)
        pos_temp.push_back(pos) ;

    unsigned long nb_pos = pos_temp.size() ;

    unsigned short int i=0;
    do
    {
        double time,temperature ;
        inptemp >> time ;
        time_temp.push_back(time);

        if(i>0 && (time_temp[i] < time_temp[i-1])) // check that the input is consistent
        {
          cout  << "Input time not monotonically increasing!\nPlease check input.\n";
          exit(1);
        }

        vector<double> new_vec(0) ;
        for(int k=0; k<nb_pos; k++)
        {
          inptemp >> temperature ;
          new_vec.push_back(temperature);
        }
        temp_inp.push_back(new_vec) ;

        i++ ;

    } while(!(inptemp.eof()));

    if(time_temp[i-1] == time_temp[i-2])
    {
        temp_inp.pop_back() ;
        time_temp.pop_back();
    }

    cout << "Number of specified temperature steps = " << i-1 << endl;
    double t_end = time_temp[i-1];
    chkinp << endl << endl << "end time (s) = " << t_end;
    
    // Return the thermal treatment as a vector of vector<double>
    // First vector contains the positions
    // Second vector contains the time stamps
    // The others contain the temperature values
    vector<vector<double>> thermal_treatment(0) ;
    thermal_treatment.push_back(pos_temp) ;
    thermal_treatment.push_back(time_temp);
    for(int k=0; k<temp_inp.size(); k++)
        thermal_treatment.push_back(temp_inp[k]);

    cout << "pos_temp = " ;
    for (int i = 0; i < pos_temp.size(); i++) {
      cout << pos_temp[i] << ", ";
    } cout << endl ;
    cout << "time_temp = " ;
    for (int i = 0; i < time_temp.size(); i++) {
      cout << time_temp[i] << ", ";
    } cout << endl ;
    for (int i = 0; i < temp_inp.size(); i++) {
      cout << "temp_inp[" << i <<  "] (@ " << time_temp[i] << "s ) = ";
      for (int j = 0; j < temp_inp[i].size(); j++) {
        cout << temp_inp[i][j] << ", ";
      } cout << endl;
    } cout << endl ;

    return thermal_treatment ;
}

vector<vector<double>> InOut::getICHydrogen(string path_exec, string file_name)
{
    ifstream inphyd;
    inphyd.open(path_exec + file_name);
    if(inphyd.fail())
    {
        cout << "Hydrogen IC file opening failed" << endl ;
        exit(1) ;
    }
    
    vector<double> pos_hyd(0);
    vector<double> hyd_inp(0);
    
    string line_pos ;
    string line_hyd ;
    
    double pos ;
    double hyd ;
    
    getline(inphyd, line_pos);
    getline(inphyd, line_hyd);
    
    istringstream stream_pos_hyd(line_pos) ;
    istringstream stream_hyd(line_hyd) ;
    
    while(stream_pos_hyd >> pos){
      pos_hyd.push_back(pos) ;
      stream_hyd >> hyd ;
      hyd_inp.push_back(hyd) ;
    }
    
    vector<vector<double>> hydrogenIC(0) ;
    hydrogenIC.push_back(pos_hyd);
    hydrogenIC.push_back(hyd_inp);
    cout << "pos_hyd = " ;
    for (int i = 0; i < pos_hyd.size(); i++) {
      cout << pos_hyd[i] << ", ";
    } cout << endl ;
      cout << "pos_hyd = " ;
    for (int i = 0; i < hyd_inp.size(); i++) {
      cout << hyd_inp[i] << ", ";
    } cout << endl << endl ;  
    return hydrogenIC ;
}

// -------------------------------- Check file writing --------------------------------
void InOut::writeSettingsInCheck(double * settings, string path_exec, string check_name)
{
  ofstream chkinp ;
  chkinp.open(path_exec + check_name, ios::out);
  // nbNodes
  chkinp << "Space discretization: " << settings[0] << " cells\n";
  // bias
  chkinp << "Bias: " << settings[1] << endl ;
  // Lenght
  chkinp << "Sample lenght: " << settings[2] << endl;
  // dt
  if(settings[3]<=0)
    chkinp << "Adaptative time step" << endl;
  else
    chkinp << "Time step: " << settings[3] << endl;
  // dtPrint
  chkinp << "Printing time step: " << settings[4] << endl ;
  // Evolution criterion
    chkinp << "Evolution criterion: " << settings[5] << endl ;
  chkinp.close();
}

void InOut::writePhysicsInCheck(double * physicalParameters, string path_exec, string check_name)
{
  ofstream chkinp ;
  chkinp.open(path_exec + check_name, std::ios_base::app);
  chkinp << "Kd0:\t"  << physicalParameters[0]  << "\ts-1\n" ;
  chkinp << "Ediss:\t"<< physicalParameters[1]  << "\teV/at\n" ;
  chkinp << "Kn0:\t"  << physicalParameters[2]  << "\ts-1\n" ;
  chkinp << "Eth0:\t" << physicalParameters[3]  << "\teV/mol\n" ;
  chkinp << "Eth1:\t" << physicalParameters[4]  << "\teV/mol/K\n" ;
  chkinp << "Eth2:\t" << physicalParameters[5]  << "\teV/mol/K2\n" ;
  chkinp << "Eth3:\t" << physicalParameters[6]  << "\teV/mol/K3\n" ;
  chkinp << "Kmob0:\t"<< physicalParameters[7]  << "\ts-1\n" ;
  chkinp << "Kth0:\t" << physicalParameters[8] << "\ts-1\n" ;
  chkinp << "Eg:\t"   << physicalParameters[9] << "\teV/at\n" ;
  chkinp << "TSSp0:\t"<< physicalParameters[10] << "\twt.ppm\n" ;
  chkinp << "Qtssp:\t"<< physicalParameters[11] << "\tJ/(mol.K)\n" ;
  chkinp << "TSSd0:\t"<< physicalParameters[12] << "\twt.ppm\n" ;
  chkinp << "Qtssd:\t"<< physicalParameters[13] << "\tJ/(mol.K)\n" ;
  chkinp << "D0:\t"   << physicalParameters[14] << "\tJ/(mol.K)\n" ;
  chkinp << "Q*:\t"   << physicalParameters[15] << "\tJ/mol\n" ;
  chkinp.close();
}

// -------------------------------- Output file writing --------------------------------
void InOut::writeOuput(HNGD hngd, string path_exec, string output_name, int nbNodes, int nbOutput, double t, double temp, int nbPosPrint, int* listPosPrint)
{
    ofstream output ;
    output.open(path_exec + output_name, std::ios_base::app);
    /*custom*/
    std::vector<double> listVector[nbOutput];
    listVector[0] = hngd.returnSample()->returnTotalContent();
    listVector[1] = hngd.returnSample()->returnSolutionContent();
    listVector[2] = hngd.returnSample()->returnHydrideContent();
    listVector[3] = hngd.returnSample()->returnTSSd();
    listVector[4] = hngd.returnSample()->returnTSSp();
    output << hngd.returnTimeStep() << "," << t << ","  ;
    for(int i=0; i<nbOutput; i++){
        for(int j=0; j<nbPosPrint; j++)
        {
            output << (listVector[i])[listPosPrint[j]] << ",";
        }
        output << ',' ;
    }
    output << "\n"   ;
    cout << "t= " << t << "s"<< endl;
    output.close();
}

void InOut :: writeInitialOutput(HNGD hngd, string path_exec, string output_name, int nbNodes, int nbOutput, int nbPosPrint, int* listPosPrint, int geometry, vector<double> customDist)
{
  ofstream output ;
  output.open(path_exec + output_name, std::ios_base::app);
  if(nbNodes>0)
  {
    std::vector<double> positionVector = hngd.returnSample()->returnPosition();
    // --------------------------------------------------------------------------------------------------------------------------------
    // ------------------------------------- THIS CODE IS REPLACED BY THE CODE BELOW FOR CONVINENCE -----------------------------------
    // Defining the positions where the values will be printed
      // for(int i=0; i<nbNodes; i++)
      //   *(listPosPrint+i) = 0 ;
      // if(nbPosPrint==nbNodes){
      //   for(int i=0; i<nbPosPrint; i++)
      //     *(listPosPrint+i) = i ;
      // }
      // else{
      //   for(int i=1; i<nbPosPrint-1; i++){
      //     while(positionVector[listPosPrint[i]] <= i*positionVector[nbNodes-1]/nbPosPrint)
      //       *(listPosPrint+i) += 1 ;
      //   }
      //   *(listPosPrint+nbPosPrint-1) = nbNodes-1;
      // }
    // ------------------------------------- THIS CODE IS REPLACED BY THE CODE BELOW FOR CONVINENCE -----------------------------------
    // --------------------------------------------------------------------------------------------------------------------------------
    
    // Initialize listPosPrint
    for (int i = 0; i < nbPosPrint; i++)
      *(listPosPrint+i) = 0 ;

  // If user didn't input customDist, nbPosPrint will equal to nbNodes  
    if (nbPosPrint == nbNodes) {
      // For each i in listPosPrint, i will equal to i
      for(int i=0; i<nbPosPrint; i++)
        *(listPosPrint+i) = i ;
    } else {
  // If user inputted customDist, nbPosPrint will not equal to nbNodes
      for(int i = 1; i < nbPosPrint-1; i++) {
        // For each i in listPosPrint, choose the closest ndCellDist depend on customDist
        // except for 0 (first must equal to 0) and nbPosPrint-1 (last must equal to nbNodes-1)
        while(positionVector[listPosPrint[i]] <= customDist[i])
          *(listPosPrint+i) += 1 ;
      }
      *(listPosPrint+nbPosPrint-1) = nbNodes-1;
    }

    cout << "listPosPrint = " ;
    for (int i = 0; i < nbPosPrint; i++){
      cout << listPosPrint[i] << ", ";
    } cout  << endl ;
    cout << "customDist = ";
    if (customDist.size()==0) {
      cout << "null" ;
    }else {
      for (int i = 0; i < customDist.size(); i++) 
        cout << customDist[i] << ", " ;
    } cout << endl ;
    cout << "ndcellDist = " ;
    for (int i = 0; i < nbPosPrint; i++){
      cout << positionVector[listPosPrint[i]] << ", " ;
    } cout << endl << endl ;
    
    output << "dt [s],Time [s]," ;
    string listOutputNames[nbOutput];
    /*custom*/
    listOutputNames[0] = "Ctot," ;
    listOutputNames[1] = "Css,"  ;
    listOutputNames[2] = "Cprec,";
    listOutputNames[3] = "TSSd," ;
    listOutputNames[4] = "TSSp," ;
    for(int j=0; j<nbOutput; j++){
      for(int i=0; i<nbPosPrint; i++) // Print the type of output as column head
        output << listOutputNames[j] ;
      output << ',' ;
    }
    output << '\n' ;
    // Print the positions as column "sub-head"
    if (geometry > 0)
        output << "Positions, [rad],";
    else
        output << "Positions, [cm],";
    for(int j=0; j<nbOutput; j++){
      for(int i=0; i<nbPosPrint; i++){
        output << positionVector[listPosPrint[i]] << "," ;
      }
      output << ',' ;
    }
  }
  else
  {
    output << "Time(s),Temp(K),TSSp,TSSd,Css,Chyd,Ctot" << endl;
  }
  output.close();
}