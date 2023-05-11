#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>

#include "InOut.hpp"
#include "HNGD.hpp"
#include "EvalEvolution.hpp"

using namespace std ;

int main(int argc, char* argv[]) {
    // The following function and variables are used to determine 
    // if the current state should be written in the output                  
    vector<double> interpolate(double t, vector<double> time_stamps, vector<vector<double>> temp_stamps);
    // This function checks if a time stamp specified in the 
    // temperature input file was passed.
    bool changeInterval(double t, double dt, vector<double> time_stamps);
    // The EvalEvolution objects determine if a given quantity
    // changed sufficiently to justify a print in the output
    EvalEvolution evalEvolTemp ; // T profile
    EvalEvolution evalEvolHyd ;  // H concen. profile

    //----------------------- DEFINE EXECUTION FOLDER AND INPUT FILES NAMES --------------------
    // Path to the folder to use /*custom*/
    string path_exec = "C:\\Users\\i1392\\Desktop\\HNGD_VScode\\" ;
    // Name of the folder containing the input files
    // The input files specified by the argument given at
    // launch.json must be placed in this folder
    string input_folder = "input_files\\" ;
    // Name of the folder cantaining the result files
    string result_folder = "result_files\\" ;
    // check command line
    cout << "path_exec = " << path_exec << endl ;
    cout << "argc = " << argc << endl ;
    cout << "argv = " ;
    for (int i = 0; i < argc; i++) {
      cout << argv[i] << " " ;
    } cout << endl;
  
    //------------------------ CHECK THE GIVEN ARGUMENTS/CASES ---------------------------------
    vector<string> cases ;
    // Default case name
    string def_name = "ex_lin" ;
    // If user didn't pass any argument to launch.json, 
    // case will be the default case
    if(argc == 1) {  
      cout << "default case = " << def_name;
      cases.push_back(def_name);
      cout << endl ;
    //argv[1] argv[2], argv[3], argv[4] ... argv[argc-1]
    }else {
      cout << "cases = " ;
      for (int i = 1; i < argc; i++) {
        cout << argv[i] << ", ";
        cases.push_back(argv[i]);
      } cout << endl ;
    } 
    cout << endl ;

    //---------------------- EXECUATE EACH CASE ------------------------------------------------
    for (int i = 0; i < cases.size(); i++) 
    { 
      // Name of the simulation case
      string name = cases[i];
      cout << "------------------------------------------------------- ";
      cout << "THIS IS CASE " << name;
      cout << " -------------------------------------------------------" << endl;
      // input_files  
      string settings_name  = input_folder + name + "_set.txt" ;
      string treatment_name = input_folder + name + "_temp.txt";
      string physics_name   = input_folder + name + "_phys.txt";
      string hydroIC_name   = input_folder + name + "_hyd.txt" ;
      // result_files
      string check_name     = result_folder + name + "_input_check.txt";
      string output_name    = result_folder + name + "_out.csv" ;

      //-------------------- INPUT READING ------------------------
      // Simulation settings contained in the *_set.txt file
      short int nbSettings = 7 ;
      double settings[nbSettings];
      vector<double> customDist ;
      InOut::getSettings(nbSettings, settings, path_exec, settings_name, check_name, customDist); 
      // Physical parameters contained in the *_phys.txt file
      short int nbPhysicalParameters = 19 ;
      double physicalParameters[nbPhysicalParameters];
      InOut::getPhysics(nbPhysicalParameters, physicalParameters, path_exec, physics_name, check_name);
      // Thermal treatment contained in the *_temp.txt file
      vector<double> time_temp(0);      // time stamps for temperature (s)
      vector<double> pos_temp(0);       // positions for temperature (cm)
      vector<vector<double>> temp_inp;  // temperature values (K)
      vector<vector<double>> thermal_treatment = InOut::getThermalTreatment(path_exec, treatment_name, check_name);
      pos_temp = thermal_treatment[0] ;
      time_temp= thermal_treatment[1] ;
      double t_end = time_temp[time_temp.size()-1] ;
      for(int k=0; k<thermal_treatment.size()-2; k++)
        temp_inp.push_back(thermal_treatment[k+2]) ;
      // Initial H profile contained in the *_hyd.txt file
      vector<vector<double>> hydrogenIC = InOut::getICHydrogen(path_exec, hydroIC_name);
      vector<double> pos_hyd = hydrogenIC[0] ;
      vector<double> hyd_inp = hydrogenIC[1] ;
      
      //---------------- SYSTEM INITIALIZATION ---------------------
      // The HNGD object collects the input information to build
      // a Sample and the objects associated with each phenomenon
      HNGD hngd(settings, physicalParameters) ;
      // Some of the settings are needed for the time loop
      int nbNodes        = settings[0];
      // int nbPosPrint     = settings[0];
      double dtPrint     = settings[4];
      double critPrint   = settings[5] ;

      // Initialize the temperature and hydrogen profiles
      double t = 0. ;
      vector<double> temp = interpolate(t, time_temp, temp_inp);
      vector<double> mem_temp = temp ;
      hngd.getInitialConditions(pos_hyd, hyd_inp, pos_temp, temp);
      
      // Associate the EvalEvolution objects to the profiles
      evalEvolHyd.setProfile(hngd.returnSample()->returnTotalContent()) ;
      evalEvolHyd.setCriterion(critPrint) ;
      evalEvolTemp.setProfile(hngd.returnSample()->returnTemperature()) ;
      evalEvolTemp.setCriterion(critPrint) ;

      // Initialize the output file /*custom*/
      ofstream output;
      output.open(path_exec + output_name, ios::out);
      if (output.fail()) {
        cout << "File opening error!\nProgram stopped.\n";
        exit(1);
      }
      const short int nbOutput = 5 ; 
      int nbPosPrint  = customDist.size() ;
      if(nbPosPrint == 0) nbPosPrint = nbNodes; 
      // If user didn't input customDist, then set ndPosPrint equal to nbNodes
      int listPosPrint[nbPosPrint] ;
      InOut::writeInitialOutput(hngd, path_exec, output_name, nbNodes, nbOutput, nbPosPrint, listPosPrint, settings[6], customDist); //TODO: less parameters
      InOut::writeOuput(hngd, path_exec, output_name, nbNodes, nbOutput, t, 0., nbPosPrint, listPosPrint);
      

      //-------------------- TIME LOOP --------------------------
      double printCountdown(0.);
      InOut::writeOuput(hngd, path_exec, output_name, nbNodes, nbOutput, t, 0., nbPosPrint, listPosPrint);
      
      // exit(21817) ; 

      do
      {
        // Interpolation of input temperature using the function "interpolate" implemented below
        t += hngd.returnTimeStep() ;
        printCountdown += hngd.returnTimeStep() ;
        temp = interpolate(t, time_temp, temp_inp);

        // Check if the temperature profile changed
        bool T_changed = false ;
        for(int k=0; k<temp.size(); k++)
        if(abs(temp[k] - mem_temp[k]) / mem_temp[k] > 0.01)
        {
            T_changed = true ;
            mem_temp = temp ;
            break ;
        }
            
        // Compute the new system state
        hngd.getInput(T_changed, pos_temp, temp);
        hngd.compute();

        // Write output if enough time has elapsed or if the temperature profile has changed
        // or if the hydrogen profile has changed of if a time stamp was reached
        if (printCountdown >= dtPrint ||
            evalEvolTemp.evaluate(hngd.returnSample()->returnTemperature()) ||
            evalEvolHyd.evaluate(hngd.returnSample()->returnTotalContent()) ||
            changeInterval(t, hngd.returnTimeStep(), time_temp))
        {
            InOut::writeOuput(hngd, path_exec, output_name, nbNodes, nbOutput, t, 0., nbPosPrint, listPosPrint);
            printCountdown = 0. ;
        }
        
      } while ( t < t_end );
    }
    // ------------------- END OF COMPUTATION -----------------------
    cout  << "The calculation was performed!\n";
    // Sound notification at the end of simulation /*custom*/
    system("\"C:/Program Files (x86)/VideoLAN/VLC/vlc.exe\"  C:\\Users\\i1392\\Desktop\\HNGD_VScode\\notif.mp3");
    return 0;
  }


vector<double> interpolate(double t, vector<double> time_stamps, vector<vector<double>> temp_stamps)
{
  // If end of simulation, no interpolation
  if(t >= time_stamps[time_stamps.size()-1])
    return temp_stamps[temp_stamps.size()-1] ;
  else
  {
    int k=0 ;
    while(t >= time_stamps[k])
      k ++ ;
    vector<double> interpolated_temp_stamps(temp_stamps[0].size()) ;
    for(int i=0; i<temp_stamps[0].size(); i++)
      interpolated_temp_stamps[i] = temp_stamps[k-1][i] + (temp_stamps[k][i] - temp_stamps[k-1][i]) * (t - time_stamps[k-1]) / (time_stamps[k] - time_stamps[k-1]) ;
    return interpolated_temp_stamps ;
  }
}


bool changeInterval(double t, double dt, vector<double> time_stamps)
{
  // If end of simulation, no interpolation
  if(t >= time_stamps[time_stamps.size()-1] || t+dt >= time_stamps[time_stamps.size()-1])
    return true ;
  int k=0 ;
  while(t >= time_stamps[k])
    k ++ ;
  return t+dt > time_stamps[k] ;
}