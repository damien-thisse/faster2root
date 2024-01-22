#include "F2R_Parameters.h"

#include<string>
#include<sstream>
#include<fstream>
#include<iostream>

using namespace std;

F2R_Parameters::F2R_Parameters(){
    ifstream F2R_ParametersFile("parameters.dat");
    
    //Default parameters
    F2R_ListOfFiles = "filename.dat";
    F2R_OutputDirectory = "Output";
    F2R_IDTableFile = "ID.dat";
    F2R_NumberOfThreads = 2;
    F2R_KeepGroups = false;

    F2R_MaxEnergyBranchNeeded = 1;

    //Test if the file exists
    if(!F2R_ParametersFile.is_open()){
        cerr << "Warning: ID Table file not found ... : " << "parameters.dat" << endl;
        cerr << "Default values will be used." << endl;
    }
    else{
        while (F2R_ParametersFile.good()){
        string oneline;
        getline(F2R_ParametersFile, oneline);
        istringstream is(oneline);
        string temp1 = "NULL";
        is>>temp1;
        if(temp1 == "final_directory:"){
            is >> F2R_OutputDirectory;
        }
        else if(temp1 == "number_of_threads:"){
            is >> F2R_NumberOfThreads;
        }
        else if(temp1 == "input_file:"){
            is >> F2R_ListOfFiles;
        }
        else if(temp1 == "id_file:"){
            is >> F2R_IDTableFile;
        }
        else if(temp1 == "keep_groups:"){
            string bb;
            is >> bb;
            if(bb == "true"){
                F2R_KeepGroups = true;
            }
        }
        else if(temp1 != "NULL"){
            cout << "WARNING: Unknown parameter " << temp1 << " will be ignored !" << endl;
        }
        }
    }
    cout << endl << "-----------Parameters--------------" << endl;
    cout << "Output directory set to : " << F2R_OutputDirectory << endl;
    cout << "Number of threads set to : " << F2R_NumberOfThreads << endl;
    cout << "Input file set to : " << F2R_ListOfFiles << endl;
    cout << "ID Table file set to : " << F2R_IDTableFile << endl;
    cout << "Ungrouping of data set to : " << F2R_KeepGroups << endl << endl;



}

F2R_Parameters::~F2R_Parameters(){

}

bool F2R_Parameters::setIDTable(){
    ifstream F2R_IDTableData(F2R_IDTableFile);
    //Test if the file exists
    if(!F2R_IDTableData.is_open()){
        cerr << "Error: ID Table file not found ... : " << F2R_IDTableFile << endl;
        return false;
    }
    //Loop on all the data in the file
    while(F2R_IDTableData.good()){
      string oneline;
      getline(F2R_IDTableData, oneline);
      istringstream is(oneline);
      int temp1, temp2;
      is>>temp1>>temp2;
      if(oneline.size()>1){
        F2R_IDTable[temp1] = temp2;
        //Look for the declaration of QDC2, 3, or 4 to adapt the branches in the tree
        if(temp2 < 5 && temp2 > F2R_MaxEnergyBranchNeeded) F2R_MaxEnergyBranchNeeded = temp2;
      }
    }
    return true;
}