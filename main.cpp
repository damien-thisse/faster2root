#include "F2R_Parameters.h"
#include<iostream>

#include"TROOT.h"


using namespace std;

int main(int argc, char** argv){

    //Instance for the first time the parameters handler class
    F2R_Parameters& parameters = F2R_Parameters::getInstance();

    //If the ID file given to the parameter class doesn't exist, stop the code
    if(!parameters.setIDTable()){
        return -1;
    }



}