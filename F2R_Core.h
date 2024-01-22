#ifndef F2R_CORE
#define F2R_CORE

#include<string>
#include"TROOT.h"

using label_type = UShort_t;
using time_type = ULong64_t;
using nrj_type = Int_t;

//Create, fill, and save the converted root file for a given FASTER file
void Convert(std::string filename);

#endif //F2R_CORE