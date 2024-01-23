//  std includes
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <dirent.h>
#include <fstream>
#include <sstream>
#include <map>
#include <thread>
#include <mutex>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <cassert>


//  root includes
#include "TFile.h"
#include "TTree.h"
#include "TThread.h"
#include "TH2.h"
#include "TH1.h"
#include "TStopwatch.h"
#include "TMath.h"

//  fasterac includes
#include "fasterac/adc.h"
#include "fasterac/adc_caras.h"
#include "fasterac/electrometer.h"
#include "fasterac/farray.h"
#include "fasterac/fast_data.h"
#include "fasterac/fasterac.h"
#include "fasterac/group.h"
#include "fasterac/jdb_hv.h"
#include "fasterac/online.h"
#include "fasterac/plas.h"
#include "fasterac/qdc.h"
#include "fasterac/qdc_caras.h"
#include "fasterac/qt2t.h"
#include "fasterac/qtdc.h"
#include "fasterac/rf.h"
#include "fasterac/rf_caras.h"
#include "fasterac/sampler.h"
#include "fasterac/sampling.h"
#include "fasterac/scaler.h"
#include "fasterac/spectro.h"
#include "fasterac/utils.h"


using namespace std;

// Definition of the standard types
using label_type = UShort_t;
using time_type = ULong64_t;
using nrj_type = Int_t;
using pu_type = Bool_t;


// Defining function templates
void init_parameters(int & nOfThreads, string & output_directory, string & input_file, string & temp_directory, bool & kGroups); // charge les paramètres depuis parameters.dat
void convert(string FileName, map<int,int> myID, const string & temp_directory, const bool kGroups); // lance la fonction de conversion
int sorterF(string FileName, const string & temp_directory, const string & output_directory, map<int,int> myID); // refait un tri en temps de l'arbre final

//******************************************************************************//
//									                                                           	//
//	 				                            MAIN	                                  //
//										                                                          //
//******************************************************************************//

int main(int argc, char** argv)
{
  // Loading the ID mapping to recognize the types of detectors
  ifstream ID_file("/local/home/dt269135/Codes/faster2root/ID.dat");
  map<int,int> myID;

  if(ID_file.is_open())
  {
    while (ID_file.good())
    {
      string oneline;
      getline(ID_file, oneline);
      istringstream is(oneline);
      int temp1, temp2;
      is>>temp1>>temp2;
      if (oneline.size()>1)
      {
        myID[temp1] = temp2;
      }
    }
  }
  else
  {
    cerr << "No ID file loaded ... " << endl;
    return -1;
  }
  ID_file.close();

  //Default values for the parameters
  int nOfThreads = 1;
  string output_directory = "~/";
  string temp_directory = "~/";
  string input_file_name = "/local/home/dt269135/Codes/faster2root/filename.dat";
  bool kGroups = false;
  //Parameters are given from an external file
  init_parameters(nOfThreads, output_directory, input_file_name, temp_directory, kGroups);

  ifstream input_file(input_file_name.c_str());
  Int_t filecounter = 0;
  if(input_file.is_open())
  {
    while(input_file.good())
    {
      string oneline;
      getline(input_file, oneline);
      if(oneline.size() > 1)
      {
        filecounter++;
      }
    }
    input_file.close();
  }
  else
  {
    cerr << "Error opening input file ..." << endl;
    return -1;
  }

  input_file.close();
  input_file.open(input_file_name.c_str());

  TThread::Initialize();
  vector<thread> threads;
  if(nOfThreads > std::thread::hardware_concurrency())
  {
    nOfThreads = std::thread::hardware_concurrency();
    cout << "Number of threads too large (hardware) -> reset to " << nOfThreads << endl;
  }
  if(nOfThreads > filecounter)
  {
    nOfThreads = filecounter;
    cout << "Number of threads too large (too few files to be processed) -> reset to " << nOfThreads << endl;
  }

  bool endOfFile = false;
  mutex myInputMutex;

  //Start of the conversion in multi-thread mode
  //Each time a file is converted, a new one is given to the thread
  for(int k = 0; k < nOfThreads; k++)
  {
    threads.emplace_back([&input_file, myID, &endOfFile, &myInputMutex, output_directory,temp_directory, kGroups](){

      while(!endOfFile)
      {
        string filename = "";
        myInputMutex.lock();
        input_file >> filename;
        if(filename == "")
        {
          endOfFile = true;
          myInputMutex.unlock();
          break;
        }
	cout << "Thread #" << std::this_thread::get_id() << ": start processing " << filename << endl;
        myInputMutex.unlock();

        convert(filename, myID, output_directory, kGroups);
        if(!kGroups)
        {
          sorterF(filename,output_directory, output_directory, myID);
        }
      }
    });
  }
  for(unsigned int i = 0; i< threads.size(); i++)
  {
    threads.at(i).join();
  }

  return 0;
}

void convert(string FileName, map<int,int> myID, const string & temp_directory, const bool kGroups)
{
  //Boolean variables used if there are QDC with several integration gates
  bool energy2 = false;
  bool energy3 = false;
  bool energy4 = false;

  for(const auto & elem : myID)
  {
    if(elem.second == 2) {energy2 = true;}
    else if(elem.second == 3) {energy2 = true;energy3 = true;}
    else if(elem.second == 4) {energy2 = true;energy3 = true; energy4 = true;}
  }

  string DATAFILENAME;
  string ROOTFILENAME;
  if (FileName.size()>1 )
  {
    //Automatically change the path and the extension of the filename for the output file
    int last_slash =  FileName.find_last_of('/');
    DATAFILENAME = FileName;
    FileName.erase(FileName.end()-5, FileName.end());
    FileName.erase(FileName.begin(), FileName.begin()+last_slash+1);
    FileName.insert(FileName.end()-4,'R');
    FileName.insert(FileName.end()-4,'_');
    ROOTFILENAME = temp_directory+FileName+".root";

    /************/
    /*  FASTER  */
    /************/

    //  file reader
    faster_file_reader_p   reader;

    //  data
    faster_data_p data;
    unsigned char alias;
    unsigned short label;
    ULong64_t clock_g;
    ULong64_t clock_a;

    //  group data
    faster_buffer_reader_p group_reader;
    char                   group_buffer [15000];
    faster_buffer_reader_p inside_group_reader;
    char                   inside_group_buffer [15000];
    unsigned short         lsize;
    unsigned short         inside_lsize;
    faster_data_p          group_data;
    faster_data_p          inside_group_data;

    //Declaration of the different types of FASTER structures
    //RF
    rf_data rf;

    //Tref
    tref_tdc TREF_tdc;

    //QDC 1 gate
    qdc_t_x1 QDC1_qdc;

    //QDC 2 gates
    qdc_t_x2 QDC2_qdc;

    //QDC 3 gates
    qdc_t_x3 QDC3_qdc;

    //QDC 4 gates
    qdc_t_x4 QDC4_qdc;

    //Trapezoidal
    trapez_spectro trapez_adc;

    //CRRC4
    crrc4_spectro crrc4_adc;

    /**********/
    /*  ROOT  */
    /**********/
    /*************/
    /* DO I keep FASTER events */
    Bool_t Keep_groups = kGroups;

    //Branches of the output ROOT tree
    std::vector<nrj_type> leaf_nrj; //energy
    std::vector<nrj_type> leaf_nrj2; //used if QDC2
    std::vector<nrj_type> leaf_nrj3; //used if QDC3
    std::vector<nrj_type> leaf_nrj4; //used if QDC4
    std::vector<time_type> leaf_t; //time
    std::vector<label_type> leaf_label; //label
    std::vector<pu_type> leaf_pu; //pile-up detection
    UShort_t leaf_mult;
    time_type leaf_group_time;
    nrj_type leaf_nrj_2; //energy
    nrj_type leaf_nrj2_2; //used if QDC2
    nrj_type leaf_nrj3_2; //used if QDC3
    nrj_type leaf_nrj4_2; //used if QDC4
    time_type leaf_t_2; //time
    label_type leaf_label_2; //label
    pu_type leaf_pu_2; //pile-up detection


    //  root tree
    TTree* tree;
    char tree_title[256];
    //  root file
    TString fName = ROOTFILENAME;
    TFile* root_file;
    //  print infos
    //printf ("\n");
    //printf ("  Transfering :\n");
    //printf ("     - read faster file '%s'\n", DATAFILENAME.c_str());
    //printf ("     - output those data to root file '%s'\n", ROOTFILENAME.c_str());
    //printf ("\n");
    //  open faster file reader
    reader = faster_file_reader_open (DATAFILENAME.c_str());
    if (reader == NULL)
    {
      printf ("error opening file %s\n", DATAFILENAME.c_str());
      return;
    }

    //Initialization of the output
    root_file = new TFile(fName.Data(), "recreate");
    sprintf (tree_title, "faster to tree root test : ungroup2tree.C");
    tree = new TTree ("DataTree", tree_title);
    // creating a tree structure where I don't break the grouping
    if(Keep_groups)
    {
      tree->Branch ("label", &leaf_label);
      tree->Branch ("nrj", &leaf_nrj);
      tree->Branch ("time", &leaf_t);
      tree->Branch ("pileup",&leaf_pu);
      if(energy2)
      {
        tree->Branch ("nrj2",&leaf_nrj2);
      }
      if(energy3)
      {
        tree->Branch ("nrj3",&leaf_nrj3);
      }
      if(energy4)
      {
        tree->Branch ("nrj4",&leaf_nrj4);
      }
      tree->Branch("multiplicity",&leaf_mult);
      tree->Branch("GroupTimeStamp",&leaf_group_time);
    }
    else
    {
      tree->Branch ("label", &leaf_label_2);
      tree->Branch ("nrj", &leaf_nrj_2);
      tree->Branch ("time", &leaf_t_2);
      tree->Branch ("pileup",&leaf_pu_2);
      if(energy2)
      {
        tree->Branch ("nrj2",&leaf_nrj2_2);
      }
      if(energy3)
      {
        tree->Branch ("nrj3",&leaf_nrj3_2);
      }
      if(energy4)
      {
        tree->Branch ("nrj4",&leaf_nrj4_2);
      }
    }
    //cout << "Starting convertion" << endl;
    //While there are date to be read...
    while ((data = faster_file_reader_next (reader)) != NULL)
    {
      alias = faster_data_type_alias (data);

      // If data are grouped
      if (alias == GROUP_TYPE_ALIAS) //in this condition, the data are grouped
      {

        leaf_nrj.clear(); //energy
        if(energy2) leaf_nrj2.clear(); //used if QDC2
        if(energy3) leaf_nrj3.clear(); //used if QDC3
        if(energy4) leaf_nrj4.clear(); //used if QDC4
        leaf_t.clear(); //time
        leaf_label.clear(); //label
        leaf_pu.clear(); //pile-up detection
        leaf_mult = 0;
        lsize = faster_data_load (data, group_buffer); // get group data
        group_reader = faster_buffer_reader_open (group_buffer, lsize); //open group reader
        while ((group_data = faster_buffer_reader_next (group_reader)) != NULL) // we loop on all the data inside the group
        {
          clock_g = faster_data_hr_clock_ns(group_data)*1000; //Time in ps

          //cout << " I have a groupe " << faster_data_label(group_data) << endl;
           leaf_group_time=clock_g;
          // I must include extra groups within groups
          if(faster_data_label(group_data)>2999)
          {
            inside_lsize = faster_data_load (group_data, inside_group_buffer);
            inside_group_reader = faster_buffer_reader_open (inside_group_buffer, inside_lsize);

            while((inside_group_data = faster_buffer_reader_next (inside_group_reader)) != NULL)
            {
              clock_g = faster_data_hr_clock_ns(inside_group_data)*1000; //Time in ps

              //RF
              if ( myID[faster_data_label(inside_group_data)] == 7)
              {
                faster_data_load (inside_group_data, &rf);
                // Getting absolute time for the hit
                clock_g = faster_data_hr_clock_ns(inside_group_data)*1000;  //Time in ps
                //if (LABR_qdc.q1_saturated!=1) // reject saturation
                {
                  leaf_label.push_back(faster_data_label (inside_group_data));
                  leaf_nrj.push_back(rf_period_ns(rf)*1000);
                  if(energy2) leaf_nrj2.push_back(0);
                  if(energy3) leaf_nrj3.push_back(0);
                  if(energy4) leaf_nrj4.push_back(0);
                  leaf_t.push_back((ULong64_t)clock_g);
                  leaf_pu.push_back(false);
                  leaf_mult++;
                }
              }

              //TRef
              if(myID[faster_data_label(inside_group_data)] == 8)
              {
                faster_data_load (inside_group_data, &TREF_tdc);
                //clock_g = 1000.*(faster_data_clock_ns (group_data) + tref_conv_dt_ns (TREF_tdc.tdc));
                clock_g = faster_data_hr_clock_ns(inside_group_data)*1000;
                leaf_pu.push_back(false);
                leaf_label.push_back(faster_data_label (inside_group_data));
                leaf_nrj.push_back(0);
                if(energy2) leaf_nrj2.push_back(0);
                if(energy3) leaf_nrj3.push_back(0);
                if(energy4) leaf_nrj4.push_back(0);
                leaf_t.push_back((ULong64_t)clock_g);
                leaf_mult++;
              }

              //trapezoidal
              else if ( myID[faster_data_label(inside_group_data)] == 5)
              {
                faster_data_load (inside_group_data, &trapez_adc);
                // Getting absolute time for the hit

                //clock_g = 1000.*(faster_data_clock_ns (group_data) + tref_conv_dt_ns (trapez_adc.tdc));
                //clock_g = faster_data_hr_clock_ns(group_data)*1000;
                clock_g = faster_data_hr_clock_ns(inside_group_data)*1000;
                if (trapez_adc.pileup!=1) leaf_pu.push_back(false);
                else leaf_pu.push_back(true);// reject pile up

                leaf_label.push_back(faster_data_label (inside_group_data));
                leaf_nrj.push_back(trapez_adc.measure);
                if(energy2) leaf_nrj2.push_back(0);
                if(energy3) leaf_nrj3.push_back(0);
                if(energy4) leaf_nrj4.push_back(0);
                leaf_t.push_back((ULong64_t)clock_g);
                leaf_mult++;
              }

              //QDC1
              else if ( myID[faster_data_label(inside_group_data)] == 1)
              {
                faster_data_load (inside_group_data, &QDC1_qdc);

                //clock_g = 1000.*(faster_data_clock_ns (group_data) + tref_conv_dt_ns (QDC1_qdc.tdc));
                clock_g = faster_data_hr_clock_ns(inside_group_data)*1000;
                if (QDC1_qdc.q1_saturated!=1) leaf_pu.push_back(false); // reject saturation
                else leaf_pu.push_back(true);

                leaf_label.push_back(faster_data_label (inside_group_data));
                leaf_nrj.push_back(QDC1_qdc.q1);
                if(energy2) leaf_nrj2.push_back(0);
                if(energy3) leaf_nrj3.push_back(0);
                if(energy4) leaf_nrj4.push_back(0);
                leaf_t.push_back((ULong64_t)clock_g);
                leaf_mult++;
              }

              //QDC2
              else if ( myID[faster_data_label (inside_group_data)] == 2)
              {
                faster_data_load (inside_group_data, &QDC2_qdc);
                // Getting absolute time for the hit
                //clock_g = (faster_data_clock_ns (group_data) + tref_conv_dt_ns (QDC2_qdc.tdc))*1000.;
                clock_g = faster_data_hr_clock_ns(inside_group_data)*1000;
                if (QDC2_qdc.q1_saturated!=1) leaf_pu.push_back(false);
                else leaf_pu.push_back(true); // reject saturation

                leaf_label.push_back(faster_data_label (inside_group_data));
                leaf_nrj.push_back(QDC2_qdc.q1);
                leaf_nrj2.push_back(QDC2_qdc.q2);
                if(energy3) leaf_nrj3.push_back(0);
                if(energy4) leaf_nrj4.push_back(0);
                leaf_t.push_back((ULong64_t)clock_g);
                leaf_mult++;
              }
              //QDC3
              else if ( myID[faster_data_label (inside_group_data)] == 3)
              {
                faster_data_load (inside_group_data, &QDC3_qdc);
                // Getting absolute time for the hit
                //clock_g = (faster_data_clock_ns (group_data) + tref_conv_dt_ns (QDC3_qdc.tdc))*1000.;
                clock_g = faster_data_hr_clock_ns(inside_group_data)*1000;
                if (QDC3_qdc.q1_saturated!=1) leaf_pu.push_back(false);
                else leaf_pu.push_back(true); // reject saturation

                leaf_label.push_back(faster_data_label (inside_group_data));
                leaf_nrj.push_back(QDC3_qdc.q1);
                leaf_nrj2.push_back(QDC3_qdc.q2);
                leaf_nrj3.push_back(QDC3_qdc.q3);
                if(energy4) leaf_nrj4.push_back(0);
                leaf_t.push_back((ULong64_t)clock_g);
                leaf_mult++;
              }

              //QDC4
              else if ( myID[faster_data_label (inside_group_data)] == 4)
              {
                faster_data_load (inside_group_data, &QDC4_qdc);
                // Getting absolute time for the hit
                //clock_g = (faster_data_clock_ns (group_data) + tref_conv_dt_ns (QDC4_qdc.tdc))*1000.;
                clock_g = faster_data_hr_clock_ns(inside_group_data)*1000;
                if (QDC4_qdc.q1_saturated!=1) leaf_pu.push_back(false);
                else leaf_pu.push_back(true); // reject saturation

                leaf_label.push_back(faster_data_label (inside_group_data));
                leaf_nrj.push_back(QDC4_qdc.q1);
                leaf_nrj2.push_back(QDC4_qdc.q2);
                leaf_nrj3.push_back(QDC4_qdc.q3);
                leaf_nrj4.push_back(QDC4_qdc.q4);
                leaf_t.push_back((ULong64_t)clock_g);
                leaf_mult++;

              }

              //CRRC4
              else if ( myID[faster_data_label (inside_group_data)] == 6)
              {
                faster_data_load (inside_group_data, &crrc4_adc);
                // Getting absolute time for the hit
                //clock_g = 1000.*(faster_data_clock_ns (group_data) + tref_conv_dt_ns (crrc4_adc.delta_t));
                clock_g = faster_data_hr_clock_ns(inside_group_data)*1000;
                leaf_pu.push_back(false);
                leaf_label.push_back(faster_data_label (inside_group_data));
                leaf_nrj.push_back(crrc4_adc.measure);
                if(energy2) leaf_nrj2.push_back(0);
                if(energy3) leaf_nrj3.push_back(0);
                if(energy4) leaf_nrj4.push_back(0);
                leaf_t.push_back((ULong64_t)clock_g);
                leaf_mult++;
              }
            }
            faster_buffer_reader_close(inside_group_reader);
          }

          //RF
          if ( myID[faster_data_label(group_data)] == 7)
          {
            faster_data_load (group_data, &rf);
            // Getting absolute time for the hit
            clock_g = faster_data_hr_clock_ns(group_data)*1000;  //Time in ps
            //if (LABR_qdc.q1_saturated!=1) // reject saturation
            {
              leaf_label.push_back(faster_data_label (group_data));
              leaf_nrj.push_back(rf_period_ns(rf)*1000);
              if(energy2) leaf_nrj2.push_back(0);
              if(energy3) leaf_nrj3.push_back(0);
              if(energy4) leaf_nrj4.push_back(0);
              leaf_t.push_back((ULong64_t)clock_g);
              leaf_pu.push_back(false);
              leaf_mult++;
            }
          }

          //TRef
          if(myID[faster_data_label(group_data)] == 8)
          {
            faster_data_load (group_data, &TREF_tdc);
            //clock_g = 1000.*(faster_data_clock_ns (group_data) + tref_conv_dt_ns (TREF_tdc.tdc));
            clock_g = faster_data_hr_clock_ns(group_data)*1000;
            leaf_pu.push_back(false);
            leaf_label.push_back(faster_data_label (group_data));
            leaf_nrj.push_back(0);
            if(energy2) leaf_nrj2.push_back(0);
            if(energy3) leaf_nrj3.push_back(0);
            if(energy4) leaf_nrj4.push_back(0);
            leaf_t.push_back((ULong64_t)clock_g);
            leaf_mult++;
          }

          //trapezoidal
          else if ( myID[faster_data_label(group_data)] == 5)
          {
            faster_data_load (group_data, &trapez_adc);
            // Getting absolute time for the hit

            //clock_g = 1000.*(faster_data_clock_ns (group_data) + tref_conv_dt_ns (trapez_adc.tdc));
            //clock_g = faster_data_hr_clock_ns(group_data)*1000;
            clock_g = faster_data_hr_clock_ns(group_data)*1000;
            if (trapez_adc.pileup!=1) leaf_pu.push_back(false);
            else leaf_pu.push_back(true);// reject pile up

            leaf_label.push_back(faster_data_label (group_data));
            leaf_nrj.push_back(trapez_adc.measure);
            if(energy2) leaf_nrj2.push_back(0);
            if(energy3) leaf_nrj3.push_back(0);
            if(energy4) leaf_nrj4.push_back(0);
            leaf_t.push_back((ULong64_t)clock_g);
            leaf_mult++;
          }

          //QDC1
          else if ( myID[faster_data_label(group_data)] == 1)
          {
            faster_data_load (group_data, &QDC1_qdc);

            //clock_g = 1000.*(faster_data_clock_ns (group_data) + tref_conv_dt_ns (QDC1_qdc.tdc));
            clock_g = faster_data_hr_clock_ns(group_data)*1000;
            if (QDC1_qdc.q1_saturated!=1) leaf_pu.push_back(false); // reject saturation
            else leaf_pu.push_back(true);

            leaf_label.push_back(faster_data_label (group_data));
            leaf_nrj.push_back(QDC1_qdc.q1);
            if(energy2) leaf_nrj2.push_back(0);
            if(energy3) leaf_nrj3.push_back(0);
            if(energy4) leaf_nrj4.push_back(0);
            leaf_t.push_back((ULong64_t)clock_g);
            leaf_mult++;
          }

          //QDC2
          else if ( myID[faster_data_label (group_data)] == 2)
          {
            faster_data_load (group_data, &QDC2_qdc);
            // Getting absolute time for the hit
            //clock_g = (faster_data_clock_ns (group_data) + tref_conv_dt_ns (QDC2_qdc.tdc))*1000.;
            clock_g = faster_data_hr_clock_ns(group_data)*1000;
            if (QDC2_qdc.q1_saturated!=1) leaf_pu.push_back(false);
            else leaf_pu.push_back(true); // reject saturation

            leaf_label.push_back(faster_data_label (group_data));
            leaf_nrj.push_back(QDC2_qdc.q1);
            leaf_nrj2.push_back(QDC2_qdc.q2);
            if(energy3) leaf_nrj3.push_back(0);
            if(energy4) leaf_nrj4.push_back(0);
            leaf_t.push_back((ULong64_t)clock_g);
            leaf_mult++;
          }
          //QDC3
          else if ( myID[faster_data_label (group_data)] == 3)
          {
            faster_data_load (group_data, &QDC3_qdc);
            // Getting absolute time for the hit
            //clock_g = (faster_data_clock_ns (group_data) + tref_conv_dt_ns (QDC3_qdc.tdc))*1000.;
            clock_g = faster_data_hr_clock_ns(group_data)*1000;
            if (QDC3_qdc.q1_saturated!=1) leaf_pu.push_back(false);
            else leaf_pu.push_back(true); // reject saturation

            leaf_label.push_back(faster_data_label (group_data));
            leaf_nrj.push_back(QDC3_qdc.q1);
            leaf_nrj2.push_back(QDC3_qdc.q2);
            leaf_nrj3.push_back(QDC3_qdc.q3);
            if(energy4) leaf_nrj4.push_back(0);
            leaf_t.push_back((ULong64_t)clock_g);
            leaf_mult++;
          }

          //QDC4
          else if ( myID[faster_data_label (group_data)] == 4)
          {
            faster_data_load (group_data, &QDC4_qdc);
            // Getting absolute time for the hit
            //clock_g = (faster_data_clock_ns (group_data) + tref_conv_dt_ns (QDC4_qdc.tdc))*1000.;
            clock_g = faster_data_hr_clock_ns(group_data)*1000;
            if (QDC4_qdc.q1_saturated!=1) leaf_pu.push_back(false);
            else leaf_pu.push_back(true); // reject saturation

            leaf_label.push_back(faster_data_label (group_data));
            leaf_nrj.push_back(QDC4_qdc.q1);
            leaf_nrj2.push_back(QDC4_qdc.q2);
            leaf_nrj3.push_back(QDC4_qdc.q3);
            leaf_nrj4.push_back(QDC4_qdc.q4);
            leaf_t.push_back((ULong64_t)clock_g);
            leaf_mult++;

          }

          //CRRC4
          else if ( myID[faster_data_label (group_data)] == 6)
          {
            faster_data_load (group_data, &crrc4_adc);
            // Getting absolute time for the hit
            //clock_g = 1000.*(faster_data_clock_ns (group_data) + tref_conv_dt_ns (crrc4_adc.delta_t));
            clock_g = faster_data_hr_clock_ns(group_data)*1000;
            leaf_pu.push_back(false);
            leaf_label.push_back(faster_data_label (group_data));
            leaf_nrj.push_back(crrc4_adc.measure);
            if(energy2) leaf_nrj2.push_back(0);
            if(energy3) leaf_nrj3.push_back(0);
            if(energy4) leaf_nrj4.push_back(0);
            leaf_t.push_back((ULong64_t)clock_g);
            leaf_mult++;
          }
        }

        faster_buffer_reader_close (group_reader);
        if(Keep_groups && leaf_mult > 0) tree->Fill();
        else if (!Keep_groups)
        {
          for(int ev = 0; ev < leaf_mult; ev++)
          {
            leaf_label_2 = leaf_label.at(ev);
            leaf_nrj_2 = leaf_nrj.at(ev);
            if(energy2) leaf_nrj2_2 = leaf_nrj2.at(ev);
            if(energy3) leaf_nrj3_2 = leaf_nrj3.at(ev);
            if(energy4) leaf_nrj4_2 = leaf_nrj4.at(ev);
            leaf_t_2 = leaf_t.at(ev);
            leaf_pu_2= leaf_pu.at(ev);
            //cout << "Label " << leaf_label_2 << " @ " << leaf_t_2 << " with " << leaf_nrj_2 << endl;
            tree->Fill();
          }
        }
      }

      // If Data are ungrouped
      if (alias != GROUP_TYPE_ALIAS)
      {
        leaf_nrj.clear(); //energy
        if(energy2) leaf_nrj2.clear(); //used if QDC2
        if(energy3) leaf_nrj3.clear(); //used if QDC3
        if(energy4) leaf_nrj4.clear(); //used if QDC4
        leaf_t.clear(); //time
        leaf_label.clear(); //label
        leaf_pu.clear(); //pile-up detection
        leaf_mult = 0;

        //clock_g = faster_data_hr_clock_ns(data)*1000;  //Time in ps

        //if (myID[faster_data_label(data)] != 0) cout << endl << "READ faster label " << faster_data_label(data) << " associated to type " << myID[faster_data_label(data)] << endl;

        //RF
        if ( myID[faster_data_label(data)] == 7)
        {
          faster_data_load (data, &rf);
          // Getting absolute time for the hit
          clock_g = faster_data_hr_clock_ns(data)*1000;  //Time in ps
          //if (LABR_qdc.q1_saturated!=1) // reject saturation
          {
            leaf_label.push_back(faster_data_label (data));
            leaf_nrj.push_back(rf_period_ns(rf)*1000);
            if(energy2) leaf_nrj2.push_back(0);
            if(energy3) leaf_nrj3.push_back(0);
            if(energy4) leaf_nrj4.push_back(0);
            leaf_t.push_back((ULong64_t)clock_g);
            leaf_pu.push_back(false);
            leaf_mult++;
          }
        }

        //TRef
        if(myID[faster_data_label(data)] == 8)
        {
          faster_data_load (data, &TREF_tdc);
          //clock_g = 1000.*(faster_data_clock_ns (data) + tref_conv_dt_ns (TREF_tdc.tdc));
          clock_g = faster_data_hr_clock_ns(data)*1000;
          leaf_pu.push_back(false);
          leaf_label.push_back(faster_data_label (data));
          leaf_nrj.push_back(0);
          if(energy2) leaf_nrj2.push_back(0);
          if(energy3) leaf_nrj3.push_back(0);
          if(energy4) leaf_nrj4.push_back(0);
          leaf_t.push_back((ULong64_t)clock_g);
          leaf_mult++;
        }

        //trapezoidal
        else if ( myID[faster_data_label(data)] == 5)
        {

          faster_data_load (data, &trapez_adc);
          // Getting absolute time for the hit

          //clock_g = 1000*((Double_t)faster_data_clock_ns (data) + (Double_t)tref_conv_dt_ns (trapez_adc.tdc)*4.);
          //cout << clock_g << endl;
          clock_g = faster_data_hr_clock_ns(data)*1000;
          //cout << clock_g << endl;

          if (trapez_adc.pileup!=1) leaf_pu.push_back(false);
          else leaf_pu.push_back(true);// reject pile up

          leaf_label.push_back(faster_data_label (data));
          leaf_nrj.push_back(trapez_adc.measure);
          if(energy2) leaf_nrj2.push_back(0);
          if(energy3) leaf_nrj3.push_back(0);
          if(energy4) leaf_nrj4.push_back(0);
          leaf_t.push_back((ULong64_t)clock_g);
          leaf_mult++;
        }

        //QDC1
        else if ( myID[faster_data_label(data)] == 1)
        {
          faster_data_load (data, &QDC1_qdc);

          //clock_g = 1000.*(faster_data_clock_ns (data) + tref_conv_dt_ns (QDC1_qdc.tdc));
          clock_g = faster_data_hr_clock_ns(data)*1000;
          if (QDC1_qdc.q1_saturated!=1) leaf_pu.push_back(false); // reject saturation
          else leaf_pu.push_back(true);

          leaf_label.push_back(faster_data_label (data));
          leaf_nrj.push_back(QDC1_qdc.q1);
          if(energy2) leaf_nrj2.push_back(0);
          if(energy3) leaf_nrj3.push_back(0);
          if(energy4) leaf_nrj4.push_back(0);
          leaf_t.push_back((ULong64_t)clock_g);
          leaf_mult++;
        }

        //QDC2
        else if ( myID[faster_data_label (data)] == 2)
        {
          faster_data_load (data, &QDC2_qdc);
          // Getting absolute time for the hit
          //clock_g = (faster_data_clock_ns (data) + tref_conv_dt_ns (QDC2_qdc.tdc))*1000.;
          clock_g = faster_data_hr_clock_ns(data)*1000;
          if (QDC2_qdc.q1_saturated!=1) leaf_pu.push_back(false);
          else leaf_pu.push_back(true); // reject saturation

          leaf_label.push_back(faster_data_label (data));
          leaf_nrj.push_back(QDC2_qdc.q1);
          leaf_nrj2.push_back(QDC2_qdc.q2);
          if(energy3) leaf_nrj3.push_back(0);
          if(energy4) leaf_nrj4.push_back(0);
          leaf_t.push_back((ULong64_t)clock_g);
          leaf_mult++;
        }
        //QDC3
        else if ( myID[faster_data_label (data)] == 3)
        {
          faster_data_load (data, &QDC3_qdc);
          // Getting absolute time for the hit
          //clock_g = (faster_data_clock_ns (data) + tref_conv_dt_ns (QDC3_qdc.tdc))*1000.;
          clock_g = faster_data_hr_clock_ns(data)*1000;
          if (QDC3_qdc.q1_saturated!=1) leaf_pu.push_back(false);
          else leaf_pu.push_back(true); // reject saturation

          leaf_label.push_back(faster_data_label (data));
          leaf_nrj.push_back(QDC3_qdc.q1);
          leaf_nrj2.push_back(QDC3_qdc.q2);
          leaf_nrj3.push_back(QDC3_qdc.q3);
          if(energy4) leaf_nrj4.push_back(0);
          leaf_t.push_back((ULong64_t)clock_g);
          leaf_mult++;
        }

        //QDC4
        else if ( myID[faster_data_label (data)] == 4)
        {
          faster_data_load (data, &QDC4_qdc);
          // Getting absolute time for the hit
          //clock_g = (faster_data_clock_ns (data) + tref_conv_dt_ns (QDC4_qdc.tdc))*1000.;
          clock_g = faster_data_hr_clock_ns(data)*1000;
          if (QDC4_qdc.q1_saturated!=1) leaf_pu.push_back(false);
          else leaf_pu.push_back(true); // reject saturation

          leaf_label.push_back(faster_data_label (data));
          leaf_nrj.push_back(QDC4_qdc.q1);
          leaf_nrj2.push_back(QDC4_qdc.q2);
          leaf_nrj3.push_back(QDC4_qdc.q3);
          leaf_nrj4.push_back(QDC4_qdc.q4);
          leaf_t.push_back((ULong64_t)clock_g);
          leaf_mult++;

        }

        //CRRC4
        else if ( myID[faster_data_label (data)] == 6)
        {
          faster_data_load (data, &crrc4_adc);
          // Getting absolute time for the hit
          //clock_g = 1000.*(faster_data_clock_ns (data) + tref_conv_dt_ns (crrc4_adc.delta_t));
          clock_g = faster_data_hr_clock_ns(data)*1000;
          leaf_pu.push_back(false);
          leaf_label.push_back(faster_data_label (data));
          leaf_nrj.push_back(crrc4_adc.measure);
          if(energy2) leaf_nrj2.push_back(0);
          if(energy3) leaf_nrj3.push_back(0);
          if(energy4) leaf_nrj4.push_back(0);
          leaf_t.push_back((ULong64_t)clock_g);
          leaf_mult++;
        }


        if(leaf_mult >0)
        {
          if(Keep_groups) tree->Fill();
          else if (!Keep_groups)
          {
            for(int ev = 0; ev < leaf_mult; ev++)
            {
              leaf_label_2 = leaf_label.at(ev);
              leaf_nrj_2 = leaf_nrj.at(ev);
              if(energy2) leaf_nrj2_2 = leaf_nrj2.at(ev);
              if(energy3) leaf_nrj3_2 = leaf_nrj3.at(ev);
              if(energy4) leaf_nrj4_2 = leaf_nrj4.at(ev);
              leaf_t_2 = leaf_t.at(ev);
              leaf_pu_2= leaf_pu.at(ev);
              //cout << "Label " << leaf_label_2 << " @ " << leaf_t_2 << " with " << leaf_nrj_2 << endl;
              tree->Fill();
            }
          }
        }
      }
    }
    //We save the output file
    faster_file_reader_close (reader);
    root_file->cd();
    tree->AutoSave();
    root_file->Close();

  }
}

//---------------------------------------------------------------------------//
//                                                                           //
//                  Fonction pour manipuler les datas                        //
//                                                                           //
//---------------------------------------------------------------------------//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Function to sort a vector and apply the same permutation to other vectors
template <typename T, typename Compare>
void getSortPermutation(
    std::vector<unsigned>& out,
    const std::vector<T>& v,
    Compare compare = std::less<T>())
    {
      out.resize(v.size());
      std::iota(out.begin(), out.end(), 0);

      std::sort(out.begin(), out.end(),
          [&](unsigned i, unsigned j){ return compare(v[i], v[j]); });
    }

template <typename T>
void applyPermutation(
    const std::vector<unsigned>& order,
    std::vector<T>& t)
{
    assert(order.size() == t.size());
    std::vector<T> st(t.size());
    for(unsigned i=0; i<t.size(); i++)
    {
        st[i] = t[order[i]];
    }
    t = st;
}

template <typename T, typename... S>
void applyPermutation(
    const std::vector<unsigned>& order,
    std::vector<T>& t,
    std::vector<S>&... s)
{
    applyPermutation(order, t);
    applyPermutation(order, s...);
}

// sort multiple vectors using the criteria of the first one
template<typename T, typename Compare, typename... SS>
void sortVectors(
    const std::vector<T>& t,
    Compare comp,
    std::vector<SS>&... ss)
{
    std::vector<unsigned> order;
    getSortPermutation(order, t, comp);
    applyPermutation(order, ss...);
}

// make less verbose for the usual ascending order
template<typename T, typename... SS>
void sortVectorsAscending(
    const std::vector<T>& t,
    std::vector<SS>&... ss)
{
    sortVectors(t, std::less<T>(), ss...);
}


//---------------------------------------------------------------------------//
//                                                                           //
//          Fonction pour trier les evenements en temps sur des datas        //
//                           non calibrees                                   //
//---------------------------------------------------------------------------//
int sorterF(string FileName, const string & temp_directory, const string & output_directory, map<int,int> myID)
{
  //Boolean variables used if there are QDC with several integration gates
  bool energy2 = false;
  bool energy3 = false;
  bool energy4 = false;
  int nbrdetect = 0;
  for(const auto & elem : myID)
  {
    if(elem.second == 2) //QDC2
    {
      energy2 = true;
    }
    else if(elem.second == 3) //QDC3
    {
      energy2 = true;
      energy3 = true;
    }
    else if(elem.second == 4)//QDC4
    {
      energy2 = true;
      energy3 = true;
      energy4 = true;
    }
    nbrdetect++;
  }

  string DATAFILENAME;
  string ROOTFILENAME;
  if (FileName.size()>1 )
  {
    DATAFILENAME = temp_directory;
    //Automatically change the path and the extension of the filename for the output file
    int last_slash =  FileName.find_last_of('/');
    FileName.erase(FileName.end()-5, FileName.end());
    FileName.erase(FileName.begin(), FileName.begin()+last_slash+1);
    FileName.insert(FileName.end()-4,'R');
    FileName.insert(FileName.end()-4,'_');
    DATAFILENAME += FileName+".root";
    ROOTFILENAME = output_directory+FileName+".root";
  }

  //  print infos
  printf ("\n");
  printf ("  Sorting :\n");
  printf ("     - read faster file '%s'\n", DATAFILENAME.c_str());
  printf ("     - output those data to root file '%s'\n", ROOTFILENAME.c_str());
  printf ("\n");

  // Declaration of new tree variables
  nrj_type enrj1,enrj2,enrj3,enrj4;
  label_type index;
  pu_type pileup;
  time_type tm;

  // Declaration of old tree variables
  nrj_type enrj1_2,enrj2_2,enrj3_2,enrj4_2;
  label_type index_2;
  pu_type pileup_2;
  time_type tm_2;

  // Another Timer
  TStopwatch timer2;

  // Autres variables
  Long64_t chainentries;
  int buffersize = 10000000;  // Size of the buffer that will be manipulated
  Long64_t nbrsave = 0;
  Long64_t nbrsave2=0;

  // Dynamic allocation of the pointeurs that will contain the data read from the file
  vector <nrj_type> buffer_enrj;
  vector <nrj_type> buffer_enrj2;
  vector <nrj_type> buffer_enrj3;
  vector <nrj_type> buffer_enrj4;
  vector <time_type> buffer_tm;
  vector <label_type> buffer_idx;
  vector <pu_type> buffer_pu;
  vector <time_type> tampon_tm;
  vector <nrj_type> tampon_enrj;
  vector <nrj_type> tampon_enrj2;
  vector <nrj_type> tampon_enrj3;
  vector <nrj_type> tampon_enrj4;
  vector <label_type> tampon_idx;
  vector <pu_type> tampon_pu;

  //-----------------------------------------------------------------------------//
  //                                                                             //
  //         Part that is used to define the parameters that belong to the       //
  //         experiment: definition of an index to recognize the differents      //
  //         parameters.                                                         //
  //                                                                             //
  //-----------------------------------------------------------------------------//

  //Starting of the chronometer
  timer2.Reset();
  timer2.Start();

  std::cout << "Opening FILE for sorting" << std::endl;

  // Je declare un TFile pour aller charger l'arbre narval
  TFile *Infile_tree = TFile::Open(DATAFILENAME.data(),"READ");

  std::cout << "Loading old tree " << std::endl;

  // Je declare un TTree pour stocker celui qui est dans le fichier
  TTree *old_tree = (TTree*) Infile_tree->Get("DataTree");

  //std::cout << "Loading Tree ..." << std::endl;
  //old_tree -> SetBranchAddress("label",&index1);   // Detector number
  //old_tree -> SetBranchAddress("nrj",&enrj1);      // Energy coded on the channel
  //old_tree -> SetBranchAddress("time",&tm1);       // Time of the hit                  // Time of the hit
  old_tree->SetBranchAddress ("label", &index_2);
  old_tree->SetBranchAddress ("nrj", &enrj1_2);
  if(energy2) old_tree->SetBranchAddress ("nrj2",&enrj2_2);
  if(energy3) old_tree->SetBranchAddress ("nrj3",&enrj3_2);
  if(energy4) old_tree->SetBranchAddress ("nrj4",&enrj4_2);
  old_tree->SetBranchAddress ("time", &tm_2);
  old_tree->SetBranchAddress ("pileup",&pileup_2);
  //std::cout << "Loaded .." << std::endl;
  // Mon vieil arbre est charge...
  std::cout << "Old tree loaded" << std::endl;

  // Je declare le fichier
  TFile *Outfile_tree = new TFile(ROOTFILENAME.data(),"RECREATE");

  // Je declare le nouvel arbre...
  // Declaration of oak, that is the tree containing all events
  std::cout<< "Building tree  " <<std::endl;
  TTree *oak = new TTree("DataTree",ROOTFILENAME.data());

  // Declaration of Branches that will contain the data
  oak->Branch ("label", &index);
  oak->Branch ("nrj", &enrj1);
  if(energy2) oak->Branch ("nrj2",&enrj2);
  if(energy3) oak->Branch ("nrj3",&enrj3);
  if(energy4) oak->Branch ("nrj4",&enrj4);
  oak->Branch ("time", &tm);
  oak->Branch ("pileup",&pileup);

  std::cout << "New Seed burried " << std::endl;
  std::cout << "and in the memory another tree will grow..."<< std::endl ;


  // I read the old tree
  chainentries = old_tree->GetEntries();
  std::cout << "There are " << (ULong64_t) chainentries << " entries to read and sort" << std::endl;

  int nbrevent2 = 0;

  // Je charge le premier tampon
  for(int jj = 0; jj < nbrdetect ; jj++)
  {
    old_tree->GetEntry(jj);


    // et je le remplis
    tampon_tm.push_back(tm_2);
    tampon_pu.push_back(pileup_2);
    tampon_enrj.push_back(enrj1_2);
    if(energy2)tampon_enrj2.push_back(enrj2_2);
    if(energy3)tampon_enrj3.push_back(enrj3_2);
    if(energy4)tampon_enrj4.push_back(enrj4_2);
    tampon_idx.push_back(index_2);
    nbrevent2++;
  }

//std::cout << "sorting the tampon " << std::endl;

// je tri le premier tampon
// I start by sorting this tampon (depending on the nubmer of NRJ values)
  if(energy4) sortVectors(tampon_tm, less<Double_t>(), tampon_tm, tampon_enrj, tampon_enrj2, tampon_enrj3, tampon_enrj4, tampon_pu, tampon_idx);
  else if (energy3 && !energy4) sortVectors(tampon_tm, less<Double_t>(), tampon_tm, tampon_enrj, tampon_enrj2, tampon_enrj3, tampon_pu, tampon_idx);
  else if (energy2 && !energy3 && !energy4) sortVectors(tampon_tm, less<Double_t>(), tampon_tm, tampon_enrj, tampon_enrj2, tampon_pu, tampon_idx);
  else if (!energy2 && !energy3 && !energy4) sortVectors(tampon_tm, less<Double_t>(), tampon_tm, tampon_enrj, tampon_pu, tampon_idx);

//std::cout << "Tampon sorted going for buffer " << std::endl;



Long64_t nbrevent = 0;
Long64_t i = nbrdetect;
while (i < chainentries)
{
  // I get the entry
  old_tree->GetEntry(i);

  //std::cout << "nbreevent " << nbrevent << std::endl;
  // Put it in the tables
  buffer_tm.push_back(tm_2);
  buffer_pu.push_back(pileup_2);
  buffer_enrj.push_back(enrj1_2);
  if(energy2)buffer_enrj2.push_back(enrj2_2);
  if(energy3)buffer_enrj3.push_back(enrj3_2);
  if(energy4)buffer_enrj4.push_back(enrj4_2);
  buffer_idx.push_back(index_2);

  nbrevent++;
  //v
  std::cout.precision(10);

  // Buffer is filled: I put it into order and go
  if(nbrevent == buffersize)
  {
    // First I concatenate the buffer to the tampon
    tampon_tm.insert( tampon_tm.end(), buffer_tm.begin(), buffer_tm.end() );
    tampon_pu.insert( tampon_pu.end(), buffer_pu.begin(), buffer_pu.end() );
    tampon_enrj.insert( tampon_enrj.end(), buffer_enrj.begin(), buffer_enrj.end() );
    if(energy2)tampon_enrj2.insert( tampon_enrj2.end(), buffer_enrj2.begin(), buffer_enrj2.end() );
    if(energy3)tampon_enrj3.insert( tampon_enrj3.end(), buffer_enrj3.begin(), buffer_enrj3.end() );
    if(energy4)tampon_enrj4.insert( tampon_enrj4.end(), buffer_enrj4.begin(), buffer_enrj4.end() );
    tampon_idx.insert( tampon_idx.end(), buffer_idx.begin(), buffer_idx.end() );

    //std::cout << "preparation buffer sans coden" << std::endl;
    // Tri des buffers
    // Now the two have been concatenated I can sort them
    if(energy4) sortVectors(tampon_tm, less<Double_t>(), tampon_tm, tampon_enrj, tampon_enrj2, tampon_enrj3, tampon_enrj4, tampon_pu, tampon_idx);
    else if (energy3 && !energy4) sortVectors(tampon_tm, less<Double_t>(), tampon_tm, tampon_enrj, tampon_enrj2, tampon_enrj3, tampon_pu, tampon_idx);
    else if (energy2 && !energy3 && !energy4) sortVectors(tampon_tm, less<Double_t>(), tampon_tm, tampon_enrj, tampon_enrj2, tampon_pu, tampon_idx);
    else if (!energy2 && !energy3 && !energy4) sortVectors(tampon_tm, less<Double_t>(), tampon_tm, tampon_enrj, tampon_pu, tampon_idx);

    //std::cout << "Apres tri " << std::endl;
    // Now the buffer is sorted I will be able to store it in the TTree
    for(int ll = 0; ll < buffersize; ll++)
    {
      // I get the information back from the buffer
      index = tampon_idx.at(ll);
      nbrsave2++;
      if(index!=0)nbrsave++;
      enrj1 = tampon_enrj.at(ll);
      if(energy2) enrj2 = tampon_enrj2.at(ll);
      if(energy3) enrj3 = tampon_enrj3.at(ll);
      if(energy4) enrj4 = tampon_enrj4.at(ll);
      tm = tampon_tm.at(ll);
      pileup = tampon_pu.at(ll);
      if(tampon_tm.at(ll) < 1.e-8) tm = 0;

      // And it goes to memory
      if (index != 0) {oak -> Fill();}
    }


    // Now I need to restore the tampon the last nbrdetect hit in the concatenated vector
    tampon_tm.erase(tampon_tm.begin(), tampon_tm.begin()+buffersize);
    tampon_pu.erase(tampon_pu.begin(), tampon_pu.begin()+buffersize);
    tampon_enrj.erase(tampon_enrj.begin(), tampon_enrj.begin()+buffersize);
    if(energy2) tampon_enrj2.erase(tampon_enrj2.begin(), tampon_enrj2.begin()+buffersize);
    if(energy3) tampon_enrj3.erase(tampon_enrj3.begin(), tampon_enrj3.begin()+buffersize);
    if(energy4) tampon_enrj4.erase(tampon_enrj4.begin(), tampon_enrj4.begin()+buffersize);
    tampon_idx.erase(tampon_idx.begin(), tampon_idx.begin()+buffersize);
    // Reset event number to start filling the buffer from the first slot
    nbrevent = 0;

    // Now I clear the buffer vectors to fill them again
    buffer_tm.clear();
    buffer_pu.clear();
    buffer_enrj.clear();
    if(energy2)buffer_enrj2.clear();
    if(energy3)buffer_enrj3.clear();
    if(energy4)buffer_enrj4.clear();
    buffer_idx.clear();
  }
  //std::cout << " Coucou " << i << std::endl;

  // Get to the next event
  i++;
}

// Normally I still have a buffer and a tampon in memory.
// I handle these last events
// First I concatenate the buffer to the tampon
tampon_tm.insert( tampon_tm.end(), buffer_tm.begin(), buffer_tm.end() );
tampon_pu.insert( tampon_pu.end(), buffer_pu.begin(), buffer_pu.end() );
tampon_enrj.insert( tampon_enrj.end(), buffer_enrj.begin(), buffer_enrj.end() );
if(energy2)tampon_enrj2.insert( tampon_enrj2.end(), buffer_enrj2.begin(), buffer_enrj2.end() );
if(energy3)tampon_enrj3.insert( tampon_enrj3.end(), buffer_enrj3.begin(), buffer_enrj3.end() );
if(energy4)tampon_enrj4.insert( tampon_enrj4.end(), buffer_enrj4.begin(), buffer_enrj4.end() );
tampon_idx.insert( tampon_idx.end(), buffer_idx.begin(), buffer_idx.end() );

// Now the two have been concatenated I can sort them
if(energy4) sortVectors(tampon_tm, less<Double_t>(), tampon_tm, tampon_enrj, tampon_enrj2, tampon_enrj3, tampon_enrj4, tampon_pu, tampon_idx);
else if (energy3 && !energy4) sortVectors(tampon_tm, less<Double_t>(), tampon_tm, tampon_enrj, tampon_enrj2, tampon_enrj3, tampon_pu, tampon_idx);
else if (energy2 && !energy3 && !energy4) sortVectors(tampon_tm, less<Double_t>(), tampon_tm, tampon_enrj, tampon_enrj2, tampon_pu, tampon_idx);
else if (!energy2 && !energy3 && !energy4) sortVectors(tampon_tm, less<Double_t>(), tampon_tm, tampon_enrj, tampon_pu, tampon_idx);

// Now they are sorted I must save data to the new TTree
for(int ll = 0; ll < tampon_tm.size(); ll++)
{
  // I get the information back from the buffer
  index = tampon_idx.at(ll);
  nbrsave2++;
  if(index!=0) nbrsave++;
  enrj1 = tampon_enrj.at(ll);
  if(energy2) enrj2 = tampon_enrj2.at(ll);
  if(energy3) enrj3 = tampon_enrj3.at(ll);
  if(energy4) enrj4 = tampon_enrj4.at(ll);
  tm = tampon_tm.at(ll);
  pileup = tampon_pu.at(ll);
  if(tampon_tm.at(ll) < 1.e-8) tm = 0;

  // And it goes to memory
  if (index != 0) {oak -> Fill();}
}

// Mon nouvel arbre est rempli
std::cout << "I sorted " << (double)nbrsave2 << " event in the new file " << std::endl;
std::cout << "I stored " << (double)nbrsave << " event in the new file " << std::endl;
// Je le sauvegarde
Outfile_tree -> cd();
oak -> AutoSave();
old_tree->Delete();
Outfile_tree -> Close();
Infile_tree->Close();

// Printing of chronometer measurement
timer2.Stop();
Double_t rtime2 = timer2.RealTime();
Double_t ctime2 = timer2.CpuTime();
std::cout << std::endl;
std::cout << "End of Time spectrum creation" << std::endl;
std::cout << "# RealTime=" << rtime2 << " seconds, CpuTime="<< ctime2 << " seconds" <<std::endl;
std::cout << std::endl;

return 1;
}


void init_parameters(int & nOfThreads, string & output_directory, string & input_file, string & temp_directory, bool & kGroups)
{
  ifstream parameters_file("/local/home/dt269135/Codes/faster2root/parameters.dat");
  //Set the default parameters here:
  nOfThreads = 1;
  output_directory = "~/";
  input_file = "/local/home/dt269135/Codes/faster2root/filename.dat";
  kGroups = false;
  //------------------------------//

  if(parameters_file.is_open())
  {
    while (parameters_file.good())
    {
      string oneline;
      getline(parameters_file, oneline);
      istringstream is(oneline);
      string temp1 = "NULL";
      is>>temp1;
      if(temp1 == "final_directory:")
      {
        is >> output_directory;
	temp_directory = output_directory;
      }
      else if(temp1 == "number_of_threads:")
      {
        is >> nOfThreads;
      }
      else if(temp1 == "input_file:")
      {
        is >> input_file;
      }
      else if(temp1 == "keep_groups:")
      {
        string bb;
        is >> bb;
        if(bb == "true")
        {
          kGroups = true;
        }
      }
      else if(temp1 != "NULL")
      {
        cout << "WARNING: Unknown parameter " << temp1 << " will be ignored !" << endl;
      }
    }
    cout << endl << "-----------Parameters--------------" << endl;
    cout << "Output directory set to : " << output_directory << endl;
    cout << "Number of threads set to : " << nOfThreads << endl;
    cout << "Input file set to : " << input_file << endl;
    cout << "Ungrouping of data set to : " << kGroups << endl << endl;
  }
  else
  {
    cout << "WARNING: No parameters file loaded !" << endl;
  }
  parameters_file.close();
}
