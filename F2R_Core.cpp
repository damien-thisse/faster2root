#include "F2R_Core.h"
#include "F2R_Parameters.h"

//-----fasterac includes-----//
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
//---------------------------//

//-----root includes-----//
#include "TFile.h"
#include "TTree.h"


using namespace std;

void Convert(string filename)
{
  F2R_Parameters& parameters = F2R_Parameters::getInstance();

  bool energy2 = false, energy3 = false, energy4 = false;
  int maxEnergyBranchNeeded = parameters.getMaxEnergyBranchNeeded();
  if(maxEnergyBranchNeeded > 1){
      energy2 = true;
  }
  if(maxEnergyBranchNeeded > 2){
      energy3 = true;
  }
  if(maxEnergyBranchNeeded > 3){
      energy4 = true;
  }

  bool keepGroups = parameters.getKeepGroups();

  map<int, int> myID = parameters.getIDTable();

  //Definition of the ROOT file name from the FASTER file name
  string DATAFILENAME;
  string ROOTFILENAME;
  int last_slash =  filename.find_last_of('/');
  DATAFILENAME = filename;
  filename.erase(filename.end()-5, filename.end()); // removes .fast
  filename.erase(filename.begin(), filename.begin()+last_slash+1); //removes the path
  // filename.insert(filename.end()-4,'R');
  // filename.insert(filename.end()-4,'_');
  ROOTFILENAME = parameters.getOutputDirectory()+filename+".root"; //add output path and .root extension

  //FASTER VARIABLE DECLARATION
  // File reader
  faster_file_reader_p reader;

  // Data handler
  faster_data_p data;
  unsigned char alias;
  unsigned short label;
  time_type clock_g;
  time_type clock_a;

  // Group data handler
      //group
  faster_buffer_reader_p group_reader;
  char                   group_buffer [15000];
  unsigned short         lsize;
  faster_data_p          group_data;
      //inside group
  faster_buffer_reader_p inside_group_reader;
  char                   inside_group_buffer [15000];
  unsigned short         inside_lsize;
  faster_data_p          inside_group_data;

  // FASTER structures for handling data
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

  //ROOT VARIABLE DECLARATION
  //Branches of the output ROOT tree
      //If keep group is true
  std::vector<nrj_type> leaf_nrj; //energy
  std::vector<nrj_type> leaf_nrj2; //used if QDC2
  std::vector<nrj_type> leaf_nrj3; //used if QDC3
  std::vector<nrj_type> leaf_nrj4; //used if QDC4
  std::vector<time_type> leaf_t; //time
  std::vector<label_type> leaf_label; //label
  std::vector<Bool_t> leaf_pu; //pile-up detection
  UShort_t leaf_mult;
  time_type leaf_group_time;
      //If data are ungrouped
  nrj_type leaf_nrj_2; //energy
  nrj_type leaf_nrj2_2; //used if QDC2
  nrj_type leaf_nrj3_2; //used if QDC3
  nrj_type leaf_nrj4_2; //used if QDC4
  time_type leaf_t_2; //time
  label_type leaf_label_2; //label
  Bool_t leaf_pu_2; //pile-up detection

  TFile* root_file = new TFile(ROOTFILENAME.data(), "RECREATE");
  TTree* root_tree = new TTree ("DataTree", "Root tree converted using F2R code");
  if(keepGroups){
      root_tree->Branch("label", &leaf_label);
      root_tree->Branch("energy", &leaf_nrj);
      root_tree->Branch("time", &leaf_t);
      root_tree->Branch("pileup",&leaf_pu);
      if(energy2){
        root_tree->Branch("energy2",&leaf_nrj2);
      }
      if(energy3){
        root_tree->Branch("energy3",&leaf_nrj3);
      }
      if(energy4){
        root_tree->Branch("energy4",&leaf_nrj4);
      }
      root_tree->Branch("multiplicity",&leaf_mult);
      root_tree->Branch("GroupTimeStamp",&leaf_group_time);
  }

  else{
      root_tree->Branch("label", &leaf_label_2);
      root_tree->Branch("nrj", &leaf_nrj_2);
      root_tree->Branch("time", &leaf_t_2);
      root_tree->Branch("pileup",&leaf_pu_2);
      if(energy2){
        root_tree->Branch("energy2",&leaf_nrj2_2);
      }
        if(energy3){
      root_tree->Branch("energy3",&leaf_nrj3_2);
      }
        if(energy4){
      root_tree->Branch("energy4",&leaf_nrj4_2);
      }
  }

  //END OF SETUP, NOW TOUGH THINGS ARE STARTING !!

  //while there are data to be read...
  while ((data = faster_file_reader_next(reader)) != NULL){
      alias = faster_data_type_alias(data);

            // If data are grouped
    if (alias == GROUP_TYPE_ALIAS){ //in this condition, the data are grouped
      leaf_nrj.clear();
      if(energy2) leaf_nrj2.clear();
      if(energy3) leaf_nrj3.clear();
      if(energy4) leaf_nrj4.clear();
      leaf_t.clear(); 
      leaf_label.clear(); 
      leaf_pu.clear();
      leaf_mult = 0;
      lsize = faster_data_load(data, group_buffer); // get group data
      group_reader = faster_buffer_reader_open(group_buffer, lsize); //open group reader

      while ((group_data = faster_buffer_reader_next (group_reader)) != NULL){ // we loop on all the data inside the group
        clock_g = faster_data_hr_clock_ns(group_data)*1000; //Time in ps

          leaf_group_time=clock_g;
        if(faster_data_label(group_data)>2999) //if it is a subgroup, then we have to loop on all the elements is this subgroup
        {
          inside_lsize = faster_data_load (group_data, inside_group_buffer);
          inside_group_reader = faster_buffer_reader_open (inside_group_buffer, inside_lsize);

          while((inside_group_data = faster_buffer_reader_next (inside_group_reader)) != NULL)
          {
            clock_g = faster_data_hr_clock_ns(inside_group_data)*1000; //Time in ps

            //RF
            if ( myID[faster_data_label(inside_group_data)] == 7){
              faster_data_load (inside_group_data, &rf);

              // Getting absolute time for the hit
              clock_g = faster_data_hr_clock_ns(inside_group_data)*1000;  //Time in ps
              leaf_label.push_back(faster_data_label (inside_group_data));
              leaf_nrj.push_back(rf_period_ns(rf)*1000); //in case of RF signal, the energy branch contains the period 
              if(energy2) leaf_nrj2.push_back(0);
              if(energy3) leaf_nrj3.push_back(0);
              if(energy4) leaf_nrj4.push_back(0);
              leaf_t.push_back((time_type)clock_g);
              leaf_pu.push_back(false);
              leaf_mult++;           
            }

            //TRef (== QDC0)
            if(myID[faster_data_label(inside_group_data)] == 8){
              faster_data_load (inside_group_data, &TREF_tdc);
              clock_g = faster_data_hr_clock_ns(inside_group_data)*1000;
              leaf_pu.push_back(false);
              leaf_label.push_back(faster_data_label (inside_group_data));
              leaf_nrj.push_back(0);
              if(energy2) leaf_nrj2.push_back(0);
              if(energy3) leaf_nrj3.push_back(0);
              if(energy4) leaf_nrj4.push_back(0);
              leaf_t.push_back((time_type)clock_g);
              leaf_mult++;
            }

            //trapezoidal
            else if ( myID[faster_data_label(inside_group_data)] == 5){
              faster_data_load (inside_group_data, &trapez_adc);
              // Getting absolute time for the hit
              clock_g = faster_data_hr_clock_ns(inside_group_data)*1000;
              if (trapez_adc.pileup!=1) leaf_pu.push_back(false);
              else leaf_pu.push_back(true);
              leaf_label.push_back(faster_data_label (inside_group_data));
              leaf_nrj.push_back(trapez_adc.measure);
              if(energy2) leaf_nrj2.push_back(0);
              if(energy3) leaf_nrj3.push_back(0);
              if(energy4) leaf_nrj4.push_back(0);
              leaf_t.push_back((time_type)clock_g);
              leaf_mult++;
            }

            //QDC1
            else if ( myID[faster_data_label(inside_group_data)] == 1){
              faster_data_load (inside_group_data, &QDC1_qdc);
              clock_g = faster_data_hr_clock_ns(inside_group_data)*1000;
              if (QDC1_qdc.q1_saturated!=1) leaf_pu.push_back(false);
              else leaf_pu.push_back(true);
              leaf_label.push_back(faster_data_label (inside_group_data));
              leaf_nrj.push_back(QDC1_qdc.q1);
              if(energy2) leaf_nrj2.push_back(0);
              if(energy3) leaf_nrj3.push_back(0);
              if(energy4) leaf_nrj4.push_back(0);
              leaf_t.push_back((time_type)clock_g);
              leaf_mult++;
            }

            //QDC2
            else if ( myID[faster_data_label (inside_group_data)] == 2){
              faster_data_load (inside_group_data, &QDC2_qdc);
              clock_g = faster_data_hr_clock_ns(inside_group_data)*1000;
              if (QDC2_qdc.q1_saturated!=1) leaf_pu.push_back(false);
              else leaf_pu.push_back(true);
              leaf_label.push_back(faster_data_label (inside_group_data));
              leaf_nrj.push_back(QDC2_qdc.q1);
              leaf_nrj2.push_back(QDC2_qdc.q2);
              if(energy3) leaf_nrj3.push_back(0);
              if(energy4) leaf_nrj4.push_back(0);
              leaf_t.push_back((time_type)clock_g);
              leaf_mult++;
            }
            //QDC3
            else if ( myID[faster_data_label (inside_group_data)] == 3){
              faster_data_load (inside_group_data, &QDC3_qdc);
              clock_g = faster_data_hr_clock_ns(inside_group_data)*1000;
              if (QDC3_qdc.q1_saturated!=1) leaf_pu.push_back(false);
              else leaf_pu.push_back(true); // reject saturation
              leaf_label.push_back(faster_data_label (inside_group_data));
              leaf_nrj.push_back(QDC3_qdc.q1);
              leaf_nrj2.push_back(QDC3_qdc.q2);
              leaf_nrj3.push_back(QDC3_qdc.q3);
              if(energy4) leaf_nrj4.push_back(0);
              leaf_t.push_back((time_type)clock_g);
              leaf_mult++;
            }

            //QDC4
            else if ( myID[faster_data_label (inside_group_data)] == 4){
              faster_data_load (inside_group_data, &QDC4_qdc);
              clock_g = faster_data_hr_clock_ns(inside_group_data)*1000;
              if (QDC4_qdc.q1_saturated!=1) leaf_pu.push_back(false);
              else leaf_pu.push_back(true); // reject saturation
              leaf_label.push_back(faster_data_label (inside_group_data));
              leaf_nrj.push_back(QDC4_qdc.q1);
              leaf_nrj2.push_back(QDC4_qdc.q2);
              leaf_nrj3.push_back(QDC4_qdc.q3);
              leaf_nrj4.push_back(QDC4_qdc.q4);
              leaf_t.push_back((time_type)clock_g);
              leaf_mult++;
            }

            //CRRC4
            else if ( myID[faster_data_label (inside_group_data)] == 6){
              faster_data_load (inside_group_data, &crrc4_adc);
              clock_g = faster_data_hr_clock_ns(inside_group_data)*1000;
              leaf_pu.push_back(false);
              leaf_label.push_back(faster_data_label (inside_group_data));
              leaf_nrj.push_back(crrc4_adc.measure);
              if(energy2) leaf_nrj2.push_back(0);
              if(energy3) leaf_nrj3.push_back(0);
              if(energy4) leaf_nrj4.push_back(0);
              leaf_t.push_back((time_type)clock_g);
              leaf_mult++;
            }
          }
          faster_buffer_reader_close(inside_group_reader); //closing the groupe reader
        }

        //RF
          if( myID[faster_data_label(group_data)] == 7){
              faster_data_load (group_data, &rf);

              // Getting absolute time for the hit
              clock_g = faster_data_hr_clock_ns(group_data)*1000;  //Time in ps
              leaf_label.push_back(faster_data_label (group_data));
              leaf_nrj.push_back(rf_period_ns(rf)*1000); //in case of RF signal, the energy branch contains the period 
              if(energy2) leaf_nrj2.push_back(0);
              if(energy3) leaf_nrj3.push_back(0);
              if(energy4) leaf_nrj4.push_back(0);
              leaf_t.push_back((time_type)clock_g);
              leaf_pu.push_back(false);
              leaf_mult++;           
          }

          //TRef (== QDC0)
          if(myID[faster_data_label(group_data)] == 8){
              faster_data_load (group_data, &TREF_tdc);
              clock_g = faster_data_hr_clock_ns(group_data)*1000;
              leaf_pu.push_back(false);
              leaf_label.push_back(faster_data_label (group_data));
              leaf_nrj.push_back(0);
              if(energy2) leaf_nrj2.push_back(0);
              if(energy3) leaf_nrj3.push_back(0);
              if(energy4) leaf_nrj4.push_back(0);
              leaf_t.push_back((time_type)clock_g);
              leaf_mult++;
          }

          //trapezoidal
          else if ( myID[faster_data_label(group_data)] == 5){
              faster_data_load (group_data, &trapez_adc);
              // Getting absolute time for the hit
              clock_g = faster_data_hr_clock_ns(group_data)*1000;
              if (trapez_adc.pileup!=1) leaf_pu.push_back(false);
              else leaf_pu.push_back(true);
              leaf_label.push_back(faster_data_label (group_data));
              leaf_nrj.push_back(trapez_adc.measure);
              if(energy2) leaf_nrj2.push_back(0);
              if(energy3) leaf_nrj3.push_back(0);
              if(energy4) leaf_nrj4.push_back(0);
              leaf_t.push_back((time_type)clock_g);
              leaf_mult++;
          }

          //QDC1
          else if ( myID[faster_data_label(group_data)] == 1){
              faster_data_load (group_data, &QDC1_qdc);
              clock_g = faster_data_hr_clock_ns(group_data)*1000;
              if (QDC1_qdc.q1_saturated!=1) leaf_pu.push_back(false);
              else leaf_pu.push_back(true);
              leaf_label.push_back(faster_data_label (group_data));
              leaf_nrj.push_back(QDC1_qdc.q1);
              if(energy2) leaf_nrj2.push_back(0);
              if(energy3) leaf_nrj3.push_back(0);
              if(energy4) leaf_nrj4.push_back(0);
              leaf_t.push_back((time_type)clock_g);
              leaf_mult++;
          }

          //QDC2
          else if ( myID[faster_data_label (group_data)] == 2){
              faster_data_load (group_data, &QDC2_qdc);
              clock_g = faster_data_hr_clock_ns(group_data)*1000;
              if (QDC2_qdc.q1_saturated!=1) leaf_pu.push_back(false);
              else leaf_pu.push_back(true);
              leaf_label.push_back(faster_data_label (group_data));
              leaf_nrj.push_back(QDC2_qdc.q1);
              leaf_nrj2.push_back(QDC2_qdc.q2);
              if(energy3) leaf_nrj3.push_back(0);
              if(energy4) leaf_nrj4.push_back(0);
              leaf_t.push_back((time_type)clock_g);
              leaf_mult++;
          }
          //QDC3
          else if ( myID[faster_data_label (group_data)] == 3){
              faster_data_load (group_data, &QDC3_qdc);
              clock_g = faster_data_hr_clock_ns(group_data)*1000;
              if (QDC3_qdc.q1_saturated!=1) leaf_pu.push_back(false);
              else leaf_pu.push_back(true); // reject saturation
              leaf_label.push_back(faster_data_label (group_data));
              leaf_nrj.push_back(QDC3_qdc.q1);
              leaf_nrj2.push_back(QDC3_qdc.q2);
              leaf_nrj3.push_back(QDC3_qdc.q3);
              if(energy4) leaf_nrj4.push_back(0);
              leaf_t.push_back((time_type)clock_g);
              leaf_mult++;
          }

          //QDC4
          else if ( myID[faster_data_label (group_data)] == 4){
              faster_data_load (group_data, &QDC4_qdc);
              clock_g = faster_data_hr_clock_ns(group_data)*1000;
              if (QDC4_qdc.q1_saturated!=1) leaf_pu.push_back(false);
              else leaf_pu.push_back(true); // reject saturation
              leaf_label.push_back(faster_data_label (group_data));
              leaf_nrj.push_back(QDC4_qdc.q1);
              leaf_nrj2.push_back(QDC4_qdc.q2);
              leaf_nrj3.push_back(QDC4_qdc.q3);
              leaf_nrj4.push_back(QDC4_qdc.q4);
              leaf_t.push_back((time_type)clock_g);
              leaf_mult++;
          }

          //CRRC4
          else if ( myID[faster_data_label (group_data)] == 6){
              faster_data_load (group_data, &crrc4_adc);
              clock_g = faster_data_hr_clock_ns(group_data)*1000;
              leaf_pu.push_back(false);
              leaf_label.push_back(faster_data_label (group_data));
              leaf_nrj.push_back(crrc4_adc.measure);
              if(energy2) leaf_nrj2.push_back(0);
              if(energy3) leaf_nrj3.push_back(0);
              if(energy4) leaf_nrj4.push_back(0);
              leaf_t.push_back((time_type)clock_g);
              leaf_mult++;
          }

      faster_buffer_reader_close (group_reader);
      if(keepGroups && leaf_mult > 0) root_tree->Fill(); //if we keep groups, we save the vectors in the tree
      else if (!keepGroups){
        for(int ev = 0; ev < leaf_mult; ev++) //else, we have to loop on all the events
        {
          leaf_label_2 = leaf_label.at(ev);
          leaf_nrj_2 = leaf_nrj.at(ev);
          if(energy2) leaf_nrj2_2 = leaf_nrj2.at(ev);
          if(energy3) leaf_nrj3_2 = leaf_nrj3.at(ev);
          if(energy4) leaf_nrj4_2 = leaf_nrj4.at(ev);
          leaf_t_2 = leaf_t.at(ev);
          leaf_pu_2= leaf_pu.at(ev);
          root_tree->Fill();
        }
      }
    }
    if(alias != GROUP_TYPE_ALIAS){
      leaf_nrj.clear(); //energy
      if(energy2) leaf_nrj2.clear(); //used if QDC2
      if(energy3) leaf_nrj3.clear(); //used if QDC3
      if(energy4) leaf_nrj4.clear(); //used if QDC4
      leaf_t.clear(); //time
      leaf_label.clear(); //label
      leaf_pu.clear(); //pile-up detection
      leaf_mult = 0;
      //RF
      if( myID[faster_data_label(data)] == 7){
          faster_data_load (data, &rf);
          clock_g = faster_data_hr_clock_ns(data)*1000;  //Time in ps
          leaf_label.push_back(faster_data_label (data));
          leaf_nrj.push_back(rf_period_ns(rf)*1000); //in case of RF signal, the energy branch contains the period 
          if(energy2) leaf_nrj2.push_back(0);
          if(energy3) leaf_nrj3.push_back(0);
          if(energy4) leaf_nrj4.push_back(0);
          leaf_t.push_back((time_type)clock_g);
          leaf_pu.push_back(false);
          leaf_mult++;           
      }

      //TRef (== QDC0)
      if(myID[faster_data_label(data)] == 8){
          faster_data_load (data, &TREF_tdc);
          clock_g = faster_data_hr_clock_ns(data)*1000;
          leaf_pu.push_back(false);
          leaf_label.push_back(faster_data_label (data));
          leaf_nrj.push_back(0);
          if(energy2) leaf_nrj2.push_back(0);
          if(energy3) leaf_nrj3.push_back(0);
          if(energy4) leaf_nrj4.push_back(0);
          leaf_t.push_back((time_type)clock_g);
          leaf_mult++;
      }

      //trapezoidal
      else if ( myID[faster_data_label(data)] == 5){
          faster_data_load (data, &trapez_adc);
          // Getting absolute time for the hit
          clock_g = faster_data_hr_clock_ns(data)*1000;
          if (trapez_adc.pileup!=1) leaf_pu.push_back(false);
          else leaf_pu.push_back(true);
          leaf_label.push_back(faster_data_label (data));
          leaf_nrj.push_back(trapez_adc.measure);
          if(energy2) leaf_nrj2.push_back(0);
          if(energy3) leaf_nrj3.push_back(0);
          if(energy4) leaf_nrj4.push_back(0);
          leaf_t.push_back((time_type)clock_g);
          leaf_mult++;
      }

      //QDC1
      else if ( myID[faster_data_label(data)] == 1){
          faster_data_load (data, &QDC1_qdc);
          clock_g = faster_data_hr_clock_ns(data)*1000;
          if (QDC1_qdc.q1_saturated!=1) leaf_pu.push_back(false);
          else leaf_pu.push_back(true);
          leaf_label.push_back(faster_data_label (data));
          leaf_nrj.push_back(QDC1_qdc.q1);
          if(energy2) leaf_nrj2.push_back(0);
          if(energy3) leaf_nrj3.push_back(0);
          if(energy4) leaf_nrj4.push_back(0);
          leaf_t.push_back((time_type)clock_g);
          leaf_mult++;
      }

      //QDC2
      else if ( myID[faster_data_label (data)] == 2){
          faster_data_load (data, &QDC2_qdc);
          clock_g = faster_data_hr_clock_ns(data)*1000;
          if (QDC2_qdc.q1_saturated!=1) leaf_pu.push_back(false);
          else leaf_pu.push_back(true);
          leaf_label.push_back(faster_data_label (data));
          leaf_nrj.push_back(QDC2_qdc.q1);
          leaf_nrj2.push_back(QDC2_qdc.q2);
          if(energy3) leaf_nrj3.push_back(0);
          if(energy4) leaf_nrj4.push_back(0);
          leaf_t.push_back((time_type)clock_g);
          leaf_mult++;
      }
      //QDC3
      else if ( myID[faster_data_label (data)] == 3){
          faster_data_load (data, &QDC3_qdc);
          clock_g = faster_data_hr_clock_ns(data)*1000;
          if (QDC3_qdc.q1_saturated!=1) leaf_pu.push_back(false);
          else leaf_pu.push_back(true); // reject saturation
          leaf_label.push_back(faster_data_label (data));
          leaf_nrj.push_back(QDC3_qdc.q1);
          leaf_nrj2.push_back(QDC3_qdc.q2);
          leaf_nrj3.push_back(QDC3_qdc.q3);
          if(energy4) leaf_nrj4.push_back(0);
          leaf_t.push_back((time_type)clock_g);
          leaf_mult++;
      }

      //QDC4
      else if ( myID[faster_data_label (data)] == 4){
          faster_data_load (data, &QDC4_qdc);
          clock_g = faster_data_hr_clock_ns(data)*1000;
          if (QDC4_qdc.q1_saturated!=1) leaf_pu.push_back(false);
          else leaf_pu.push_back(true); // reject saturation
          leaf_label.push_back(faster_data_label (data));
          leaf_nrj.push_back(QDC4_qdc.q1);
          leaf_nrj2.push_back(QDC4_qdc.q2);
          leaf_nrj3.push_back(QDC4_qdc.q3);
          leaf_nrj4.push_back(QDC4_qdc.q4);
          leaf_t.push_back((time_type)clock_g);
          leaf_mult++;
      }

      //CRRC4
      else if ( myID[faster_data_label (data)] == 6){
          faster_data_load (data, &crrc4_adc);
          clock_g = faster_data_hr_clock_ns(data)*1000;
          leaf_pu.push_back(false);
          leaf_label.push_back(faster_data_label (data));
          leaf_nrj.push_back(crrc4_adc.measure);
          if(energy2) leaf_nrj2.push_back(0);
          if(energy3) leaf_nrj3.push_back(0);
          if(energy4) leaf_nrj4.push_back(0);
          leaf_t.push_back((time_type)clock_g);
          leaf_mult++;
      }

      if(leaf_mult > 0){
        if(keepGroups) root_tree->Fill();
        else if (!keepGroups){
          for(int ev = 0; ev < leaf_mult; ev++){
            leaf_label_2 = leaf_label.at(ev);
            leaf_nrj_2 = leaf_nrj.at(ev);
            if(energy2) leaf_nrj2_2 = leaf_nrj2.at(ev);
            if(energy3) leaf_nrj3_2 = leaf_nrj3.at(ev);
            if(energy4) leaf_nrj4_2 = leaf_nrj4.at(ev);
            leaf_t_2 = leaf_t.at(ev);
            leaf_pu_2= leaf_pu.at(ev);
            root_tree->Fill();
          }
        }
      }
    }
  }
  //We save the output file
  faster_file_reader_close (reader);
  root_file->cd();
  root_tree->AutoSave();
  root_file->Close();
}