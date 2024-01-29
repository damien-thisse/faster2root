#include "F2R_Parameters.h"
#include "F2R_Core.h"
#include<iostream>

#include<fstream>
#include<thread>
#include<mutex>
#include<pthread.h>

using namespace std;

int main(){
    ROOT::EnableThreadSafety();
    //Instance for the first time the parameters handler class
    F2R_Parameters& parameters = F2R_Parameters::getInstance();

    //If the ID file given to the parameter class doesn't exist, stop the code
    if(!parameters.setIDTable()){
        return -1;
    }

    ifstream input_file(parameters.getListOfFiles().c_str());
    int filecounter = 0;
    if(input_file.is_open()){
        while(input_file.good()){
            string oneline;
            getline(input_file, oneline);
            if(oneline.size() > 1){
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
    input_file.open(parameters.getListOfFiles().c_str());
    int nOfThreads = parameters.getNumberOfThread();
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
    bool keepGroups = parameters.getKeepGroups();

    thread_data data_for_MT(parameters.getListOfFiles().c_str(), endOfFile, keepGroups);

    //pthread_t threads_process[nOfThreads];
    vector<thread> threads_process;
    mutex myInputMutex;
    // for(int i = 0; i < nOfThreads; i++){
    //     pthread_create(&threads_process[i], NULL, process, &data_for_MT);
    // }
    // while(!endOfFile)
    // {
    //     string filename = "";
    //     input_file >> filename;
    //     if(filename == ""){
    //         endOfFile = true;
    //         break;
    //     }

    //     cout << "Will Convert : " << filename << endl;

    //     Convert(filename);
    //     if(!keepGroups){
    //         cout << "Will Sort : " << filename << endl;
    //         Sort(filename);
    //     }
    // }

    // Start of the conversion in multi-thread mode
    // Each time a file is converted, a new one is given to the thread
    for(int k = 0; k < nOfThreads; k++){
        threads_process.emplace_back([&input_file, &myInputMutex, &endOfFile, &keepGroups](){
            while(!endOfFile)
            {
                string filename = "";
                myInputMutex.lock();
                input_file >> filename;
                if(filename == ""){
                    endOfFile = true;
                    myInputMutex.unlock();
                    break;
                }
                cout << "Thread #" << std::this_thread::get_id() << ": start processing " << filename << endl;
                myInputMutex.unlock();

                Convert(filename);
                if(!keepGroups){
                    Sort(filename);
                }
            }
        });
    }
    // for (int i = 0; i < nOfThreads; i++) {
    //     pthread_join(threads_process[i], NULL);
    // }

     for (int i = 0; i < nOfThreads; i++) {
        threads_process[i].join();
    }
    return 0;
}