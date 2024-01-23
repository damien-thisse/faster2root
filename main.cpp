#include "F2R_Parameters.h"
#include "F2R_Core.h"
#include<iostream>

#include"TROOT.h"
#include<fstream>
#include<thread>
#include<pthread.h>

using namespace std;

int main(int argc, char** argv){
    setenv("ROOT_ENABLE_IMT", "1", 1);
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
    vector<pthread_t> threads;
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
    pthread_mutex_t myInputMutex;

    //Start of the conversion in multi-thread mode
    //Each time a file is converted, a new one is given to the thread
    for(int k = 0; k < nOfThreads; k++){
        threads.emplace_back([&input_file, &myInputMutex, &endOfFile](){
            F2R_Parameters& parameters = F2R_Parameters::getInstance();
            while(!endOfFile)
            {
                string filename = "";
                pthread_mutex_lock(&myInputMutex);
                input_file >> filename;
                if(filename == ""){
                    endOfFile = true;
                    pthread_mutex_unlock(&myInputMutex);
                    break;
                }
                cout << "Thread #" << std::this_thread::get_id() << ": start processing " << filename << endl;
                pthread_mutex_unlock(&myInputMutex);

                Convert(filename);
                if(!parameters.getKeepGroups()){
                    Sort(filename);
                }
            }
        });
    }
    for (int i = 0; i < nOfThreads; i++) {
        pthread_join(threads[i], NULL);
    }
    return 0;
}