#ifndef F2R_PARAMETERS
#define F2R_PARAMETERS

#include<map>
#include<string>

class F2R_Parameters
{
public:
    // Fonction statique pour obtenir l'instance unique de la classe
    static F2R_Parameters& getInstance() {
        // Vérifier si l'instance existe déjà
        // Si non, en créer une nouvelle
        if (instance == nullptr) {
            instance = new F2R_Parameters();
        }
        // Retourner l'instance unique
        return *instance;
    }

    // Autres méthodes de la classe
    bool setIDTable();

    int getMaxEnergyBranchNeeded() const {return F2R_MaxEnergyBranchNeeded;}
    std::string getOutputDirectory() const {return F2R_OutputDirectory;}
    bool getKeepGroups() const {return F2R_KeepGroups;}
    std::map<int,int> getIDTable() const {return F2R_IDTable;}
    
private:
    // Constructeur privé pour empêcher l'instanciation directe
    F2R_Parameters();

    // Destructeur privé pour éviter la suppression accidentelle
    ~F2R_Parameters();

    // Empêcher la copie de l'objet
    F2R_Parameters(const F2R_Parameters&) = delete;
    F2R_Parameters& operator=(const F2R_Parameters&) = delete;

    // Instance unique de la classe
    static F2R_Parameters* instance;

    //Variables membres de la classe
    std::map<int,int> F2R_IDTable;
    std::string F2R_ListOfFiles;
    std::string F2R_OutputDirectory;
    std::string F2R_IDTableFile;
    int F2R_NumberOfThreads;
    bool F2R_KeepGroups;

    int F2R_MaxEnergyBranchNeeded;
};

#endif //F2R_PARAMETERS