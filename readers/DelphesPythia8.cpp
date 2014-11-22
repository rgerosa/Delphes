#include <stdexcept>
#include <iostream>
#include <sstream>
#include <vector>

#include <signal.h>

#include "Pythia8/Pythia.h"

#include "TROOT.h"
#include "TApplication.h"

#include "TFile.h"
#include "TObjArray.h"
#include "TStopwatch.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TLorentzVector.h"

#include "modules/Delphes.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"

#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootProgressBar.h"

using namespace std;

static bool interrupted = false;

void SignalHandler(int sig){
  interrupted = true;
}

//**********************************************************
//*** Method to filter LHE event on the fly and store info
//**********************************************************
bool lhe_event_preselection(vector< vector<float> >* LHE_event, float Mjj_cut, int SkimFullyHadronic, vector<ExRootTreeBranch*> branchVector);

//********************************************************************
//*** Input Converter for Delphes (frpm Pythia8 to delphes particles)
//********************************************************************

void ConvertInput(Long64_t eventCounter, Pythia8::Pythia* pythia,
		  ExRootTreeBranch *branch, DelphesFactory *factory,
		  TObjArray *allParticleOutputArray, TObjArray *stableParticleOutputArray, TObjArray *partonOutputArray,
		  TStopwatch *readStopWatch, TStopwatch *procStopWatch);
//*******************************************************************************************************************
// main code : ./DelphesPythia8 <delphes_card> <lhe file> <output root file> <mjj cut> <filter fully hadronic FS> <starting event> <total event>
//*******************************************************************************************************************

int main(int argc, char *argv[]){

    
  // Initialize objects
  char appName[] = "DelphesPythia8";
  std::stringstream message;
  std::string inputFile;

  TFile *outputFile = 0;

  TStopwatch readStopWatch, procStopWatch;

  ExRootTreeWriter *treeWriter  = 0;
  ExRootTreeWriter *treeHepMC   = 0;    
  ExRootTreeBranch *branchEvent = 0;
  std::vector<ExRootTreeBranch*> branchGen;
  ExRootConfReader *confReader = 0;
  Delphes *modularDelphes = 0;
  DelphesFactory *factory = 0;

  TObjArray *stableParticleOutputArray = 0, *allParticleOutputArray = 0, *partonOutputArray = 0;

  Long64_t eventCounter, errorCounter;

  Pythia8::Pythia *pythia = 0;

  if(argc < 4){
    std::cout << "------------------------------------Manual----------------------------------------" <<std::endl;
    std::cout << "- Usage: " << appName << " config_file" << " lhe_file" << " output_file" << " Mjj_cut" << " start" << " number" << " signal" <<std::endl;
    std::cout << "- config_file          ->  configuration file in Tcl format" << std::endl;
    std::cout << "- input_file           ->  lhe file for Pythia8" << std::endl;
    std::cout << "- output_file          ->  output file in ROOT format" << std::endl;
    std::cout << "- Mjj_cut  (optional)  ->  cut on Mjj in GeV -- default = 0 GeV" << std::endl;
    std::cout << "- filter   (optional)  ->  flag to filter fully hadronic events at LHE level -- default = 1" << std::endl;
    std::cout << "- start    (optional)  ->  number of starting event" << std::endl;
    std::cout << "- number   (optional)  ->  number of total events to be processed" << std::endl;
    std::cout << "- signal   (optional)  ->  1 if the sample is a graviton signal sample" << std::endl;
    std::cout << "----------------------------------------------------------------------------------" << std::endl << std::endl;  
    return 1;
  }

  signal(SIGINT,SignalHandler);
  gROOT->SetBatch();

  int appargc = 1;
  char* appargv [] = {appName};
  TApplication app(appName, &appargc,appargv); // open a TApplication process

  try{

    inputFile  = argv[2]; // input file name for LHE
    outputFile = TFile::Open(argv[3], "RECREATE"); // output file name

    if(outputFile == NULL or outputFile == 0){
      message << "can't create output file " << argv[3];
      throw runtime_error(message.str());
    }

    int fsignal = 0;        
    if( argc >= 9){
      fsignal = atoi(argv[8]);
    }

    //--- create output tree ---
    treeWriter = new ExRootTreeWriter(outputFile, "Delphes");

    //--- create gen (lhe level) branch --- 
    // done by Simone --> std::vector for each particle type  
    // new branches for the latino tree type flat tree is modified (added) by Deniz

    //--- lep number
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_lep_number"));      

    //--- leptons
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_lep_pt1"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_lep_eta1"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_lep_phi1"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_lep_pid1"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_lep_pt2"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_lep_eta2"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_lep_phi2"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_lep_pid2"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_lep_pt3"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_lep_eta3"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_lep_phi3"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_lep_pid3"));
    //--- neutrinos
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_nu_pt1"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_nu_eta1"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_nu_phi1"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_nu_pid1"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_nu_pt2"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_nu_eta2"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_nu_phi2"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_nu_pid2"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_nu_pt3"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_nu_eta3"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_nu_phi3"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_nu_pid3"));

    //--- gen partons infos
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_par_pt1"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_par_eta1"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_par_phi1"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_par_pid1"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_par_pt2"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_par_eta2"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_par_phi2"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_par_pid2"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_par_pt3"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_par_eta3"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_par_phi3"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_par_pid3"));
        
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_par_number"));     
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_nu_number"));     
    //--- gen W infos
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_W_pt"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_W_eta"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_W_phi"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_W_mass"));
    branchGen.push_back(treeWriter->NewFloatBranch("lhe_W_pid"));

       
        
    //--- deals with the HepMc output of Pythia8 ---> no need to store it
    treeHepMC = new ExRootTreeWriter();
    branchEvent = treeHepMC->NewBranch("Event",HepMCEvent::Class());

    //----- Delphes init ----- --> Card reader       
    confReader = new ExRootConfReader;
    confReader->ReadFile(argv[1]);

    modularDelphes = new Delphes("Delphes");
    modularDelphes->SetConfReader(confReader); // read the config file
    modularDelphes->SetTreeWriter(treeWriter); // set the output tree

    factory = modularDelphes->GetFactory();
    allParticleOutputArray    = modularDelphes->ExportArray("allParticles"); 
    stableParticleOutputArray = modularDelphes->ExportArray("stableParticles");
    partonOutputArray         = modularDelphes->ExportArray("partons");

    modularDelphes->InitTask();

    //-----------------------------
    //----- Initialize pythia -----
    //-----------------------------

    pythia = new Pythia8::Pythia;        

    if(pythia == NULL or pythia == 0){
      throw runtime_error("can't create Pythia instance");
    }


    // Mjj cut set to zero as default, starting event and number of events
    std::string sSeed = "0";
    float Mjj_cut     = 0;
    int   skimFullyHadronic  = 1;
    if (argc >= 5){
      Mjj_cut = atof(argv[4]);
    }

    if (argc >= 6){
      skimFullyHadronic = atoi(argv[5]);
    }
                
    int startEvent = 0, nEvent = -1;
    if (argc >= 7) {
      startEvent = atoi(argv[6]);
    }

    if (argc >= 8) {
      nEvent = atoi(argv[7]);                                                     
    }

    //--- Initialize Les Houches Event File run. List initialization information.
    std::string sfile       = "Beams:LHEF ="+inputFile;
    std::string sRandomSeed = "Random:seed = "+sSeed;
    //--- random seed from start event number
    pythia->readString("Random:setSeed = on");
    pythia->readString("HadronLevel:Hadronize = on"); // turn on the hadronize module
    pythia->readString(sRandomSeed.c_str());          // random seed set
    pythia->readString("Beams:frameType = 4");        
    pythia->readString(sfile.c_str());
    pythia->init();

    //------------------------------------------------------------
    //----- initialize fast lhe file reader for preselection -----
    //------------------------------------------------------------

    ifstream inputLHE (inputFile.c_str(), ios::in); // read the input LHE file
    int count          =-1;
    int skippedCounter = 0;
    char buffer[256];
    std::vector<std::vector<float> > LHE_event;

    //--- start from specified event
    while(count < startEvent){
      inputLHE.getline(buffer,256);
      while(strcmp(buffer, "<event>") != 0){
	inputLHE.getline(buffer,258);
      }
      count++;
    }

    if(pythia->LHAeventSkip(startEvent)){
      std::cout << "### skipped first " << startEvent << " events" << std::endl;
    }

    ExRootProgressBar progressBar(-1);

    // Loop over all events
    errorCounter = 0;
    eventCounter = 0;
    modularDelphes->Clear();
    readStopWatch.Start();
    while(!interrupted){
      if(eventCounter >= nEvent && nEvent != -1){
	break;
      }
      LHE_event.clear();  

      //--- end of file check
      if(strcmp(buffer, "</LesHouchesEvents>") == 0)        {
	std::cout << "### end of LHE file reached! ### " << std::endl;
	break;
      }
      //--- reads and store the event
      else if(strcmp(buffer, "<event>") == 0){   

	int nPart;
	float info_tmp;
	inputLHE >> nPart;
	for(int i=0; i<5; i++)
	  inputLHE >> info_tmp;		
	for(int i=0; i<nPart; i++){
		
	  std::vector<float> LHE_particle;
	  LHE_particle.clear();                
	  for(int j=0; j<13; j++){
	    inputLHE >> info_tmp;
	    LHE_particle.push_back(info_tmp);
	  }
	  LHE_event.push_back(LHE_particle);
	}
	//--- process only selected events
	if(lhe_event_preselection(&LHE_event, Mjj_cut, skimFullyHadronic, branchGen)){
	  if(!pythia->next()){
	    //--- If failure because reached end of file then exit event loop
	    if (pythia->info.atEndOfFile()){
	      std::cerr << "Aborted since reached end of Les Houches Event File" << std::endl;
	      break;
	    }
	    //--- keep trace of faulty events
	    errorCounter++;
	  }
	  readStopWatch.Stop();

	  //--- delphes simulation fase
	  procStopWatch.Start();
	  ConvertInput(eventCounter,pythia,branchEvent,factory,allParticleOutputArray,stableParticleOutputArray,partonOutputArray,&readStopWatch,&procStopWatch);
	  modularDelphes->ProcessTask();
	  procStopWatch.Stop();
		   
	  //--- filling the output tree
	  treeWriter->Fill();
    
	  //--- logistic 
	  treeWriter->Clear();
	  modularDelphes->Clear();
	  readStopWatch.Start();
	}
	else{		
	  if(pythia->LHAeventSkip(1)) skippedCounter++;
	  else std::cout << "### ERROR: couldn't skip event" << endl;
	}
		
	eventCounter++;
	progressBar.Update(eventCounter, eventCounter);
      }
      inputLHE.getline(buffer,256);
    }	

    progressBar.Update(eventCounter, eventCounter, kTRUE);
    progressBar.Finish();

    std::cout << "--------------------Statistics---------------------" <<std::endl;
    std::cout << "-#######  read events:        " << eventCounter << std::endl; 
    std::cout << "-#######  failed events:      " << errorCounter << std::endl;
    std::cout << "-#######  skipped events:     " << skippedCounter << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
        
    modularDelphes->FinishTask();
    treeWriter->Write();
    
    std::cout << std::endl <<  "** Exiting..." << std::endl;

    delete pythia;
    delete confReader;
	    
    return 0;
  }
  catch(runtime_error &e){
    if(treeWriter) delete treeWriter;
    if(outputFile) delete outputFile;
    std::cerr << "** ERROR: " << e.what() << std::endl;
    return 1;
  }

  return 0 ;

}


//*******************************************************

bool lhe_event_preselection(vector< vector<float> >* LHE_event, float Mjj_cut, int SkimFullyHadronic, vector<ExRootTreeBranch*> branchVector){

  std::vector<TLorentzVector> outPartons;
  int leptons   = 0;
  int Mjj_check = 0;
    
  //---search for final state partons in the event---
  for(size_t iPart = 0; iPart < LHE_event->size(); iPart++){

    std::vector<float> particle = LHE_event->at(iPart);
    TLorentzVector tmp4vect;
    tmp4vect.SetPxPyPzE(particle.at(6), particle.at(7), particle.at(8), particle.at(9));
    if(particle.at(1) == 1 && ((abs(particle.at(0)) > 0 && abs(particle.at(0)) < 7) || abs(particle.at(0)) == 21)){
      outPartons.push_back(tmp4vect);
    }
    if(particle.at(1) == 1 && (abs(particle.at(0)) == 11 || abs(particle.at(0)) == 13 || abs(particle.at(0)) == 15)){
      leptons++;
    }
  }

  //---reject fully hadronic events --> fully hadronic events are filtered out by simone
  if(SkimFullyHadronic and leptons < 1){
    return false;
  }

  //---apply the VBF preselection cut on Mjj---
  for(size_t iPartons = 0; iPartons < outPartons.size(); iPartons++){
    TLorentzVector tmp4vect = outPartons.at(iPartons);
    for(size_t jQuark = iPartons+1; jQuark < outPartons.size(); jQuark++){
      tmp4vect += outPartons.at(jQuark);
      if(tmp4vect.M() > Mjj_cut){
	Mjj_check = 1;
      }
    }
  }

  if( outPartons.size() > 2 && Mjj_check == 0 ){
    return false;
  }

  //---if preselection is passed save lhe infos---
  //---changed by deniz
    
  int nlepton = 0, nparton = 0, nneutrino = 0 ;
	
  vector<float> *n_lep_gen, *n_par_gen, *n_nu_gen;
  vector<float> *lep_pt1, *lep_eta1, *lep_phi1, *lep_flv1;
  vector<float> *lep_pt2, *lep_eta2, *lep_phi2, *lep_flv2;
  vector<float> *lep_pt3, *lep_eta3, *lep_phi3, *lep_flv3;

  vector<float> *nu_pt1, *nu_eta1, *nu_phi1, *nu_flv1;
  vector<float> *nu_pt2, *nu_eta2, *nu_phi2, *nu_flv2;
  vector<float> *nu_pt3, *nu_eta3, *nu_phi3, *nu_flv3;
    
  vector<float> *p_pt1, *p_eta1, *p_phi1, *p_flv1;
  vector<float> *p_pt2, *p_eta2, *p_phi2, *p_flv2;
  vector<float> *p_pt3, *p_eta3, *p_phi3, *p_flv3;
    

    
  vector<float> *W_pt, *W_phi, *W_eta, *W_m, *W_pid;


  n_lep_gen = (std::vector<float>*)((branchVector.at(0))->NewFloatEntry());  
  lep_pt1 = (vector<float>*)((branchVector.at(1))->NewFloatEntry());
  lep_eta1 = (vector<float>*)((branchVector.at(2))->NewFloatEntry());
  lep_phi1 = (vector<float>*)((branchVector.at(3))->NewFloatEntry());
  lep_flv1 = (vector<float>*)((branchVector.at(4))->NewFloatEntry());    
  lep_pt2 = (vector<float>*)((branchVector.at(5))->NewFloatEntry());
  lep_eta2 = (vector<float>*)((branchVector.at(6))->NewFloatEntry());
  lep_phi2 = (vector<float>*)((branchVector.at(7))->NewFloatEntry());
  lep_flv2 = (vector<float>*)((branchVector.at(8))->NewFloatEntry());  
  lep_pt3 = (vector<float>*)((branchVector.at(9))->NewFloatEntry());
  lep_eta3 = (vector<float>*)((branchVector.at(10))->NewFloatEntry());
  lep_phi3 = (vector<float>*)((branchVector.at(11))->NewFloatEntry());
  lep_flv3 = (vector<float>*)((branchVector.at(12))->NewFloatEntry());  
  nu_pt1 = (vector<float>*)((branchVector.at(13))->NewFloatEntry());
  nu_eta1 = (vector<float>*)((branchVector.at(14))->NewFloatEntry());
  nu_phi1 = (vector<float>*)((branchVector.at(15))->NewFloatEntry());
  nu_flv1 = (vector<float>*)((branchVector.at(16))->NewFloatEntry());
  nu_pt2 = (vector<float>*)((branchVector.at(17))->NewFloatEntry());
  nu_eta2 = (vector<float>*)((branchVector.at(18))->NewFloatEntry());
  nu_phi2 = (vector<float>*)((branchVector.at(19))->NewFloatEntry());
  nu_flv2 = (vector<float>*)((branchVector.at(20))->NewFloatEntry());
  nu_pt3 = (vector<float>*)((branchVector.at(21))->NewFloatEntry());
  nu_eta3 = (vector<float>*)((branchVector.at(22))->NewFloatEntry());
  nu_phi3 = (vector<float>*)((branchVector.at(23))->NewFloatEntry());
  nu_flv3 = (vector<float>*)((branchVector.at(24))->NewFloatEntry());
  p_pt1 = (vector<float>*)((branchVector.at(25))->NewFloatEntry());
  p_eta1 = (vector<float>*)((branchVector.at(26))->NewFloatEntry());
  p_phi1 = (vector<float>*)((branchVector.at(27))->NewFloatEntry());
  p_flv1 = (vector<float>*)((branchVector.at(28))->NewFloatEntry());
  p_pt2 = (vector<float>*)((branchVector.at(29))->NewFloatEntry());
  p_eta2 = (vector<float>*)((branchVector.at(30))->NewFloatEntry());
  p_phi2 = (vector<float>*)((branchVector.at(31))->NewFloatEntry());
  p_flv2 = (vector<float>*)((branchVector.at(32))->NewFloatEntry());
  p_pt3 = (vector<float>*)((branchVector.at(33))->NewFloatEntry());
  p_eta3 = (vector<float>*)((branchVector.at(34))->NewFloatEntry());
  p_phi3 = (vector<float>*)((branchVector.at(35))->NewFloatEntry());
  p_flv3 = (vector<float>*)((branchVector.at(36))->NewFloatEntry());
  n_par_gen = (vector<float>*)((branchVector.at(37))->NewFloatEntry());
  n_nu_gen = (vector<float>*)((branchVector.at(38))->NewFloatEntry());
  W_pt      = (vector<float>*)((branchVector.at(39))->NewFloatEntry());
  W_eta     = (vector<float>*)((branchVector.at(40))->NewFloatEntry());
  W_phi     = (vector<float>*)((branchVector.at(41))->NewFloatEntry());
  W_m       = (vector<float>*)((branchVector.at(42))->NewFloatEntry());
  W_pid     = (vector<float>*)((branchVector.at(43))->NewFloatEntry());
   
  //---loop on lhe events particle searching for W's---
  for(size_t iPart = 0; iPart < LHE_event->size(); iPart++){ 
    std::vector<float> particle = LHE_event->at(iPart);
    TLorentzVector tmp4vect;
    tmp4vect.SetPxPyPzE(particle.at(6), particle.at(7), particle.at(8), particle.at(9));   
    if(particle.at(1) == 2 && abs(particle.at(0)) == 24){ 
      W_pt->push_back(tmp4vect.Pt());
      W_eta->push_back(tmp4vect.Eta());
      W_phi->push_back(tmp4vect.Phi());
      W_m->push_back(tmp4vect.M());
      W_pid->push_back(particle.at(0));
    }

    //---Partons
    if(particle.at(1) == 1 && ( (particle.at(0) > -7 && particle.at(0) < 7) || abs(particle.at(0)) == 21))
      {
	nparton++;
	if(nparton==1){
	  p_pt1->push_back(tmp4vect.Pt());
	  p_eta1->push_back(tmp4vect.Eta());
	  p_phi1->push_back(tmp4vect.Phi());
	  p_flv1->push_back(particle.at(0));
	}
	if(nparton==2){
	  p_pt2->push_back(tmp4vect.Pt());
	  p_eta2->push_back(tmp4vect.Eta());
	  p_phi2->push_back(tmp4vect.Phi());
	  p_flv2->push_back(particle.at(0));
	}
	if(nparton==3){
	  p_pt3->push_back(tmp4vect.Pt());
	  p_eta3->push_back(tmp4vect.Eta());
	  p_phi3->push_back(tmp4vect.Phi());
	  p_flv3->push_back(particle.at(0));
	}
         
      }

    //---charged leptons
    if(particle.at(1) == 1 && (particle.at(0) == 11 || particle.at(0) == -11 || particle.at(0) == 13 || particle.at(0) == -13 || particle.at(0) == 15 || particle.at(0) == -15)){
      nlepton++;
      if(nlepton == 1){
	lep_pt1->push_back(tmp4vect.Pt());
	lep_eta1->push_back(tmp4vect.Eta());
	lep_phi1->push_back(tmp4vect.Phi());
	lep_flv1->push_back(particle.at(0));
      }
      if(nlepton == 2){
	lep_pt2->push_back(tmp4vect.Pt());
	lep_eta2->push_back(tmp4vect.Eta());
	lep_phi2->push_back(tmp4vect.Phi());
	lep_flv2->push_back(particle.at(0));
				
      }
      if(nlepton == 3){
	lep_pt3->push_back(tmp4vect.Pt());
	lep_eta3->push_back(tmp4vect.Eta());
	lep_phi3->push_back(tmp4vect.Phi());
	lep_flv3->push_back(particle.at(0));
				
      }
    }

    //---neutrinos
    if(particle.at(1) == 1 && (particle.at(0) == 12 || particle.at(0) == -12 || particle.at(0) == 14 || particle.at(0) == -14 || particle.at(0) == 16 || particle.at(0) == -16)){
      nneutrino++;
      if(nneutrino ==1){
	nu_pt1->push_back(tmp4vect.Pt());
	nu_eta1->push_back(tmp4vect.Eta());
	nu_phi1->push_back(tmp4vect.Phi());
	nu_flv1->push_back(particle.at(0));
      }
      if(nneutrino ==2){
	nu_pt2->push_back(tmp4vect.Pt());
	nu_eta2->push_back(tmp4vect.Eta());
	nu_phi2->push_back(tmp4vect.Phi());
	nu_flv2->push_back(particle.at(0));
      }
      if(nneutrino ==3){
	nu_pt3->push_back(tmp4vect.Pt());
	nu_eta3->push_back(tmp4vect.Eta());
	nu_phi3->push_back(tmp4vect.Phi());
	nu_flv3->push_back(particle.at(0));
      }
    }
        
  }

  n_lep_gen->push_back(nlepton);
  n_par_gen->push_back(nparton);
  n_nu_gen->push_back(nneutrino);

  return true;
  
}

//***************************************************

void ConvertInput(Long64_t eventCounter, Pythia8::Pythia* pythia,
		  ExRootTreeBranch *branch, DelphesFactory *factory,
		  TObjArray *allParticleOutputArray, TObjArray *stableParticleOutputArray, TObjArray *partonOutputArray,
		  TStopwatch *readStopWatch, TStopwatch *procStopWatch){

  HepMCEvent *element  = 0;
  Candidate *candidate = 0;
  TDatabasePDG *pdg = 0; 
  TParticlePDG *pdgParticle = 0;
  Int_t pdgCode;

  Int_t pid, status;
  Double_t px, py, pz, e, mass;
  Double_t x, y, z, t;

  // event information
  element = static_cast<HepMCEvent*>(branch->NewEntry());

  element->Number = eventCounter;

  element->ProcessID = pythia->info.code();
  element->MPI       = 1;
  element->Weight    = pythia->info.weight();
  element->Scale     = pythia->info.QRen();
  element->AlphaQED  = pythia->info.alphaEM();
  element->AlphaQCD  = pythia->info.alphaS();

  element->ID1 = pythia->info.id1();
  element->ID2 = pythia->info.id2();
  element->X1  = pythia->info.x1();
  element->X2  = pythia->info.x2();
  element->ScalePDF = pythia->info.QFac();
  element->PDF1 = pythia->info.pdf1();
  element->PDF2 = pythia->info.pdf2();

  element->ReadTime = readStopWatch->RealTime();
  element->ProcTime = procStopWatch->RealTime();

  pdg = TDatabasePDG::Instance();

  for(int i = 0; i < int(pythia->event.size()); ++i){

    Pythia8::Particle &particle = pythia->event[i];
        
    pid    = particle.id();
    status = pythia->event.statusHepMC(i);
    px = particle.px(); 
    py = particle.py(); 
    pz = particle.pz(); 
    e  = particle.e(); 
    mass = particle.m();
    x = particle.xProd(); 
    y = particle.yProd(); 
    z = particle.zProd(); 
    t = particle.tProd();

    candidate = factory->NewCandidate();

    candidate->PID = pid;
    pdgCode = TMath::Abs(candidate->PID);
    candidate->Status = status;

    candidate->M1 = particle.mother1() - 1;
    candidate->M2 = particle.mother2() - 1;

    candidate->D1 = particle.daughter1() - 1;
    candidate->D2 = particle.daughter2() - 1;

    pdgParticle = pdg->GetParticle(pid);
    candidate->Charge = pdgParticle ? Int_t(pdgParticle->Charge()/3.0) : -999;
    candidate->Mass = mass;

    candidate->Momentum.SetPxPyPzE(px, py, pz, e);

    candidate->Position.SetXYZT(x, y, z, t);

    allParticleOutputArray->Add(candidate);

    if(status == 1){ // stable particles are Pythia8 status 1, allParticle are all the pythia8 particles
      stableParticleOutputArray->Add(candidate);
    }
    else if(pdgCode <= 5 || pdgCode == 21 || pdgCode == 15){ // only partons + gluon 
      partonOutputArray->Add(candidate);
    }
  }
}



