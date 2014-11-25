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
bool lhe_event_preselection(vector< vector<float> >* LHE_event, float Mjj_cut, int SkimFullyHadronic,  DelphesFactory *factory, TObjArray* LHEparticlesArray);

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

    
  // usage of the code
  if(argc < 4){
    std::cout << "------------------------------------Manual----------------------------------------" <<std::endl;
    std::cout << "- Usage: " << "DelphesPythia8 " << " config_file" << " lhe_file" << " output_file" << " Mjj_cut" << " start" << " number" << " signal" <<std::endl;
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
  char* appargv [] = {"DelphesPythia8"};
  TApplication app(appargv[0], &appargc,appargv); // open a TApplication process

  // Initialize objects
  std::string inputFile;
  TFile *outputFile = 0;

  // parsing input parameters
  ExRootTreeWriter *treeWriter = 0;
  try{

    inputFile  = argv[2]; // input file name for LHE
    outputFile = TFile::Open(argv[3], "RECREATE"); // output file name

    if(outputFile == NULL or outputFile == 0){
      std::stringstream message;
      message << "can't create output file " << argv[3];
      throw runtime_error(message.str());
    }

    //--- create output tree ---
    treeWriter = new ExRootTreeWriter(outputFile, "Delphes");

    // Mjj cut set to zero as default, starting event and number of events
    std::string sSeed = "0";
    float       Mjj_cut     = 0;
    int         skimFullyHadronic  = 1;
    int         startEvent = 0, nEvent = -1;

    if (argc >= 5) Mjj_cut = atof(argv[4]);
    
    if (argc >= 6) skimFullyHadronic = atoi(argv[5]);
                    
    if (argc >= 7) startEvent = atoi(argv[6]);
    
    if (argc >= 8) nEvent = atoi(argv[7]);                                                     
    
    //--- deals with the HepMc output of Pythia8 ---> no need to store it
    ExRootTreeWriter *treeHepMC = new ExRootTreeWriter();
    ExRootTreeBranch *branchEvent = treeHepMC->NewBranch("Event",HepMCEvent::Class());

    //----- Delphes init ----- --> Card reader       
    ExRootConfReader *confReader = new ExRootConfReader;
    confReader->ReadFile(argv[1]);

    Delphes *modularDelphes = new Delphes("Delphes");
    modularDelphes->SetConfReader(confReader); // read the config file
    modularDelphes->SetTreeWriter(treeWriter); // set the output tree

    DelphesFactory *factory = modularDelphes->GetFactory();

    TObjArray *allParticleOutputArray    = modularDelphes->ExportArray("allParticles"); 
    TObjArray *stableParticleOutputArray = modularDelphes->ExportArray("stableParticles");
    TObjArray *partonOutputArray         = modularDelphes->ExportArray("partons");
    TObjArray *LHEparticlesArray         = modularDelphes->ExportArray("LHEParticles");

    modularDelphes->InitTask();

    //------------------------------------------------------------
    //----- initialize fast lhe file reader for preselection -----
    //------------------------------------------------------------

    ifstream inputLHE (inputFile.c_str(), ios::in); // read the input LHE file
    int count          =-1;
    int skippedCounter = 0;
    char buffer[256];
    std::vector<std::vector<float> > LHE_event;

    //--- start from specified event --> read the LHE file
    while(count < startEvent){
      inputLHE.getline(buffer,256);
      while(strcmp(buffer, "<event>") != 0){
	inputLHE.getline(buffer,258);
      }
      count++;
    }

    //-----------------------------
    //----- Initialize pythia -----
    //-----------------------------
    TStopwatch readStopWatch, procStopWatch;
    Long64_t eventCounter, errorCounter;

    Pythia8::Pythia *pythia = new Pythia8::Pythia;        

    if(pythia == NULL or pythia == 0){
      throw runtime_error("can't create Pythia instance");
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
      if(eventCounter >= nEvent && nEvent != -1) break;      
      LHE_event.clear();  

      //--- end of file check
      if(strcmp(buffer, "</LesHouchesEvents>") == 0) {
	std::cout << "### end of LHE file reached! ### " << std::endl;
	break;
      }

      //--- reads and store the event
      else if(strcmp(buffer, "<event>") == 0){   
	int nPart;
	float info_tmp;
	inputLHE >> nPart;
	for(int i=0; i<5; i++) inputLHE >> info_tmp;		
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
	if(lhe_event_preselection(&LHE_event,Mjj_cut,skimFullyHadronic,factory,LHEparticlesArray)){
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

bool lhe_event_preselection(vector< vector<float> >* LHE_event, float Mjj_cut, int SkimFullyHadronic,  DelphesFactory *factory, TObjArray* LHEparticlesArray){

  std::vector<TLorentzVector> outPartons;
  int leptons   = 0;
  int Mjj_check = 0;

  TParticlePDG *pdgParticle = 0;
  TDatabasePDG *pdg = 0; 
  Candidate *candidate = 0;

  pdg = TDatabasePDG::Instance();

  //---search for final state partons in the event---
  for(size_t iPart = 0; iPart < LHE_event->size(); iPart++){

    std::vector<float> particle = LHE_event->at(iPart);
    TLorentzVector tmp4vect;
    tmp4vect.SetPxPyPzE(particle.at(6), particle.at(7), particle.at(8), particle.at(9));

    // store output parton information
    if(particle.at(1) == 1 && ((abs(particle.at(0)) > 0 && abs(particle.at(0)) < 7) || abs(particle.at(0)) == 21)){
      outPartons.push_back(tmp4vect);
    }
    // store number of leptons
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

  //---loop on lhe events particle searching for W's---
  for(size_t iPart = 0; iPart < LHE_event->size(); iPart++){ 

    std::vector<float> particle = LHE_event->at(iPart);
    TLorentzVector tmp4vect;
    tmp4vect.SetPxPyPzE(particle.at(6), particle.at(7), particle.at(8), particle.at(9));   

    // for each particle in the LHE we can store the info in the treeWriter and the making branches in the dumper
    candidate = factory->NewCandidate();
    candidate->PID    = particle.at(0);
    candidate->Status = particle.at(1);

    // mother and  daughters are not set
    candidate->M1 = -1;
    candidate->M2 = -1;
    candidate->D1 = -1;
    candidate->D2 = -1;

    pdgParticle = pdg->GetParticle(candidate->PID);
    candidate->Charge = pdgParticle ? Int_t(pdgParticle->Charge()/3.0) : -999;
    // store mass and 4V 
    candidate->Mass = tmp4vect.M();
    candidate->Momentum.SetPxPyPzE(tmp4vect.Px(),tmp4vect.Py(),tmp4vect.Pz(),tmp4vect.E());
    candidate->Position.SetXYZT(tmp4vect.X(),tmp4vect.Y(),tmp4vect.Z(),tmp4vect.T());

    LHEparticlesArray->Add(candidate);
  }
  
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



