/*
root -l examples/Example1.C\(\"delphes_output.root\"\)

or just:

./GeneralExample delphes_output.root

N.B. you must touch GeneralExample.cpp in order to convince Make that there is anything to do to build this!
*/

//------------------------------------------------------------------------------

void GeneralExample(const char *inputFile, const char *outputFile)
{
  //  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);
  
  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  
  // Get pointers to branches used in this analysis
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  //  TClonesArray *branchPuppiJet = treeReader->UseBranch("PuppiJet");
  TClonesArray *branchRho = treeReader->UseBranch("Rho");
  TClonesArray *branchGlobalRho = treeReader->UseBranch("GlobalRho");
  TClonesArray *branchNPU = treeReader->UseBranch("NPU");
  
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");

  // Constituents will be 0 otherwise
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchEFlowTower = treeReader->UseBranch("EFlowTower");
  TClonesArray *branchEFlowMuon = treeReader->UseBranch("EFlowMuon");
  TClonesArray *branchGenParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchBeamSpotParticle = treeReader->UseBranch("BeamSpotParticle");
  //  TClonesArray *branchGenParticleWithPU = treeReader->UseBranch("ParticleWithPU");

  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
  TClonesArray *branchGenMissingET = treeReader->UseBranch("GenMissingET");
  TClonesArray *branchPileUpJetIDMissingET = treeReader->UseBranch("PileUpJetIDMissingET");
  TClonesArray *branchPuppiMissingET = treeReader->UseBranch("PuppiMissingET");


  bool verbose = true;
  bool listJetTowers = false;
  bool listMET = true;
  bool listRho = false;

  bool treeMET = true;

  TFile *f;
  TTree *t;
    //  TH1F *hdmet, hdpujidmet, hdpuppimet;
    //  TH1F *hfdmet, hfdpujidmet, hfdpuppimet;
    //  TH2F *h2met, h2pujidmet, h2puppimet;
  float genmet, met, puppimet, pujidmet;
  float genmetx, metx, puppimetx, pujidmetx;
  float genmety, mety, puppimety, pujidmety;


  if (treeMET) {
    f = new TFile(outputFile,"RECREATE");
    t = new TTree("t","t");
    t->Branch("genmet",&genmet,"genmet/F");
    t->Branch("met",&met,"met/F");
    t->Branch("pujidmet",&pujidmet,"pujidmet/F");
    t->Branch("puppimet",&puppimet,"puppimet/F");
    t->Branch("genmetx",&genmetx,"genmetx/F");
    t->Branch("metx",&metx,"metx/F");
    t->Branch("pujidmetx",&pujidmetx,"pujidmetx/F");
    t->Branch("puppimetx",&puppimetx,"puppimetx/F");
    t->Branch("genmety",&genmety,"genmety/F");
    t->Branch("mety",&mety,"mety/F");
    t->Branch("pujidmety",&pujidmety,"pujidmety/F");
    t->Branch("puppimety",&puppimety,"puppimety/F");

  }

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
  
    if (listMET||verbose||entry%5000==0) cout << "Event " << entry << " / " << numberOfEntries << endl;
    
    for (int i = 0 ;  i < branchGenMissingET->GetEntries() ; i++) {
      MissingET *m = (MissingET*) branchGenMissingET->At(i);
      if (verbose || listMET) cout << "Gen MissingET: " << m->MET << endl;
      genmet = m->MET;
      genmetx = m->MET*cos(m->Phi);
      genmety = m->MET*sin(m->Phi);
    }

    for (int i = 0 ;  i < branchMissingET->GetEntries() ; i++) {
      MissingET *m = (MissingET*) branchMissingET->At(i);
      if (verbose || listMET) cout << "MissingET: " << m->MET << endl;
      met = m->MET;
      metx = m->MET*cos(m->Phi);
      mety = m->MET*sin(m->Phi);

    }

    for (int i = 0 ;  i < branchPileUpJetIDMissingET->GetEntries() ; i++) {
      MissingET *m = (MissingET*) branchPileUpJetIDMissingET->At(i);
      if (verbose || listMET) cout << "MissingET using PileUpJetID: " << m->MET << endl;
      pujidmet = m->MET;
      pujidmetx = m->MET*cos(m->Phi);
      pujidmety = m->MET*sin(m->Phi);
    }

    for (int i = 0 ;  i < branchPuppiMissingET->GetEntries() ; i++) {
      MissingET *m = (MissingET*) branchPuppiMissingET->At(i);
      if (verbose || listMET) cout << "Puppi MissingET: " << m->MET << endl;
      puppimet = m->MET;
      puppimetx = m->MET*cos(m->Phi);
      puppimety = m->MET*sin(m->Phi);
    }

    if (treeMET) t->Fill();

    for (int i = 0 ; i < branchRho->GetEntries() ; i++) {
      Rho *rho = (Rho*) branchRho->At(i);
      if (verbose || listRho) cout << "  Rho (" << rho->Edges[0] << "-" << rho->Edges[1] << "): " << rho->Rho << endl;
    }

    for (int i = 0 ; i < branchGlobalRho->GetEntries() ; i++) {
      Rho *rho = (Rho*) branchGlobalRho->At(i);
      if (verbose || listRho) cout << "  GlobalRho (" << rho->Edges[0] << "-" << rho->Edges[1] << "): " << rho->Rho << endl;
    }


    //    cout << "before scalarHT" << endl;
    // I have cheated and recorded the true number of pileup vertices in a "ScalarHT" object!
    ScalarHT *NPU = (ScalarHT*) branchNPU->At(0);
    int nPUvertices_true = (int)NPU->HT;
    if (verbose) cout << "  Number of true pileup vertices: " << nPUvertices_true << endl;

    // One particle from primary vertex
    if (verbose && branchBeamSpotParticle) {
      cout << branchBeamSpotParticle->GetEntries() << endl;
      GenParticle *part = (GenParticle*) branchBeamSpotParticle->At(0);
      cout <<  "  True primary vertex X Y Z T: " << part->X << " " << part->Y << " " << part->Z << " " << part->T << endl;
    }

    // Status code 3 (+high pt leptons, b+t quarks, etc)  particle collection
    if (verbose && branchGenParticle) {
      for (int i = 0 ; i < branchGenParticle->GetEntries() ; i++ ) {
	GenParticle *part = (GenParticle*) branchGenParticle->At(i);
	cout << "     Status code" << part->Status << " generator particle PID Pt Eta Phi Z T (at origin) "  << part->PID << " "
	     << part->PT << " " << part->Eta << " " << part->Phi << " " << part->Z << " " << part->T << endl;
      }
    }
      
    if (verbose) {

      for (int i = 0 ; i < branchElectron->GetEntries() ; i++) {
	Electron *ele = (Electron*) branchElectron->At(i);
	cout << "    Electron " << i << ": PT Eta Phi Isolation " << ele->PT << " " << ele->Eta << " " << ele->Phi << " " << ele->IsolationVar << endl;
	GenParticle *part = (GenParticle*) ele->Particle.GetObject();
	if (part) {
	  cout << "      Electron matches to generated particle with Status code" << part->Status << " generator particle PID Pt Eta Phi Z T (at origin) "  << part->PID << " "
	       << part->PT << " " << part->Eta << " " << part->Phi << " " << part->Z << " " << part->T << endl;
	} else {
	  cout << "      Electron does not match to a generated particle" << endl;
	}
      }

      for (int i = 0 ; i < branchPhoton->GetEntries() ; i++) {
        Photon *pho = (Photon*) branchPhoton->At(i);
        cout << "    Photon " << i << ": PT Eta Phi Isolation T " << pho->PT << " " << pho->Eta << " " << pho->Phi << " " << pho->IsolationVar << " " << pho->TOuter << endl;
      }

      for (int i = 0 ; i < branchMuon->GetEntries() ; i++) {
        Muon *mu = (Muon*) branchMuon->At(i);
        cout << "    Muon " << i << ": PT Eta Phi Isolation " << mu->PT << " " << mu->Eta << " " << mu->Phi << " " << mu->IsolationVar << endl;
      }

      for (int i = 0 ; i < branchGenJet->GetEntries() ; i++) {
        Jet *jet = (Jet*) branchGenJet->At(i);
        if (jet->PT > 30.) {
          cout << "  Gen Jet " << i << endl;
          cout << "    pT: " << jet->PT << endl;
          cout << "    Eta: " << jet->Eta << endl;
          cout << "    Phi: " << jet->Phi << endl;
	  /*
          cout << "    Area: " << jet->AreaP4().Pt() << endl;
          cout << "    Constituents: " << jet->Constituents.GetEntries() << endl;
          if (branchGenParticleWithPU) { // now same colection
            for (int j = 0 ; j < jet->Constituents.GetEntries() ; j++) {
              TObject *obj = jet->Constituents[j];
              if (obj && obj->IsA() == GenParticle::Class()) {
                GenParticle *part = static_cast<GenParticle *> ( obj ) ;
                cout << "     Jet constituent Pt Eta Phi Z T (at origin) " << part->PT << " " << part->Eta << " " << part->Phi << " " << part->Z << " " << part->T << endl;
              } else {
                cout << "     Jet constituent is not a particle (?)" << endl;
              }
            }
          }
	  */
        }
      }

      // Loop over jets
      for (int i = 0 ; i < branchJet->GetEntries() ; i++) {
	Jet *jet = (Jet*) branchJet->At(i);
	if (jet->PT > 30.) {
	  cout << "  Jet " << i << endl;
	  cout << "    pT: " << jet->PT << endl;
	  cout << "    Eta: " << jet->Eta << endl;
	  cout << "    BTag: " << bool(jet->BTagAlgo&1) << " | " << bool(jet->BTagAlgo&2) << endl;
          cout << "    TauTag: " << jet->TauTag << endl;
	  //	  cout << "    Area: " << jet->AreaP4().Pt() << endl;
	  cout << "    Jet Pileup ID" << endl;
	  //	  cout << "      Beta*: " << jet->BetaStar << endl;
	  cout << "      Fractional pT in annuli (<0.1, 0.1-0.2, ..., 0.4-0.5) " << jet->FracPt[0] << " " << jet->FracPt[1] << " " << jet->FracPt[2] << " " << jet->FracPt[3] << " " << jet->FracPt[4] << endl;
	  //	  cout << "      <dR^2>: " << jet->MeanSqDeltaR << endl;
	  //          cout << "      NNeutrals: " << jet->NNeutrals << endl;
	  //	  cout << "      NCharged: " << jet->NCharged << endl;
	  cout << "    Number of constituents: " << jet->Constituents.GetEntries() << endl;
	  if (listJetTowers && branchEFlowTrack && branchEFlowTower && branchEFlowMuon) {
	    for (int j = 0 ; j < jet->Constituents.GetEntries() ; j++) {
	      TObject *obj = jet->Constituents[j];
	      if (obj && obj->IsA() == Tower::Class()) {
		Tower *tow = static_cast<Tower *> ( obj ) ;
		cout << "     Jet constituent Et Eta Phi Time (at calo) " << tow->ET << " " << tow->Eta << " " << tow->Phi << " " << tow->TOuter << endl;
	      } else {
		//		cout << "  not a tower - could check if it's a track instead (cf. Example3.C)" << endl;
	      }
	    }
	  }
	}
      }
    } // verbose 

  } // event

  if (treeMET) {
    f->cd();
    t->Write();
    f->Close();
  }

}
