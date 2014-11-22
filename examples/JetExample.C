/*
root -l examples/Example1.C\(\"delphes_output.root\"\)

or just:

./JetExample delphes_output.root

N.B. you must touch JetExample.cpp in order to convince Make that there is anything to do to build this!
*/

//------------------------------------------------------------------------------

void JetExample(const char *inputFile,const char *outputFile)
{
  //  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);
  
  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  
  // Get pointers to branches used in this analysis
  TClonesArray *branchRawJetNoPU = treeReader->UseBranch("RawJetNoPU");
  TClonesArray *branchRawJet = treeReader->UseBranch("RawJet");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchPuppiJet = treeReader->UseBranch("PuppiJet");
  TClonesArray *branchSubtractedPuppiJet = treeReader->UseBranch("SubtractedPuppiJet");

  TClonesArray *branchRho = treeReader->UseBranch("Rho");
  TClonesArray *branchPuppiRho = treeReader->UseBranch("PuppiRho");
  TClonesArray *branchNPU = treeReader->UseBranch("NPU");

  TClonesArray *branchGenJetWithPU = treeReader->UseBranch("GenJetWithPU");


  // Constituents will be 0 otherwise
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchEFlowTower = treeReader->UseBranch("EFlowTower");
  TClonesArray *branchEFlowMuon = treeReader->UseBranch("EFlowMuon");
  TClonesArray *branchGenParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchBeamSpotParticle = treeReader->UseBranch("BeamSpotParticle");
  TClonesArray *branchGenParticleWithPU = treeReader->UseBranch("ParticleWithPU");

  TClonesArray *branchPuppiWeightedParticles = treeReader->UseBranch("PuppiWeightedParticles");

  bool verbose = false;
  bool puppiConstituents = false;
  bool doResTree = true;
  bool doFakeTree = true;
  bool doPuppiResTree = true;
  bool doPuppiFakeTree = true;

  float rt_geneta, rt_genpt, rt_recpt, rt_dr;
  float ft_receta, ft_genpt, ft_recpt, ft_dr;
  float prt_geneta, prt_genpt, prt_recpt, prt_dr;
  float pft_receta, pft_genpt, pft_recpt, pft_dr;
  TFile *f;
  TTree *rt, *ft, *prt, *pft;

  //  float ft_BetaStar, ft_MeanSqDeltaR;

  if (doResTree||doFakeTree||doPuppiResTree||doPuppiFakeTree) {
    f = new TFile(outputFile,"RECREATE");
  }

  if (doResTree) {
    rt = new TTree("rt","rt");
    rt->Branch("geneta",&rt_geneta,"geneta/F");
    rt->Branch("genpt",&rt_genpt,"genpt/F");
    rt->Branch("recpt",&rt_recpt,"recpt/F");
    rt->Branch("dr",&rt_dr,"dr/F");
  }

  if (doPuppiResTree) {
    prt = new TTree("prt","prt");
    prt->Branch("geneta",&prt_geneta,"geneta/F");
    prt->Branch("genpt",&prt_genpt,"genpt/F");
    prt->Branch("recpt",&prt_recpt,"recpt/F");
    prt->Branch("dr",&prt_dr,"dr/F");
  }

  if (doFakeTree) {
    ft = new TTree("ft","ft");
    ft->Branch("receta",&ft_receta,"geneta/F");
    ft->Branch("genpt",&ft_genpt,"genpt/F");
    ft->Branch("recpt",&ft_recpt,"recpt/F");
    ft->Branch("dr",&ft_dr,"dr/F");
    //    ft->Branch("BetaStar",&ft_BetaStar,"BetaStar/F");
    //    ft->Branch("MeanSqDeltaR",&ft_MeanSqDeltaR,"MeanSqDeltaR/F");
  }

  if (doPuppiFakeTree) {
    pft = new TTree("pft","pft");
    pft->Branch("receta",&pft_receta,"geneta/F");
    pft->Branch("genpt",&pft_genpt,"genpt/F");
    pft->Branch("recpt",&pft_recpt,"recpt/F");
    pft->Branch("dr",&pft_dr,"dr/F");
  }




  // Book histograms

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
  
    if (puppiConstituents||verbose||entry%5000==0) cout << "Event " << entry << " / " << numberOfEntries << endl;
    
    if (puppiConstituents) {

      if (branchEFlowTrack && branchEFlowTower) {
	for (int i = 0 ; i < branchEFlowTrack->GetEntries() ; i++) {
	  Track *track = (Track*) branchEFlowTrack->At(i);
	  if (track->IsRecoPU == 0) {
	    cout << " Input Primary EFlow Track PT Eta Phi " << track->PT << " " << track->Eta << " " << track->Phi << endl;
	  } else {
            cout << " Input PileUp EFlow Track PT Eta Phi " << track->PT << " " << track->Eta << " "<< track->Phi << endl;
	  }
	}
        for (int i = 0 ; i < branchEFlowTower->GetEntries() ; i++) {
          Tower *tower = (Tower*) branchEFlowTower->At(i);
	  cout << " Input EFlow Tower ET Eta Phi " << tower->ET << " " << tower->Eta << " " << tower->Phi << endl;
	}
      }

      if (branchPuppiWeightedParticles) {
        for (int i = 0 ; i < branchPuppiWeightedParticles->GetEntries() ; i++) {
          GenParticle *part = (GenParticle*) branchPuppiWeightedParticles->At(i);
	  cout << " PUPPI weighted particle PT Eta Phi M " << part->PT << " " << part->Eta << " " << part->Phi << " " << part->Mass << endl;
        }
      }

      if (branchGenParticle) {
        for (int i = 0 ; i < branchGenParticle->GetEntries() ; i++) {
          GenParticle *part = (GenParticle*) branchGenParticle->At(i);
	  if (part->IsPU == 0) {
	    cout << " PV GenParticle PT Eta Phi M " << part->PT << " " << part->Eta << " " << part->Phi << " " << part->Mass << endl;
	  } else {
            cout << " Pileup (vtx " << part->IsPU << ") GenParticle PT Eta Phi M " << part->PT << " " << part->Eta << " " << part->Phi << " " << part->Mass << endl;
	  }
        }
      }

    }



    for (int i = 0 ; i < branchRho->GetEntries() ; i++) {
      Rho *rho = (Rho*) branchRho->At(i);
      if (verbose) cout << "  Rho (" << rho->Edges[0] << "-" << rho->Edges[1] << "): " << rho->Rho << endl;
    }

    for (int i = 0 ; i < branchPuppiRho->GetEntries() ; i++) {
      Rho *rho = (Rho*) branchPuppiRho->At(i);
      if (verbose) cout << "  PUPPI Rho (" << rho->Edges[0] << "-" << rho->Edges[1] << "): " << rho->Rho << endl;
    }


    //    cout << "before scalarHT" << endl;
    // I have cheated and recorded the true number of pileup vertices in a "ScalarHT" object!
    ScalarHT *NPU = (ScalarHT*) branchNPU->At(0);
    int nPUvertices_true = (int)NPU->HT;
    if (verbose) cout << "  Number of true pileup vertices: " << nPUvertices_true << endl;

    if (doResTree) {
      // some very simple variables in a very simple tree
      for (int i = 0 ; i < branchGenJet->GetEntries() ; i++) {
	Jet *genjet = (Jet*) branchGenJet->At(i);
	TLorentzVector genp4 = genjet->P4();
	rt_genpt = genjet->PT;
	rt_geneta = genjet->Eta;
	rt_dr = 9999.;
	rt_recpt = -1.;
	for (int j = 0 ; j < branchJet->GetEntries() ; j++ ) {
	  Jet *jet = (Jet*) branchJet->At(j);
	  TLorentzVector p4 = jet->P4();
	  float tempdr = p4.DeltaR(genp4);
	  if (tempdr < rt_dr) {
	    rt_dr = tempdr;
	    rt_recpt = jet->PT;
	  }
	}
	rt->Fill(); // per genjet
      }
    } // doResTree

    if (doFakeTree) {
      // some very simple variables in a very simple tree 
      for (int i = 0 ; i < branchJet->GetEntries() ; i++) {
        Jet *jet = (Jet*) branchJet->At(i);
        TLorentzVector p4 = jet->P4();
        ft_recpt = jet->PT;
        ft_receta = jet->Eta;
        ft_dr = 9999.;
        ft_genpt = -1.;
	//	ft_BetaStar = jet->BetaStar;
	//	ft_MeanSqDeltaR = jet->MeanSqDeltaR;
        for (int j = 0 ; j < branchGenJet->GetEntries() ; j++ ) {
          Jet *genjet = (Jet*) branchGenJet->At(j);
          TLorentzVector genp4 = genjet->P4();
          float tempdr = p4.DeltaR(genp4);
          if (tempdr < ft_dr) {
            ft_dr = tempdr;
            ft_genpt = jet->PT;
          }
        }
      ft->Fill(); // per recjet
      }
    } // doFakeTree
                                                                                                                                                                 
    if (doPuppiResTree) {
      for (int i = 0 ; i < branchGenJet->GetEntries() ; i++) {
        Jet *genjet = (Jet*) branchGenJet->At(i);
        TLorentzVector genp4 = genjet->P4();
        prt_genpt = genjet->PT;
        prt_geneta = genjet->Eta;
        prt_dr = 9999.;
        prt_recpt = -1.;
        for (int j = 0 ; j < branchPuppiJet->GetEntries() ; j++ ) {
          Jet *jet = (Jet*) branchPuppiJet->At(j);
          TLorentzVector p4 = jet->P4();
          float tempdr = p4.DeltaR(genp4);
          if (tempdr < prt_dr) {
            prt_dr = tempdr;
            prt_recpt = jet->PT;
          }
        }
        prt->Fill(); // per genjet
      }
    } // doPuppiResTree                                                                                                                                                                   
    if (doPuppiFakeTree) {
      for (int i = 0 ; i < branchPuppiJet->GetEntries() ; i++) {
	Jet *jet = (Jet*) branchPuppiJet->At(i);
        TLorentzVector p4 = jet->P4();
        pft_recpt = jet->PT;
        pft_receta = jet->Eta;
        pft_dr = 9999.;
        pft_genpt = -1.;
        for (int j = 0 ; j < branchGenJet->GetEntries() ; j++ ) {
          Jet *genjet = (Jet*) branchGenJet->At(j);
          TLorentzVector genp4 = genjet->P4();
          float tempdr = p4.DeltaR(genp4);
          if (tempdr < pft_dr) {
            pft_dr = tempdr;
            pft_genpt = jet->PT;
          }
        }
        pft->Fill(); // per puppijet
      }
    }  

    //    cout << "before branchBeamSpotParticle" << endl;
    // One particle from primary vertex
    if (verbose && branchBeamSpotParticle) {
      cout << branchBeamSpotParticle->GetEntries() << endl;
      GenParticle *part = (GenParticle*) branchBeamSpotParticle->At(0);
      cout <<  "  True primary vertex X Y Z T: " << part->X << " " << part->Y << " " << part->Z << " " << part->T << endl;
      //      cout << "     Generator particle PID Pt Eta Phi X Y Z T (at origin) "  << part->PID << " "
      //	   << part->PT << " " << part->Eta << " " << part->Phi << " " << part->X << " " << part->Y << " " << part->Z << " " << part->T << endl;
    }
    //    cout << "after branchBeamSpotParticle" << endl;


    // Status code 3 particle collection
    if (verbose && branchGenParticle) {
      for (int i = 0 ; i < branchGenParticle->GetEntries() ; i++ ) {
	GenParticle *part = (GenParticle*) branchGenParticle->At(i);
	cout << "     Status code 3 generator particle PID Pt Eta Phi Z T (at origin) "  << part->PID << " "
	     << part->PT << " " << part->Eta << " " << part->Phi << " " << part->Z << " " << part->T << endl;
      }
    }
      
    if (verbose) {
      // Loop over jets
      for (int i = 0 ; i < branchJet->GetEntries() ; i++) {
	//	cout << "  Jet " << i << "/" << branchJet->GetEntries() << endl;
	Jet *jet = (Jet*) branchJet->At(i);
	//	cout << jet->PT << endl;
	//	continue;
	if (jet->PT > 100.) {
	  cout << "  Jet " << i << endl;
	  cout << "    pT: " << jet->PT << endl;
	  cout << "    Eta: " << jet->Eta << endl;
	  cout << "    Area: " << jet->AreaP4().Pt() << endl;
	  //	  cout << "    Time at calo: " << jet->t0 << endl;
	  //	  cout << "    Time at calo smeared by 1ps 10ps 20ps 30ps 40ps: " << jet->t1 << " " << jet->t10 << " " << jet->t20 << " " << jet->t30 << " " << jet->t40 << endl;
	  //	  cout << "    Beta*:" << jet->BetaStar << endl;
	  //	  cout << "    NCharged: " << jet->NCharged << endl;
	  cout << "    Number of constituents: " << jet->Constituents.GetEntries() << endl;
	  if (branchEFlowTrack && branchEFlowTower && branchEFlowMuon) {
	    for (int j = 0 ; j < jet->Constituents.GetEntries() ; j++) {
	      TObject *obj = jet->Constituents[j];
	      if (obj->IsA() == Tower::Class()) {
		Tower *tow = static_cast<Tower *> ( obj ) ;
		cout << "     Jet constituent Et Eta Phi Time (at calo) " << tow->ET << " " << tow->Eta << " " << tow->Phi << " " << tow->TOuter << endl;
		if (jet->P4().DeltaR(tow->P4()) > 0.5) cout << "        !!!! dr=" << jet->P4().DeltaR(tow->P4()) << endl;
	      } else {
		cout << "  not a tower - could check if it's a track instead (cf. Example3.C)" << endl;
	      }
	    }
	  }
	}
      }
      
      // Loop over jets
      if (branchRawJet) {
	for (int i = 0 ; i < branchRawJet->GetEntries() ; i++) {
	  Jet *jet = (Jet*) branchRawJet->At(i);
	  if (jet->PT > 100.) {
	    cout << "  Raw Jet " << i << endl;
	    cout << "    pT: " << jet->PT << endl;
	    cout << "    Eta: " << jet->Eta << endl;
	    cout << "    Area: " << jet->AreaP4().Pt() << endl;
	  }
	}
      }
      if (branchRawJetNoPU) {
	for (int i = 0 ; i < branchRawJetNoPU->GetEntries() ; i++) {
	  Jet *jet = (Jet*) branchRawJetNoPU->At(i);
	  if (jet->PT > 30.) {
	    cout << "  Raw Jet (no PU mixing) " << i << endl;
	    cout << "    pT: " << jet->PT << endl;
	    cout << "    Eta: " << jet->Eta << endl;
	    cout << "    Area: " << jet->AreaP4().Pt() << endl;
	  }
	}
      }
	
      for (int i = 0 ; i < branchGenJet->GetEntries() ; i++) {
	Jet *jet = (Jet*) branchGenJet->At(i);
	if (jet->PT > 30.) {
	  cout << "  Gen Jet " << i << endl;
	  cout << "    pT: " << jet->PT << endl;
	  cout << "    Eta: " << jet->Eta << endl;
	  cout << "    Area: " << jet->AreaP4().Pt() << endl;
	  cout << "    Constituents: " << jet->Constituents.GetEntries() << endl;
	  if (branchGenParticleWithPU) { // now same colection
	    for (int j = 0 ; j < jet->Constituents.GetEntries() ; j++) {
	      TObject *obj = jet->Constituents[j];
	      if (obj->IsA() == GenParticle::Class()) {
		GenParticle *part = static_cast<GenParticle *> ( obj ) ;
		cout << "     Jet constituent Pt Eta Phi Z T (at origin) " << part->PT << " " << part->Eta << " " << part->Phi << " " << part->Z << " " << part->T << endl;
	      } else {
		cout << "     Jet constituent is not a particle (?)" << endl;
	      }
	    }
	  }
	}
      }

      for (int i = 0 ; i < branchPuppiJet->GetEntries() ; i++) {
        Jet *jet = (Jet*) branchPuppiJet->At(i);
        if (jet->PT > 30.) {
          cout << "  PUPPI Jet " << i << endl;
          cout << "    pT: " << jet->PT << endl;
          cout << "    Eta: " << jet->Eta << endl;
	}
      }
      
      for (int i = 0 ; i < branchSubtractedPuppiJet->GetEntries(); i++) {
	Jet *jet = (Jet*) branchSubtractedPuppiJet->At(i);
	  if (jet->PT > 30.) {
	    cout << "  Subtracted PUPPI Jet " << i << endl;
	    cout << "    pT: " << jet->PT << endl;
	    cout << "    Eta: " << jet->Eta << endl;
	  }
      }
      
      if (branchGenJetWithPU) {
	for (int i = 0 ; i < branchGenJetWithPU->GetEntries() ; i++) {
	  Jet *jet = (Jet*) branchGenJetWithPU->At(i);
	  if (jet->PT > 100.) {
	    cout << "  Gen Jet (with PU included) " << i << endl;
	    cout << "    pT: " << jet->PT << endl;
	    cout << "    Eta: " << jet->Eta << endl;
	    cout << "    Area: " << jet->AreaP4().Pt() << endl;
	    cout << "    Constituents: " << jet->Constituents.GetEntries() << endl;
	    if (branchGenParticleWithPU) {
	      for (int j = 0 ; j < jet->Constituents.GetEntries() ; j++) {
		TObject *obj = jet->Constituents[j];
		if (obj->IsA() == GenParticle::Class()) {
		  GenParticle *part = static_cast<GenParticle *> ( obj ) ;
		  cout << "     Jet constituent Pt Eta Phi Z T (at origin) " << part->PT << " " << part->Eta << " " << part->Phi << " " << part->Z << " " << part->T << endl;
		} else {
		  cout << "     Jet constituent is not a particle (?)" << endl;
		}
	      }
	    }
	  }
	}
      }


    } // verbose 
  } // event

  if (doResTree||doFakeTree||doPuppiResTree||doPuppiFakeTree) {
    f->cd();
    if (doResTree) rt->Write();
    if (doFakeTree) ft->Write();
    if (doPuppiResTree) prt->Write();
    if (doPuppiFakeTree) pft->Write();
    f->Close();
  }

}
