// -*- C++ -*-
#include "Rivet/Analysis.hh"

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh" //this gives us the tracks for the particles in the event??
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"

#include "Rivet/Projections/SmearedMET.hh"// Wrapper projection for smearing missing (transverse) energy/momentum with detector resolutions.
#include "Rivet/Projections/MissingMomentum.hh"
//#include "Rivet/Projections/Jet.hh"


//David:: will this allow me to use containBottom??? Eleni:: nice question !!!
#include "Rivet/Jet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Analyses/MC_JetAnalysis.hh"

#include "Rivet/Particle.fhh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/ParticleBase.hh"


namespace Rivet {
  
  // @brief Add a short analysis description here
  class MonoH_Truth : public Analysis {
  
  public:
    
    //declare the histograms here first
    Histo1DPtr ptlj0;
    Histo1DPtr ptlj1;
    Histo1DPtr etalj0; 
    Histo1DPtr etalj1;
    Histo1DPtr numev_evw;
    Histo1DPtr bjet1_pt;
    Histo1DPtr bjet2_sublpt;
    Histo1DPtr hist_pt_trk;
    Histo1DPtr vEt_miss;
    Histo1DPtr vEt_miss_phi;
    Histo1DPtr vEt_miss_eta;
    Histo1DPtr missing_Et_phi;
    Histo1DPtr missing_Et_eta;
    Histo1DPtr missing_Et;
    Histo1DPtr track_pt;
    Histo1DPtr bjet2_lpt;
    Histo1DPtr dh_m        ;
    Histo1DPtr dh_eta      ;        
    Histo1DPtr dh_phi      ;
    Histo1DPtr dh_m_bb     ;
    Histo1DPtr dh_m_bb_jets;
    Histo1DPtr dh_m_fat_jet;
    Histo1DPtr ptlj0_all_cuts  ; 
    Histo1DPtr ptlj1_all_cuts  ;
    Histo1DPtr etalj0_all_cuts ;
    Histo1DPtr etalj1_all_cuts ;

    Histo1DPtr bjet1_pt_all_cuts    ;  
    Histo1DPtr bjet2_sublpt_all_cuts;  
    Histo1DPtr bjet2_mass_all_cuts   ;  
    Histo1DPtr dh_m_MET500        ;
    Histo1DPtr dh_eta_MET500      ;        
    Histo1DPtr dh_phi_MET500      ;
    Histo1DPtr dh_m_bb_MET500     ;
    Histo1DPtr dh_m_fat_jet_MET500;

    Histo1DPtr FatJetsSize       ;
    Histo1DPtr bjetsSize         ;

    Histo1DPtr dh_m_fat_jet_all_cuts;
    
    // Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MonoH_Truth);
    
    
    // ---------------------------------------------------------
    // init method:: here you book histograms and initialise projections 
    // before the run
    // ---------------------------------------------------------
    void init() {
      
      auto vfs = VisibleFinalState(Cuts::abseta < 4.9);
      declare(vfs, "VFS");

      // Initialise and register projections
      
      auto fs = FinalState(Cuts::abseta < 5 && Cuts::pT > 100*MeV);
      declare(fs, "FS");
      
      //small R-jets
      declare(FastJets(fs, FastJets::ANTIKT, 0.4), "Jets");
      
      //large R-jets
      declare(FastJets(fs, FastJets::ANTIKT, 1.0), "Fat_Jets");      
      
      //Track Jets
      declare(FastJets(fs, FastJets::TRACKJET, 0.2), "Track_Jets");


      declare(MissingMomentum(fs), "MissingPt");
      declare(MissingMomentum(fs), "MissingET");
      declare(VisibleFinalState(Cuts::abseta <4.5 && Cuts::pT > 500*MeV), "Particle_Tracks");
      //ChargedFinalState

      // Book histograms
      ptlj0                    = bookHisto1D("ptlj0", 100, -0.5, 1500,"", "$p_T^{leading~jet}[GeV]$", "Sum of event weights");
      ptlj1                    = bookHisto1D("ptlj1", 100, -0.5, 1500,"", "$p_T^{subleading~jet}[GeV]$", "Sum of event weights");
      etalj0                   = bookHisto1D("etalj0", 40, -4,4,"","$eta^{leading~jet}$", "Sum of event weights" );
      etalj1                   = bookHisto1D("etalj1", 40, -4, 4,"", "$eta^{subleading~jet}$", "Sum of event weights" );
      ptlj0_all_cuts           = bookHisto1D( "ptlj0_all_cuts", 100, -0.5, 1500,"", "$p_T^{leading~jet}[GeV]$ after cuts", "Sum of event weights");
      ptlj1_all_cuts           = bookHisto1D( "ptlj1_all_cuts", 100, -0.5, 1500,"", "$p_T^{subleading~jet}[GeV]$ after cuts", "Sum of event weights");
      etalj0_all_cuts          = bookHisto1D("etalj0_all_cuts", 40, -4,4,"","$eta^{leading~jet}$ after cuts", "Sum of event weights" );
      etalj1_all_cuts          = bookHisto1D("etalj1_all_cuts", 40, -4, 4,"", "$eta^{subleading~jet}$ after cuts", "Sum of event weights" );

      numev_evw                = bookHisto1D("numev_evw", 50, 1e-9, 1.3e-9 ,"", "event number ", "Event weight");
             
      bjet1_pt                 = bookHisto1D("bjet1_pt", 70,0,7, "","Delta R(Lead Fat, Lead b) (1 bj/event)", "Sum of event weight");
      bjet2_sublpt             = bookHisto1D("bjet2_sublpt", 70, -0.5, 7,"", "Delta R(Lead Fat, Sublead b)(2 bj/event)", "Sum of event weight");
      bjet2_lpt                = bookHisto1D("bjet2_lpt", 70,-0.5, 7,"", "Delta R (Lead Fat, Lead b)(2 bj/event)", "Sum of event weight");

      bjet1_pt_all_cuts        = bookHisto1D("bjet1_pt_all_cuts", 70,0,7, "","Delta R(Lead Fat, Lead b) (1 bj/event) (after cuts) ", "Sum of event weight");
      bjet2_sublpt_all_cuts    = bookHisto1D("bjet2_sublpt_all_cuts", 70, -0.5, 7,"", "Delta R(Lead Fat, Lead b)(2 bj/event)  (after cuts)", "Sum of event weight");
      bjet2_mass_all_cuts       = bookHisto1D("bjet2_mass_all_cuts", 100,0, 1000,"", "$M_{bb}$  (after cuts) [GeV]", "Event number");

      missing_Et               = bookHisto1D("missing_Et", 50,-0.5, 1500,"", "$P^{miss}_T$[GeV] (calc.manually)", "Sum of event weight");
      missing_Et_phi           = bookHisto1D("missing_Et_phi", 50,0, 7,"", "phi of missing 4momentum (calc.manually)", "Sum of event weight");
      missing_Et_eta           = bookHisto1D("missing_Et_eta", 50,-5, 5,"", "eta of missing 4momentum (calc.manually)", "Sum of event weight");
      vEt_miss                 = bookHisto1D("vEt_miss", 50,-0.5, 1500,"", "vec(P$^{miss}_T)[GeV] $ (using projection)", "Sum of event weight"); 
      vEt_miss_phi             = bookHisto1D("vEt_miss_phi", 50,0,7,"",   "phi of missing momentum (using projection)", "Sum of event weight");  
      vEt_miss_eta             = bookHisto1D("vEt_miss_eta", 50,-5, 5,"", "eta of missing momentum (using projection)", "Sum of event weight");  

      track_pt                 = bookHisto1D("track_pt", 50, -0.5, 1500,"", "P^{trk}_T[GeV]$", "Sum of event weights");
      
      //histograms for the reconstruction of the dark Higgs mass 
      dh_m                     = bookHisto1D("dh_m", 100,0, 1000,"", "$M_{dh,reco} [GeV]$", "Events number");				 
      dh_m_bb                  = bookHisto1D("dh_m_bb", 100,0 , 1000, "", "$M_{bb}(track jets) [GeV]$", "Events number");	    		 
      dh_m_bb_jets             = bookHisto1D("dh_m_bb_jets", 100,0 , 1000, "", "$M_{bb} (jets) [GeV]$", "Events number");	    		 
      dh_m_fat_jet             = bookHisto1D("dh_m_fat_jet", 100, 0, 1000, "", "$M_{fat~jet} [GeV]$", "Events number");			 
      dh_eta                   = bookHisto1D("dh_eta", 50,-5, 5,"", "$eta_{dh,reco} [GeV]$", "Events number");				 
      dh_phi                   = bookHisto1D("dh_phi", 50,0,7,"", "$phi_{dh,recp} [GeV]$", "Events number");					 
			                                                                                                                      
      dh_m_MET500              = bookHisto1D("dh_inv_m_MET500", 100,0, 1000,"" ,"$M_{dh,reco}(met>500)[GeV]$ ", "Events number");		 
      dh_eta_MET500            = bookHisto1D("dh_eta_MET500", 50,-5, 5,"", "$eta_{dh,reco} [GeV]$", "Events number");			 
      dh_phi_MET500            = bookHisto1D("dh_phi_MET500", 50,0,7,"", "$phi_{dh,recp} [GeV]$", "Events number");				 
			                                                                                                                      
      dh_m_bb_MET500           = bookHisto1D("dh_inv_m_bb_MET500", 100,0 , 1000,"", "$M_{bb}(met>500) [GeV]$", "Events number");	    	 
      dh_m_fat_jet_MET500      = bookHisto1D("dh_inv_m_fat_jet_MET500", 100, 0, 1000,"", "$M_{fat~jet}(met>500) [GeV]$", "Events number");	 
    			                                                                                                                      
      FatJetsSize              = bookHisto1D("FatJetsSize", 10,0, 10,"" ,"Fat jets Number", "Events number");				 
      bjetsSize                = bookHisto1D("bjetsSize", 10,0, 10,"" ,"b jets nunmber", "Events number");                                    
      
      dh_m_fat_jet_all_cuts    = bookHisto1D("dh_m_fat_jet_all_cuts", 100,0, 1000,"", "$M_{dh,leading fat jet} (all cuts) [GeV]$", "Events number");
      
    }// end of init function 
    
    

    
    // ---------------------------------------------------------
    // This method loops in the events so we can perform the per-event analysis
    // ---------------------------------------------------------
 
    //defining the counters that are used to calculate the efficiency.
    double numEvent_met500 = 0;
    double numEvent_PT_miss_trk30 = 0;
    double numEvent_min_Dphi = 0;
    double numEvent_max_Dphi = 0;
    double numEvent_fat_jets = 0;
    double numEvent_bjet_veto = 0;
    double numEvent_nonbjet_veto = 0;
    double numEvent_HT = 0;
    double numEvent_mass_reco = 0;
    double numEvent_1bjet = 0;
    double numEvent_2bjet = 0;
    double _2bjetcounter = 0;

    void analyze(const Event& event) {
      
      //Apply Projections for small R jets, large R jets and track jets 
      const Jets jets = applyProjection<FastJets>(event, "Jets").jetsByPt((Cuts::pT > 20*GeV && Cuts::abseta < 2.5) || 
									  ( Cuts::abseta > 2.5 && Cuts::abseta < 4.5 && Cuts::pT > 30*GeV ));
      const Jets fat_jets = applyProjection<FastJets>(event, "Fat_Jets").jetsByPt(Cuts::pT > 200*GeV && Cuts::abseta < 2.0);
      const Jets track_jets = applyProjection<FastJets>(event, "Track_Jets").jetsByPt(Cuts::pT > 10*GeV && Cuts::abseta < 2.5);

      //const Vector3& vectorPt = applyProjection<MissingMomentum>(event, "MissingPt").vectorPt();
      //const Vector3& vector_Et_miss = -vectorPt;
      
      const MissingMomentum& met2 = applyProjection<MissingMomentum>(event, "MissingET"); 
      FourMomentum met1 = met2.missingMomentum();
 
      numev_evw->fill(event.weight());// fill this plot with event weight of all events for debugging purposes 
      
      //   fill histos for dark higgs mass 
      FatJetsSize -> fill(fat_jets.size());
      
      if(fat_jets.size()>=1) {
	dh_m->fill(fat_jets[0].mass());
	dh_eta->fill(fat_jets[0].eta());
	dh_phi->fill(fat_jets[0].phi());
	dh_m_fat_jet -> fill(fat_jets[0].mass());
      }
      Jets hbjets;
      foreach(const Jet & myjets, track_jets){ 
	if(myjets.containsBottom()) hbjets.push_back(myjets);
      }
      bjetsSize -> fill(hbjets.size() );
      
      Jets notrJets;
      foreach(const Jet & myjets, jets){
        if(myjets.containsBottom()) notrJets.push_back(myjets);
      }
      if(notrJets.size()>=2) dh_m_bb_jets->fill((notrJets[0].mom()+notrJets[1].mom()).mass());	    


      if(hbjets.size()>=2){
	dh_m->fill(add(hbjets[0].mom(),hbjets[1].mom()).mass());
	dh_eta->fill(add(hbjets[0].mom(),hbjets[1].mom()).eta());
	dh_phi->fill(add(hbjets[0].mom(),hbjets[1].mom()).phi());
	dh_m_bb->fill(add(hbjets[0].mom(),hbjets[1].mom()).mass());	    
      }        

      vEt_miss->fill(met2.vectorEt().mod(), event.weight());
      vEt_miss_eta->fill(met1.eta(), event.weight());
      vEt_miss_phi->fill(met1.phi(), event.weight());

      if (jets.size()>0){
	ptlj0->fill(jets[0].pT(), event.weight());
	etalj0->fill(jets[0].eta(), event.weight());
	if(jets.size()>1){
	ptlj1->fill(jets[1].pT(), event.weight());
	etalj1->fill(jets[1].eta(), event.weight());
	}
      }
      
      if(hbjets.size()>0 && fat_jets.size()>0) 
      	bjet1_pt->fill(deltaR(fat_jets[0],hbjets[0]), event.weight());  
      if(hbjets.size()>1 && fat_jets.size()>0) {
	bjet2_sublpt->fill(deltaR(fat_jets[0],hbjets[1]), event.weight());
	bjet2_lpt->fill(deltaR(fat_jets[0],hbjets[0]), event.weight());
      }

      // pTmiss      
      Particles vfs_particles = apply<VisibleFinalState>(event,"VFS").particles();
      FourMomentum pTmiss;      
      foreach ( const Particle & p, vfs_particles ) {        
	pTmiss -= p.momentum(); 
      }     
      double eTmiss = pTmiss.pT();      
      double eTmissPhi = pTmiss.phi();
      
      missing_Et->fill(eTmiss, event.weight());
      missing_Et_phi->fill(eTmissPhi, event.weight());
      missing_Et_eta->fill(pTmiss.eta(), event.weight());
      
      
      //=====================================================================================
      //            C U T S
      
      if(eTmiss<500)vetoEvent;  
      numEvent_met500 +=1;
      
      
      //   fill histos for dark higgs mass 
      
      if(fat_jets.size()>=1) {
	dh_m_MET500->fill(fat_jets[0].mass());
	dh_m_fat_jet_MET500->fill(fat_jets[0].mass());
	dh_eta_MET500->fill(fat_jets[0].eta());
	dh_phi_MET500->fill(fat_jets[0].phi());
      }
      Jets hbjets_met;
      foreach(const Jet & jets, track_jets){ 
	if(jets.containsBottom()) hbjets_met.push_back(jets);
      }
      if(hbjets_met.size()>=2){
	dh_m_MET500->fill((hbjets_met[0].mom()+hbjets_met[1].mom()).mass());
	dh_m_bb_MET500->fill((hbjets_met[0].mom()+hbjets_met[1].mom()).mass());
	dh_eta_MET500->fill((hbjets_met[0].mom()+hbjets_met[1].mom()).eta());
	dh_phi_MET500->fill((hbjets_met[0].mom()+hbjets_met[1].mom()).phi());
      }     
      

      
      // applying projection charged finalstate and particles in order to calculate missing track pT
      const VisibleFinalState& Particle_Tracks = apply<VisibleFinalState>(event, "Particle_Tracks");
      const Particles& event_parts = Particle_Tracks.particles();
      Vector3 pt_miss_trk;
      
      //summing up all the track pT
      foreach(const Particle&  p, event_parts){
	pt_miss_trk += p.momentum().pTvec();  	
      }

      track_pt->fill(pt_miss_trk.mod());

     
      
      double min_Dphi= 999.00;
      int maxindex = 3;
      if(jets.size()<3){
	maxindex = jets.size();}
      for (int jetindex = 0; jetindex< maxindex; jetindex++){
	if (deltaPhi( met1.phi(), jets.at(jetindex).momentum().pTvec().phi()) < min_Dphi){
	  min_Dphi= jets.at(jetindex).momentum().pTvec().phi();
	  
	} 
      }
      
      if (deltaPhi( met1.phi(), min_Dphi) < pi/9) vetoEvent;
     
	numEvent_min_Dphi += 1;
      
      //applies the Delta phi >pi/2 cut
      //if (deltaPhi( met1.phi(), pt_miss_trk.phi()) > PI/2) vetoEvent;//check the interval of the deltaphi function
      //numEvent_max_Dphi +=1;
      
      
      //require at least one large R jet. plot number of fat jets i have
      if (fat_jets.size()<1) vetoEvent;
      numEvent_fat_jets +=1;
      
      
      //applying the b-jet veto, requiring all non-associated jets to not contains any bottom
      foreach(const Jet & j, jets){
	if (j.abseta()<2.5)
	  if (deltaR(fat_jets[0].momentum(),j.momentum()) >= 1) 
	    if(j.containsBottom())
	      vetoEvent;
      }
      
      numEvent_bjet_veto += 1;
     
      
      //Ht ratio(0.57) cut , Ht is the pt of the large R leading jet + the sum of the non associated jet pt
      //Summing the non associated jet pt to get HT for  nonassociated jet 
      double HT_non_associated_jet = 0;
      foreach(const Jet & j, jets){
	if (deltaR(fat_jets[0],j) >= 1){
	  HT_non_associated_jet += j.pT();
	}
      }
      
      //defining the new variables that will be used to the HT-rato cut.
      double HT;
      double HT_ratio;
      HT= HT_non_associated_jet + fat_jets[0].pt();
      HT_ratio = HT_non_associated_jet / HT;
      
      if(HT_ratio >0.57) vetoEvent;
      numEvent_HT += 1;
      
      
      // takes all the particles and sort them by Pt, then acording to the code liked above
      // if we take all_particles.size() we will get the number of tracks
      //const Particles& all_particles = Particle_Tracks.particlesByPt();
      
      if (fat_jets[0].mass()<50*GeV && fat_jets[0].mass()>270*GeV) vetoEvent;
      numEvent_mass_reco +=1;
      
      //Finding b-jets. do this 
      Jets bjets, no_bjets;
      foreach(const Jet & jets, track_jets) {
	if(jets.containsBottom()){
	  bjets.push_back(jets);}

      } 
      //*if(bjets.size()>=2){
      //_2bjetcounter +=1;
	// }
      
      // plotting the mass of the leading fat jet
      dh_m_fat_jet_all_cuts->fill(fat_jets[0].mass());
      
      //const size_t njets = jets.size();
      //const size_t nbjets = bjets.size();     
      
      //here we check if we have jets (and if we have 2 or more) in the event and we take the leading (and subleading jet). 
      if (jets.size()>0){
	ptlj0_all_cuts->fill(jets[0].pT(), event.weight());
	etalj0_all_cuts->fill(jets[0].eta(), event.weight());
	if(jets.size()>1){
	ptlj1_all_cuts->fill(jets[1].pT(), event.weight());
	etalj1_all_cuts->fill(jets[1].eta(), event.weight());
	}
      }
      
      
      // making cut for 
      switch(bjets.size()){ 
      case (1) : 
	bjet1_pt_all_cuts->fill(deltaR(fat_jets[0],bjets[0]), event.weight());  
        numEvent_1bjet +=1;
        break;
      case (2) : 
	bjet2_sublpt_all_cuts->fill(deltaR(fat_jets[0],bjets[1]), event.weight());
	bjet2_mass_all_cuts->fill(add(bjets[0].mom(),bjets[1].mom()).mass());
	numEvent_2bjet +=1;
	break;
      };
      
      //Filling the vector et missing thing that acctually is the -vector pt...
      
    }//end of alanyze  
    
    
    
    // ---------------------------------------------------------
    // In this method we normalise histograms etc., after the run 
    // ---------------------------------------------------------
    
    void finalize() {
      
      //printing the efficiency of the different cuts
      cout<< "Table of the efficiency of the different cuts: \n";
      cout<< "the efficiency for Et_miss >500 GeV cut: "<<numEvent_met500/numEvents()<< "\n";
      cout<< "the efficiency for Missing Pt >30 GeV cut: "<<numEvent_PT_miss_trk30/numEvents()<< "\n";
      cout<< "the efficiency for Delta phi > pi/9 cut: "<<numEvent_min_Dphi/numEvents()<< "\n";
      cout<< "the efficiency for Delta phi < pi/2 cut: "<<numEvent_max_Dphi/numEvents()<< "\n";
      cout<< "the efficiency for nr fat jets >= 1 cut: "<<numEvent_fat_jets/numEvents()<< "\n";
      cout<< "the efficiency for the b-jet veto cut : "<< numEvent_bjet_veto/numEvents()<< "\n";
      cout<< "the efficiency for the HT-ratio cut: " <<numEvent_HT/numEvents()<< "\n";
      cout<< "the efficiency for the reconstructed mass cut: " <<numEvent_mass_reco/numEvents()<< "\n";
      cout<< "the efficiency for the 1 bjet cut: " <<numEvent_1bjet/numEvents()<< "\n";
      cout<< "the efficiency for the 2 bjet cut: " <<numEvent_2bjet/numEvents()<< "\n";


      // printing a "table" of how many events there are after every cut is made
      cout<< "Table of # Events left:\n";
      cout<< "No cuts: "<< numEvents() << "\n";
      cout<< "Missing Et >500 GeV cut: " << numEvent_met500 << "\n";
      cout<< "Missing Pt >30 GeV cut: "<< numEvent_PT_miss_trk30<< "\n";
      cout<< "Delta phi > pi/9 cut: "<< numEvent_min_Dphi << "\n";
      cout<< "Delta phi < pi/2 cut: "<< numEvent_max_Dphi << "\n";
      cout<< "nr fat jets > 0 cut: " <<numEvent_fat_jets << "\n";
      cout<< "b-jet veto cut: " << numEvent_bjet_veto << "\n";
      cout<< "HT -ratio cut: "<<numEvent_HT << "\n";
      cout<< "reconstructed mass cut: "<< numEvent_mass_reco << "\n";
      cout<< "1 bjet cut: "<< numEvent_1bjet << "\n";
      cout<< "2 bjet cut: "<< numEvent_2bjet << "\n";
      cout<< "number of 2 b-jets" << _2bjetcounter <<endl;

      normalize(ptlj0);                   
      normalize(ptlj1);                   
      normalize(etalj0);                  
      normalize(etalj1);                  
      normalize(ptlj0_all_cuts);          
      normalize(ptlj1_all_cuts);          
      normalize(etalj0_all_cuts);         
      normalize(etalj1_all_cuts);         
      normalize(bjet1_pt);              
      normalize(bjet2_sublpt);            
      normalize(bjet2_lpt);               
      normalize(bjet1_pt_all_cuts);       
      normalize(bjet2_sublpt_all_cuts);   
      normalize(bjet2_mass_all_cuts);      
      normalize(missing_Et);              
      normalize(missing_Et_phi);          
      normalize(missing_Et_eta);          
      normalize(vEt_miss);                
      normalize(vEt_miss_phi);            
      normalize(vEt_miss_eta);            
      normalize(track_pt);                
      normalize(dh_m);                    
      normalize(dh_m_bb);                 
      normalize(dh_m_bb_jets);              
      normalize(dh_m_fat_jet);            
      normalize(dh_eta);                  
      normalize(dh_phi);                  
      normalize(dh_m_MET500);             
      normalize(dh_eta_MET500);           
      normalize(dh_phi_MET500);           
      normalize(dh_m_bb_MET500);          
      normalize(dh_m_fat_jet_MET500);                                         
      normalize(dh_m_fat_jet_all_cuts);



    }//end of finalize method
    
    
    
    
  };//end of public
  
  
  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MonoH_Truth);

}
