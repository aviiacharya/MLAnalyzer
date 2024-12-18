#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

using std::vector;

const unsigned nJets = 50; //TODO: use cfg level nJets_

vector<float> v_top_goodvertices_;
//vector<float> v_top_jetIsSignal_;        //removed
vector<float> v_top_jet_pt_;
vector<float> v_top_jet_eta_;
vector<float> v_top_jet_phi_;
vector<float> v_top_jet_energy_;
vector<float> v_top_jet_m0_;           //thisJet->mass() Line 315 (rest mass m0)
vector<float> v_top_jet_dR_;
//vector<float> v_top_jet_M_;  //Added Invariant mass (M)
vector<float> v_gen_top_mass_;
vector<float> v_gen_top_pt_;
vector<float> v_gen_top_eta_;
vector<float> v_gen_top_phi_;
vector<float> v_reco_top_mass_;
vector<float> v_reco_top_pt_;
vector<float> v_reco_top_eta_;
vector<float> v_reco_top_phi_;



std::vector<std::vector<int> > seljet_genpart_collid;
std::vector<std::vector<int> > seljet_genpart_pdgid;
std::vector<std::vector<int> > seljet_genpart_charge;

std::vector<std::vector<float> > seljet_genpart_px;
std::vector<std::vector<float> > seljet_genpart_py;
std::vector<std::vector<float> > seljet_genpart_pz;
std::vector<std::vector<float> > seljet_genpart_energy;

std::vector<std::vector<int> > seljet_genpart_status;

std::vector<std::vector<int> > seljet_genpart_motherpdgid;
std::vector<std::vector<int> > seljet_genpart_dau1pdgid;
std::vector<std::vector<int> > seljet_genpart_dau2pdgid;

vector<vector<float>> v_top_jetPFCandE_;
vector<vector<float>> v_top_jetPFCandPx_;
vector<vector<float>> v_top_jetPFCandPy_;
vector<vector<float>> v_top_jetPFCandPz_;
vector<vector<int>> v_top_jetPFCandType_;

vector<vector<float>> v_top_jetSV_PtRel_;
vector<vector<float>> v_top_jetSV_ERel_;
vector<vector<float>> v_top_jetSV_PhiRel_;
vector<vector<float>> v_top_jetSV_EtaRel_;
vector<vector<float>> v_top_jetSV_DeltaR_;
vector<vector<float>> v_top_jetSV_Pt_;
vector<vector<float>> v_top_jetSV_Eta_;
vector<vector<float>> v_top_jetSV_Phi_;
vector<vector<float>> v_top_jetSV_Mass_;

vector<vector<float>> v_top_jetSV_ntracks_;
vector<vector<float>> v_top_jetSV_chi2_;
vector<vector<float>> v_top_jetSV_ndf_;
vector<vector<float>> v_top_jetSV_normchi2_;

vector<vector<float>> v_top_jetSV_dxy_;
vector<vector<float>> v_top_jetSV_dxyerr_;
vector<vector<float>> v_top_jetSV_dxysig_;

vector<vector<float>> v_top_jetSV_d3d_;
vector<vector<float>> v_top_jetSV_d3derr_;
vector<vector<float>> v_top_jetSV_d3dsig_;
vector<vector<float>> v_top_jetSV_costhetasvpv_;

// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet_dijet_top ( TTree* tree, edm::Service<TFileService> &fs ) {
    
    tree->Branch("goodvertices",  &v_top_goodvertices_);
    //tree->Branch("jet_IsSignal",     &v_top_jetIsSignal_);
    tree->Branch("jet_Pt",        &v_top_jet_pt_);
    tree->Branch("jet_Eta",       &v_top_jet_eta_);
    tree->Branch("jet_Phi",       &v_top_jet_phi_);   // Just added
    tree->Branch("jet_Energy",    &v_top_jet_energy_);
    tree->Branch("jet_m0",        &v_top_jet_m0_);
   // tree->Branch("jet_M",         &v_top_jet_M_);     //Added tot mass of daughter particles -> Invariant mass
    tree->Branch("jet_dR",        &v_top_jet_dR_);
    
    tree->Branch("gen_top_mass",      &v_gen_top_mass_);
    tree->Branch("gen_top_pt",        &v_gen_top_pt_);
    tree->Branch("gen_top_eta",       &v_gen_top_eta_);
    tree->Branch("gen_top_phi",      &v_gen_top_phi_);

    tree->Branch("reco_top_mass",      &v_reco_top_mass_);
    tree->Branch("reco_top_pt",        &v_reco_top_pt_);
    tree->Branch("reco_top_eta",       &v_reco_top_eta_);
    tree->Branch("reco_top_phi",      &v_reco_top_phi_);
    
    tree->Branch("seljet_genpart_collid", &seljet_genpart_collid);
    tree->Branch("seljet_genpart_pdgid", &seljet_genpart_pdgid);
    tree->Branch("seljet_genpart_charge", &seljet_genpart_charge);
    
    tree->Branch("seljet_genpart_px", &seljet_genpart_px);
    tree->Branch("seljet_genpart_py", &seljet_genpart_py);
    tree->Branch("seljet_genpart_pz", &seljet_genpart_pz);
    tree->Branch("seljet_genpart_energy", &seljet_genpart_energy);
    
    tree->Branch("seljet_genpart_status", &seljet_genpart_status);
    
    tree->Branch("seljet_genpart_motherpdgid", &seljet_genpart_motherpdgid);
    tree->Branch("seljet_genpart_dau1pdgid", &seljet_genpart_dau1pdgid);
    tree->Branch("seljet_genpart_dau2pdgid", &seljet_genpart_dau2pdgid);
    
    tree->Branch("jetPFCandE",        &v_top_jetPFCandE_);
    tree->Branch("jetPFCandPx",       &v_top_jetPFCandPx_);
    tree->Branch("jetPFCandPy",       &v_top_jetPFCandPy_);
    tree->Branch("jetPFCandPz",       &v_top_jetPFCandPz_);
    tree->Branch("jetPFCandType",     &v_top_jetPFCandType_);
    
    tree->Branch("jetSV_PtRel",     &v_top_jetSV_PtRel_);
    tree->Branch("jetSV_ERel",      &v_top_jetSV_ERel_);
    tree->Branch("jetSV_PhiRel",    &v_top_jetSV_PhiRel_);
    tree->Branch("jetSV_EtaRel",  &v_top_jetSV_EtaRel_);
    tree->Branch("jetSV_DeltaR",  &v_top_jetSV_DeltaR_);
    tree->Branch("jetSV_Pt",      &v_top_jetSV_Pt_);
    tree->Branch("jetSV_Eta",     &v_top_jetSV_Eta_);
    tree->Branch("jetSV_Phi",     &v_top_jetSV_Phi_);
    tree->Branch("jetSV_Mass",    &v_top_jetSV_Mass_);
    
    tree->Branch("jetSV_ntracks",   &v_top_jetSV_ntracks_);
    tree->Branch("jetSV_chi2",      &v_top_jetSV_chi2_);
    tree->Branch("jetSV_ndf",       &v_top_jetSV_ndf_);
    tree->Branch("jetSV_normchi2",  &v_top_jetSV_normchi2_);
    
    tree->Branch("jetSV_dxy",       &v_top_jetSV_dxy_);
    tree->Branch("jetSV_dxyerr",    &v_top_jetSV_dxyerr_);
    tree->Branch("jetSV_dxysig",    &v_top_jetSV_dxysig_);
    
    tree->Branch("jetSV_d3d",       &v_top_jetSV_d3d_);
    tree->Branch("jetSV_d3derr",    &v_top_jetSV_d3derr_);
    tree->Branch("jetSV_d3dsig",    &v_top_jetSV_d3dsig_);
    tree->Branch("jetSV_costhetasvpv",  &v_top_jetSV_costhetasvpv_);
} // branchesEvtSel_jet_dijet_top()



// Run jet selection _____________________________________________________//
bool RecHitAnalyzer::runEvtSel_jet_dijet_top( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
    
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByToken( genParticleCollectionT_, genParticles );
    
    edm::Handle<reco::PFJetCollection> jets;
    iEvent.getByToken(jetCollectionT_, jets);
    
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vertexCollectionT_, vertices);
    
    vJetIdxs.clear();
    
    int nJet = 0;
    
    std::vector<TLorentzVector> gen_tops,bdau,wdau,reco_tops,reco_W;
    if (isBoostedTop_) { //if its boosted top sample
        for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
            
            int id = iGen->pdgId();
            std::cout << "Gen Particle id: " << id << std::endl;
            if(abs(id) != 6 || iGen->numberOfDaughters()!=2) continue;
            int iw=-1;
            int ib=-1;
            if (abs(iGen->daughter(0)->pdgId())==24 && abs(iGen->daughter(1)->pdgId())==5)
            {
                iw=0;ib=1;
            }
            else
            {
                if(abs(iGen->daughter(1)->pdgId())==24 && abs(iGen->daughter(0)->pdgId())==5)
                {
                    iw=1;ib=0;
                }
                else continue;
            }
            const reco::Candidate *d = iGen->daughter(iw);
            const reco::Candidate *b = iGen->daughter(ib);
            while(d->numberOfDaughters() == 1) d = d->daughter(0);
            if(!(abs(d->daughter(0)->pdgId()) > 6 || abs(d->daughter(1)->pdgId()) > 6)) continue; //Initially 10,PDGID should be less tha 7 for quarks
            
            const reco::Candidate* W_dau1 = d->daughter(0);   //daughters of W
            const reco::Candidate* W_dau2 = d->daughter(1);
            
            TLorentzVector the_top,the_w,the_b;   //Gen level
            the_top.SetPtEtaPhiE(iGen->pt(),iGen->eta(),iGen->phi(),iGen->energy());
            the_w.SetPtEtaPhiE(d->pt(),d->eta(),d->phi(),d->energy());
            the_b.SetPtEtaPhiE(b->pt(),b->eta(),b->phi(),b->energy());
            
            gen_tops.push_back(the_top);
            wdau.push_back(the_w);
            bdau.push_back(the_b);
            
            TLorentzVector vec_W_dau1, vec_W_dau2, reco_w, reco_top;  //Reco level
            vec_W_dau1.SetPtEtaPhiE(W_dau1->pt(), W_dau1->eta(), W_dau1->phi(), W_dau1->energy());
            vec_W_dau2.SetPtEtaPhiE(W_dau2->pt(), W_dau2->eta(), W_dau2->phi(), W_dau2->energy());
            reco_w = vec_W_dau1 + vec_W_dau2;
            
            reco_top = reco_w + the_b;
            
            reco_tops.push_back(reco_top);
            reco_W.push_back(reco_w);
            
            
        } //gen particle loop

        v_gen_top_mass_.clear();
        v_gen_top_pt_.clear();
    	v_gen_top_eta_.clear(); 
        v_gen_top_phi_.clear();
        v_reco_top_mass_.clear();
        v_reco_top_pt_.clear();
        v_reco_top_eta_.clear();
        v_reco_top_phi_.clear();
        //add vectors................
        
        // Loop over jets
        std::cout << "Number of tops in event: " << gen_tops.size() << std::endl;
        for ( unsigned ihad=0; ihad<gen_tops.size(); ihad++)
        {
            TLorentzVector the_top = gen_tops[ihad];
	    TLorentzVector reco_top = reco_tops[ihad];

            // Calculate properties
            float gen_top_mass = the_top.M();
            float gen_top_pt = the_top.Pt();
            float gen_top_eta = the_top.Eta();
            float gen_top_phi = the_top.Phi();
    
            v_gen_top_mass_.push_back(gen_top_mass);
            v_gen_top_pt_.push_back(gen_top_pt);
            v_gen_top_eta_.push_back(gen_top_eta);
            v_gen_top_phi_.push_back(gen_top_phi);
            
            float reco_top_mass = reco_top.M();
            float reco_top_pt = reco_top.Pt();
            float reco_top_eta = reco_top.Eta();
            float reco_top_phi = reco_top.Phi();
            
            v_reco_top_mass_.push_back(reco_top_mass);
            v_reco_top_pt_.push_back(reco_top_pt);
            v_reco_top_eta_.push_back(reco_top_eta);
            v_reco_top_phi_.push_back(reco_top_phi);
            

            float jet_P = 0;
            float jet_E = 0;
            float top_invM = 0;
            float minJetPt_ = 0.0;
            float maxJetEta_ = 2.4;
           // std::cout << "Number of jets in event: " << jets->size() << std::endl;
            
            for ( unsigned iJ = 0; iJ != jets->size(); ++iJ )
            {
                
                reco::PFJetRef iJet( jets, iJ );
                TLorentzVector vjet;
                vjet.SetPtEtaPhiE(iJet->pt(),iJet->eta(),iJet->phi(),iJet->energy());
                
                jet_P = jet_P + iJet->pt();
                jet_E = jet_E + iJet->energy();
                if ( std::abs(iJet->pt()) < minJetPt_ ) continue;
                if ( std::abs(iJet->eta()) > maxJetEta_) continue;
              /*  if ( std::abs(iJet->pt()) < minJetPt_ )
 
               	{ std::cout << " >> iJte_pt < minJetPt:"  << std::endl; continue;
                }

                if ( std::abs(iJet->eta()) > maxJetEta_) 

                { std::cout << " >> iJte_eta < maxJetEta:"  << std::endl; 
                  std::cout << " >> iJte_eta: "<< iJet->eta()  << std::endl; continue;
                }*/
                if ( debug ) std::cout << " >> jet[" << iJ << "]Pt:" << iJet->pt() << " jetE:" << iJet->energy() << " jet_m0:" << iJet->mass() << std::endl;
               // if (gen_tops[ihad].DeltaR(vjet)>0.8) continue;
               // if (wdau[ihad].DeltaR(vjet)>0.8) continue;
               // if (bdau[ihad].DeltaR(vjet)>0.8) continue;
                
                if ( debug ) std::cout << " >> jet[" << iJ << "]Pt:" << iJet->pt() << " jetE:" << iJet->energy() << " jet_m0:" << iJet->mass() << std::endl;
                
                //top_invM = sqrt(jet_E*jet_E - jet_P*jet_P);
                
                vJetIdxs.push_back(iJ);
                
                nJet++;
                break;// This should allow two hardonic tops
            } // jets
            
            
            if ( (nJets_ > 0) && (nJet >= nJets_) ) break;
        } // hadronic tops
    } // isBoostedTop
    
    if ( debug ) {
        for(int thisJetIdx : vJetIdxs)
            std::cout << " >> vJetIdxs:" << thisJetIdx << std::endl;
    }
    
    if ( (nJets_ > 0) && (nJet != nJets_) ){
        if ( debug ) std::cout << " Fail jet multiplicity:  " << nJet << " < " << nJets_ << std::endl;
        return false;
    }
    
    if ( vJetIdxs.size() == 0){
        if ( debug ) std::cout << " No passing jets...  " << std::endl;
        return false;
    }
    
    if ( debug ) std::cout << " >> has_jet_dijet: passed" << std::endl;
    return true;
    
    unsigned int nMatchedJets = 0;
    
    
    
} // runEvtSel_jet_dijet_top()

// Fill branches and histograms _____________________________________________________//
void RecHitAnalyzer::fillEvtSel_jet_dijet_top ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {
    
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vertexCollectionT_, vertices);
    
    edm::Handle<reco::PFJetCollection> jets;
    iEvent.getByToken(jetCollectionT_, jets);
    
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByToken( genParticleCollectionT_, genParticles );
    
    edm::Handle<reco::VertexCompositePtrCandidateCollection> secVertices;
    iEvent.getByToken(secVertexCollectionT_, secVertices);
    
    v_top_goodvertices_.clear();
   // v_top_jetIsSignal_.clear();
    v_top_jet_pt_.clear();
    v_top_jet_eta_.clear();
    v_top_jet_phi_.clear();
    v_top_jet_energy_.clear();
    v_top_jet_m0_.clear();
    v_top_jet_dR_.clear();
   // v_top_jet_M_.clear();
   // v_gen_top_mass_.clear();
   // v_gen_top_pt_.clear();
   // v_gen_top_eta_.clear();
    
    
    seljet_genpart_collid.clear();
    seljet_genpart_pdgid.clear();
    seljet_genpart_charge.clear();
    
    seljet_genpart_px.clear();
    seljet_genpart_py.clear();
    seljet_genpart_pz.clear();
    seljet_genpart_energy.clear();
    
    seljet_genpart_status.clear();
    
    seljet_genpart_motherpdgid.clear();
    seljet_genpart_dau1pdgid.clear();
    seljet_genpart_dau2pdgid.clear();
    
    v_top_jetPFCandE_.clear();
    v_top_jetPFCandPx_.clear();
    v_top_jetPFCandPy_.clear();
    v_top_jetPFCandPz_.clear();
    v_top_jetPFCandType_.clear();
    
    v_top_jetSV_PtRel_.clear();
    v_top_jetSV_ERel_.clear();
    v_top_jetSV_PhiRel_.clear();
    v_top_jetSV_EtaRel_.clear();
    v_top_jetSV_DeltaR_.clear();
    v_top_jetSV_Pt_.clear();
    v_top_jetSV_Eta_.clear();
    v_top_jetSV_Phi_.clear();
    v_top_jetSV_Mass_.clear();
    
    v_top_jetSV_ntracks_.clear();
    v_top_jetSV_chi2_.clear();
    v_top_jetSV_ndf_.clear();
    v_top_jetSV_normchi2_.clear();
    
    v_top_jetSV_dxy_.clear();
    v_top_jetSV_dxyerr_.clear();
    v_top_jetSV_dxysig_.clear();
    
    v_top_jetSV_d3d_.clear();
    v_top_jetSV_d3derr_.clear();
    v_top_jetSV_d3dsig_.clear();
    v_top_jetSV_costhetasvpv_.clear();
    
    unsigned int goodVertices = 0;
    
    if ( debug ) std::cout << " >>>>>>>>>>>>>>>>>>>> evt:" << std::endl;
    
    if (vertices.isValid())
        if (vertices->size() > 0)
            for (auto v : *vertices)
                if (v.ndof() >= 4 && !v.isFake())
                    ++goodVertices;
    if ( debug ) std::cout << "\t" << " good vertices in the event (PU) = " << goodVertices << std::endl;
    
    //if ( debug ) std::cout << "\t" << " JETS IN THE EVENT = " << vTauJets.size() << " | Selection requires minpT = " << minJetPt_ << " and maxEta = "<< maxJetEta_ << std::endl;
    
    v_top_goodvertices_.push_back(goodVertices);
    //v_top_jet_M_.push_back( thisJet->top_invM() );
    
    
    for (int iJ : vJetIdxs){
        
        reco::PFJetRef thisJet( jets, iJ );
        if ( debug ) std::cout << " >> Jet[" << iJ << "] Pt:" << thisJet->pt() << std::endl;
        
        //v_top_jetIsSignal_.push_back( 0 );
        v_top_jet_pt_.push_back( thisJet->pt() );
        v_top_jet_eta_.push_back( thisJet->eta() );
        v_top_jet_energy_.push_back( thisJet->energy() );
        v_top_jet_m0_.push_back( thisJet->mass() );
        v_top_jet_phi_.push_back( thisJet->phi() );
        //v_gen_top_mass_.push_back( the_top->mass() );
        //v_gen_top_pt_.push_back( the_top->pt() );
        //v_gen_top_eta_.push_back( the_top->eta() );
        
        
        TLorentzVector TLVJet(thisJet->px(),thisJet->py(),thisJet->pz(),thisJet->energy());
        double cosTheta = TLVJet.CosTheta();
        if (cosTheta*cosTheta >=0)
            TLVJet.SetPx(0.0001);
        
        std::vector<int> genpart_collid;
        std::vector<int> genpart_pdgid;
        std::vector<int> genpart_charge;
        
        std::vector<float> genpart_px;
        std::vector<float> genpart_py;
        std::vector<float> genpart_pz;
        std::vector<float> genpart_energy;
        
        std::vector<int> genpart_status;
        
        std::vector<int> genpart_motherpdgid;
        std::vector<int> genpart_dau1pdgid;
        std::vector<int> genpart_dau2pdgid;
        
        std::vector<reco::GenParticle>::const_iterator genpartIterator      = (genParticles.product())->begin();
        std::vector<reco::GenParticle>::const_iterator genpartIteratorEnd   = (genParticles.product())->end();
        for ( ; genpartIterator != genpartIteratorEnd; genpartIterator++)
        {
            
            TLorentzVector TLVgenpart(genpartIterator->px(),genpartIterator->py(),genpartIterator->pz(),genpartIterator->energy());
            cosTheta = TLVgenpart.CosTheta();
            if (cosTheta*cosTheta >=0)
                TLVgenpart.SetPx(0.0001);
            
            
            if (TLVJet.DeltaR(TLVgenpart)<0.8)
            {
                float dR = TLVJet.DeltaR(TLVgenpart);
                v_top_jet_dR_.push_back(dR);
                genpart_collid.push_back(genpartIterator->collisionId());
                genpart_pdgid.push_back(genpartIterator->pdgId());
                genpart_charge.push_back(genpartIterator->charge());
                
                genpart_px.push_back(genpartIterator->px());
                genpart_py.push_back(genpartIterator->py());
                genpart_pz.push_back(genpartIterator->pz());
                genpart_energy.push_back(genpartIterator->energy());
                
                genpart_status.push_back(genpartIterator->status());
                
                if (genpartIterator->numberOfMothers()>0)
                {
                    genpart_motherpdgid.push_back(genpartIterator->mother(0)->pdgId());
                }
                else
                {
                    genpart_motherpdgid.push_back(-9999);
                }
                
                switch (genpartIterator->numberOfDaughters())
                {
                    case 0:
                        genpart_dau1pdgid.push_back(-9999);
                        genpart_dau2pdgid.push_back(-9999);
                        break;
                        
                    case 1:
                        genpart_dau1pdgid.push_back(genpartIterator->daughter(0)->pdgId());
                        genpart_dau2pdgid.push_back(-9999);
                        break;
                        
                    default:
                        genpart_dau1pdgid.push_back(genpartIterator->daughter(0)->pdgId());
                        genpart_dau2pdgid.push_back(genpartIterator->daughter(1)->pdgId());
                        break;
                }//switch
            }//DeltaR condition
        }//genpart loop
        seljet_genpart_collid.push_back(genpart_collid);
        seljet_genpart_pdgid.push_back(genpart_pdgid);
        seljet_genpart_charge.push_back(genpart_charge);
        
        seljet_genpart_px.push_back(genpart_px);
        seljet_genpart_py.push_back(genpart_py);
        seljet_genpart_pz.push_back(genpart_pz);
        seljet_genpart_energy.push_back(genpart_energy);
        
        seljet_genpart_status.push_back(genpart_status);
        
        seljet_genpart_motherpdgid.push_back(genpart_motherpdgid);
        seljet_genpart_dau1pdgid.push_back(genpart_dau1pdgid);
        seljet_genpart_dau2pdgid.push_back(genpart_dau2pdgid);
        
        // jet constituents
        vector<float> pfcand_px;
        vector<float> pfcand_py;
        vector<float> pfcand_pz;
        vector<float> pfcand_energy;
        vector<int> pfcand_type;
        
        //std::vector<reco::PFCandidatePtr> jetConstituents = thisJet->getPFConstituents();
        unsigned int nConstituents = thisJet->getPFConstituents().size();
        std::cout << " >>>>>> Jet " << iJ << " has "<< nConstituents << " constituents" << std::endl;
        if ( debug ) std::cout << " >>>>>> Jet " << iJ << " has "<< nConstituents << " constituents" << std::endl;
        for ( unsigned int j = 0; j < nConstituents; j++ ) {
            const reco::PFCandidatePtr jetPFCand = thisJet->getPFConstituent( j );
            //if ( debug ) std::cout << " >>>>>>>>>>>  Jet constituent " << j << "-> E:" << jetPFCand->energy() << " px:" << jetPFCand->px() << " py:" << jetPFCand->py() << " pz:" << jetPFCand->pz() << std::endl;
            //std::cout << " >>>>>>>>>>>  Jet constituent " << j << "-> E:" << jetPFCand->energy() << " px:" << jetPFCand->px() << " py:" << jetPFCand->py() << " pz:" << jetPFCand->pz() << std::endl;
            
            pfcand_px.push_back(jetPFCand->px());
            pfcand_py.push_back(jetPFCand->py());
            pfcand_pz.push_back(jetPFCand->pz());
            pfcand_energy.push_back(jetPFCand->energy());
            pfcand_type.push_back((int)jetPFCand->particleId());
        }//jet constituents loop
        iJ++;
        
        v_top_jetPFCandE_.push_back( pfcand_energy );
        v_top_jetPFCandPx_.push_back( pfcand_px );
        v_top_jetPFCandPy_.push_back( pfcand_py );
        v_top_jetPFCandPz_.push_back( pfcand_pz );
        v_top_jetPFCandType_.push_back( pfcand_type );
        
        std::vector<const reco::VertexCompositePtrCandidate*> jetSVs;
        for (const auto &sv : *secVertices){
            if (reco::deltaR(sv, *thisJet) < 0.8) { //0.4 before
                jetSVs.push_back(&sv);
            }
        }
        
        // sort by dxy significance
        const auto &pv = vertices->at(0);
        std::sort(jetSVs.begin(), jetSVs.end(), [&](const reco::VertexCompositePtrCandidate *sv1, const reco::VertexCompositePtrCandidate *sv2){
            return vertexDxy(*sv1, pv).significance() > vertexDxy(*sv2, pv).significance();
        });
        
        
        vector<float> sv_PtRel;
        vector<float> sv_ERel;
        vector<float> sv_PhiRel;
        vector<float> sv_EtaRel;
        vector<float> sv_DeltaR;
        vector<float> sv_Pt;
        vector<float> sv_Eta;
        vector<float> sv_Phi;
        vector<float> sv_Mass;
        
        vector<float> sv_ntracks;
        vector<float> sv_chi2;
        vector<float> sv_ndf;
        vector<float> sv_normchi2;
        
        vector<float> sv_dxy;
        vector<float> sv_dxyerr;
        vector<float> sv_dxysig;
        
        vector<float> sv_d3d;
        vector<float> sv_d3derr;
        vector<float> sv_d3dsig;
        vector<float> sv_costhetasvpv;
        
        
        float etasign = thisJet->eta()>0 ? 1 : -1;
        
        for (const auto *sv : jetSVs){
            
            std::cout << "=================================================Secondary vertex pt:" << sv->pt() << " eta:" << sv->eta() << "  phi:" << sv->phi() << std::endl;
            
            sv_PtRel.push_back(sv->pt()/thisJet->pt());
            sv_ERel.push_back(sv->energy()/thisJet->energy());
            sv_PhiRel.push_back(reco::deltaPhi(*sv, *thisJet));
            sv_EtaRel.push_back(etasign * (sv->eta() - thisJet->eta()));
            sv_DeltaR.push_back(reco::deltaR(*sv, *thisJet));
            sv_Pt.push_back(sv->pt());
            sv_Eta.push_back(sv->eta());
            sv_Phi.push_back(sv->phi());
            sv_Mass.push_back(sv->mass());
            
            //sv properties
            sv_ntracks.push_back(sv->numberOfDaughters());
            sv_chi2.push_back(sv->vertexChi2());
            sv_ndf.push_back(sv->vertexNdof());
            sv_normchi2.push_back(catchInfs(sv->vertexNormalizedChi2()));
            
            
            const auto &dxy = vertexDxy(*sv, pv);
            sv_dxy.push_back(dxy.value());
            sv_dxyerr.push_back(dxy.error());
            sv_dxysig.push_back(dxy.significance());
            
            
            const auto &d3d = vertexD3d(*sv, pv);
            sv_d3d.push_back(d3d.value());
            sv_d3derr.push_back(d3d.error());
            sv_d3dsig.push_back(d3d.significance());
            sv_costhetasvpv.push_back(vertexDdotP(*sv, pv));
            
        }
        
        v_top_jetSV_PtRel_.push_back(sv_PtRel);
        v_top_jetSV_ERel_.push_back(sv_ERel);
        v_top_jetSV_PhiRel_.push_back(sv_PhiRel);
        v_top_jetSV_EtaRel_.push_back(sv_EtaRel);
        v_top_jetSV_DeltaR_.push_back(sv_DeltaR);
        v_top_jetSV_Pt_.push_back(sv_Pt);
        v_top_jetSV_Eta_.push_back(sv_Eta);
        v_top_jetSV_Phi_.push_back(sv_Phi);
        v_top_jetSV_Mass_.push_back(sv_Mass);
        
        v_top_jetSV_ntracks_.push_back(sv_ntracks);
        v_top_jetSV_chi2_.push_back(sv_chi2);
        v_top_jetSV_ndf_.push_back(sv_ndf);
        v_top_jetSV_normchi2_.push_back(sv_normchi2);
        
        v_top_jetSV_dxy_.push_back(sv_dxy);
        v_top_jetSV_dxyerr_.push_back(sv_dxyerr);
        v_top_jetSV_dxysig_.push_back(sv_dxysig);
        
        v_top_jetSV_d3d_.push_back(sv_d3d);
        v_top_jetSV_d3derr_.push_back(sv_d3derr);
        v_top_jetSV_d3dsig_.push_back(sv_d3dsig);
        v_top_jetSV_costhetasvpv_.push_back(sv_costhetasvpv);
        
        
    }//jet loop
    
} // fillEvtSel_jet_dthisJet_top()
