/*************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

////////////////////////////////////////////////////////////////////////////
//
//  This class is used to reconstruct the neutral Xi(1530) resonance for
//  Run2 data.
//  This class essentially combines charged Xi candidates from the Xi Vert-
//  exer with primary charged pions.
//
//  author: Bong-Hwi Lim (bong-hwi.lim@cern.ch)
//  Origin of the Code: Beomkyu KIM (kimb@cern.ch)
//
//  Last Modified Date: 2018/09/25
//
////////////////////////////////////////////////////////////////////////////


#include "TFile.h"
#include "TChain.h"
#include "TSystem.h"
#include "THistManager.h"
#include "AliAnalysisTaskXi1530.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliGenEventHeader.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliMultSelection.h"
#include "AliESDcascade.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliMultiplicity.h"

// Some constants
const Double_t        pi = TMath::Pi();
const Double_t  pionmass = AliPID::ParticleMass(AliPID::kPion);
//const Double_t  Ximass = AliPID::ParticleMass(AliPID::kCascade);
const Double_t  Ximass = 1.32171;
enum { kData=1, kLS, kMixing, kAllType}; //P=Positive charge, N=Negative


AliAnalysisTaskXi1530RunTable::AliAnalysisTaskXi1530RunTable() :
fCollisionType(kUnknownCollType)
{;}

AliAnalysisTaskXi1530RunTable::AliAnalysisTaskXi1530RunTable(Int_t runnumber)
{
    // Need to be fixed
    if (runnumber>=256504 && runnumber<=260014) fCollisionType=kPP;//LHC16kl
    else fCollisionType=kUnknownCollType;
}
AliAnalysisTaskXi1530RunTable::~AliAnalysisTaskXi1530RunTable()
{;}

//___________________________________________________________________
AliAnalysisTaskXi1530::AliAnalysisTaskXi1530()
:AliAnalysisTaskSE("AliAnalysisTaskXi1530"),
    fOption(),
    goodtrackindices(),
    goodcascadeindices(),
    fEMpool ()
{
    DefineInput (0, TChain::Class());
    DefineOutput (1, TList::Class());
}
//___________________________________________________________________
AliAnalysisTaskXi1530::AliAnalysisTaskXi1530(
    const char *name,
    const char *option
)
    :AliAnalysisTaskSE(name),
    fOption(option),
    goodtrackindices(),
    goodcascadeindices(),
    fEMpool ()
{
    DefineInput (0, TChain::Class());
    DefineOutput (1, TList::Class());
}
AliAnalysisTaskXi1530::AliAnalysisTaskXi1530
(
 const AliAnalysisTaskXi1530& ap
 )
: fOption(ap.fOption)
, goodtrackindices(ap.goodtrackindices)
, goodcascadeindices(ap.goodcascadeindices)
, fEMpool (ap.fEMpool)
{
}
//___________________________________________________________________
AliAnalysisTaskXi1530& AliAnalysisTaskXi1530::operator =
(
 const AliAnalysisTaskXi1530& ap
 )
{
    // assignment operator
    
    this->~AliAnalysisTaskXi1530();
    new(this) AliAnalysisTaskXi1530(ap);
    return *this;
}
//___________________________________________________________________
AliAnalysisTaskXi1530::~AliAnalysisTaskXi1530()
{
    delete fOutput;
    delete fTrackCuts;
    delete fPIDResponse;
    delete fRunTable;
}
//________________________________________________________________________
void AliAnalysisTaskXi1530::UserCreateOutputObjects()
{
    std::cout << "User Create Object" << std::endl;
    // Histograms container
    fOutput = new TList();
    fOutput->SetOwner(kTRUE);
    
    // TrackCuts for Xi1530--------------------------------------------------
    fTrackCuts = new AliESDtrackCuts();
    fTrackCuts -> GetStandardITSTPCTrackCuts2011(kTRUE,kTRUE);
    fTrackCuts -> SetEtaRange(-0.8,0.8);
    fTrackCuts -> SetPtRange(0.15, 1e20);
    // ----------------------------------------------------------------------
    
    fHistos = new THistManager("Xi1530hists");
    
    auto binType = AxisStr("Type",{"DATA","Mixing"});
    if (IsAA) binCent = AxisFix("Cent",10,0,100);
    //else binCent = AxisFix("Cent",300,0,300);
    else binCent = AxisVar("Cent",{0,1,5,10,15,20,30,40,50,70,100}); // 0 ~ -1, overflow, 0 ~ +1, underflow
    auto binPt   = AxisFix("Pt",200,0,20);
    auto binMass = AxisFix("Mass",2000,0.5,2.5);
    
    CreateTHnSparse("hInvMass","InvMass",4,{binType,binCent,binPt,binMass},"s");
    CreateTHnSparse("hInvMass_dXi","InvMass",3,{binCent,binPt,binMass},"s");
    CreateTHnSparse("hMult","Multiplicity",1,{binCent},"s");
    
    vector<TString> ent = {"All","PS","PSpileup","Goodz","Goodzcut"};
    auto hNofEvt = fHistos->CreateTH1("hEventNumbers","",ent.size(), 0, ent.size());
    for(auto i=0u;i<ent.size();i++) hNofEvt->GetXaxis()->SetBinLabel(i+1,ent.at(i).Data());
    
    fHistos -> CreateTH2("hPhiEta","",180,0,2*pi,40,-2,2);
    fHistos -> CreateTH2("hPhiEta_Xi","",180,0,2*pi,40,-2,2);
    
    binZ = AxisVar("Z",{-10,-5,-3,-1,1,3,5,10});
    fHistos->CreateTH1("hMul","",200,0,200,"s");
    fHistos->CreateTH1("hZvtx","",600,-30,30,"s");
    
    // MC true inv mass distribution
    CreateTHnSparse("hInvMassMCXi1530","InvMass",3,{binCent,binPt,binMass},"s");
    
    // To get Trigger efficiency in each trk/V0M Multiplicity region
    
    vector<TString> mcent = {"All","IsINELg0","tracklet in |Eta|<1","CINT7 triggered","AliMultiSelection"};
    vector<Int_t> ntrklet = {0,5,10,15,20,25,30,35,40,50};  // # of Tracklet bin
    vector<Int_t> nmult = {0,1,5,10,15,20,30,40,50,70,100}; // V0M Multiplicity bin
    
    auto hNofEvtMC = fHistos->CreateTH1("htotalEvent_MC","",mcent.size(), 0, mcent.size());
    auto htrkINELg0 = fHistos->CreateTH1("htriggered_INELg0_tracklet","",ntrklet.size(), 0, ntrklet.size());
    auto htrkCINT7 = fHistos->CreateTH1("htriggered_CINT7_tracklet","",ntrklet.size(), 0, ntrklet.size());
    auto hmultCINT7 = fHistos->CreateTH1("htriggered_CINT7_VOM","",nmult.size(), 0, nmult.size());
    /*
    for(auto i=0u;i<mcent.size();i++) hNofEvtMC->GetXaxis()->SetBinLabel(i+1,mcent.at(i).Data());
    htrkINELg0->GetXaxis()->SetBinLabel(1,"0 to Inf (MB)");
    for(auto i=1u;i<ntrklet.size();i++) htrkINELg0->GetXaxis()->SetBinLabel(i+1,Form("%d to %d",ntrklet.at(i),ntrklet.at(i)+5));
    htrkCINT7->GetXaxis()->SetBinLabel(1,"0 to Inf (MB)");
    for(auto i=1u;i<ntrklet.size();i++) htrkCINT7->GetXaxis()->SetBinLabel(i+1,Form("%d to %d",ntrklet.at(i),ntrklet.at(i)+5));
    hmultCINT7->GetXaxis()->SetBinLabel(1,"0 - 100 % (MB)");
    for(auto i=1u;i<nmult.size();i++) hmultCINT7->GetXaxis()->SetBinLabel(i+1,Form("%d - %d%%",nmult.at(i),nmult.at(i+1)));
    */
    CreateTHnSparse("hMult_MC","Multiplicity",1,{binCent},"s");
    CreateTHnSparse("hMult_MC_selected","Multiplicity",1,{binCent},"s");
    
    // QA Histograms--------------------------------------------------
    // T P C   P I D
    // before
    fHistos -> CreateTH2("hTPCPIDLambdaProton","",1000,0,20,1000,-5,5);
    fHistos -> CreateTH2("hTPCPIDLambdaPion","",1000,0,20,1000,-5,5);
    fHistos -> CreateTH2("hTPCPIDBachelorPion","",1000,0,20,1000,-5,5);
    fHistos -> CreateTH2("hTPCPIDXi1530Pion","",1000,0,20,1000,-5,5);
    // after
    fHistos -> CreateTH2("hTPCPIDLambdaProton_cut","",1000,0,20,1000,-5,5);
    fHistos -> CreateTH2("hTPCPIDLambdaPion_cut","",1000,0,20,1000,-5,5);
    fHistos -> CreateTH2("hTPCPIDBachelorPion_cut","",1000,0,20,1000,-5,5);
    fHistos -> CreateTH2("hTPCPIDXi1530Pion_cut","",1000,0,20,1000,-5,5);
    
    // D C A
    // between daughters
    // before
    fHistos->CreateTH1("hDCADist_Lambda_BTW_Daughters","",300,0,3,"s");
    fHistos->CreateTH1("hDCADist_Xi_BTW_Daughters","",300,0,3,"s");
    // after
    fHistos->CreateTH1("hDCADist_Lambda_BTW_Daughters_cut","",300,0,3,"s");
    fHistos->CreateTH1("hDCADist_Xi_BTW_Daughters_cut","",300,0,3,"s");
    // to PV
    // before
    fHistos->CreateTH1("fDCADist_lambda_to_PV","",500,0,0.5,"s");
    fHistos->CreateTH1("hDCADist_Xi_to_PV","",500,0,0.5,"s");
    fHistos->CreateTH1("hDCADist_LambdaProton_to_PV","",500,0,0.5,"s");
    fHistos->CreateTH1("hDCADist_LambdaPion_to_PV","",500,0,0.5,"s");
    fHistos->CreateTH1("hDCADist_BachelorPion_to_PV","",500,0,0.5,"s");
    // after
    fHistos->CreateTH1("fDCADist_lambda_to_PV_cut","",500,0,0.5,"s");
    fHistos->CreateTH1("hDCADist_Xi_to_PV_cut","",500,0,0.5,"s");
    fHistos->CreateTH1("hDCADist_LambdaProton_to_PV_cut","",500,0,0.5,"s");
    fHistos->CreateTH1("hDCADist_LambdaPion_to_PV_cut","",500,0,0.5,"s");
    fHistos->CreateTH1("hDCADist_BachelorPion_to_PV_cut","",500,0,0.5,"s");
    
    // C P A
    // before
    fHistos->CreateTH1("hCosPA_lambda","",150,0.85,1.0,"s");
    fHistos->CreateTH1("hCosPA_Xi","",150,0.85,1.0,"s");
    // after
    fHistos->CreateTH1("hCosPA_lambda_cut","",150,0.85,1.0,"s");
    fHistos->CreateTH1("hCosPA_Xi_cut","",150,0.85,1.0,"s");
    
    // M a s s   W i n d o w
    // before
    fHistos->CreateTH1("hMass_Xi","",200,1.2,1.4,"s");
    // after
    fHistos->CreateTH1("hMass_Xi_cut","",200,1.2,1.4,"s");
    
    fEMpool.resize(binCent.GetNbins(),vector<eventpool> (binZ.GetNbins()));
    PostData(1, fHistos->GetListOfHistograms());
}

//________________________________________________________________________
void AliAnalysisTaskXi1530::UserExec(Option_t *)
{
    std::cout << "AliAnalysisTaskXi1530:: UserExec" << std::endl;
    // Pointer to a event----------------------------------------------------
    AliVEvent *event = InputEvent();
    if (!event)
    {
        std::cout << "ERROR: Could not retrieve event" << std::endl;
        return;
    }
    // ----------------------------------------------------------------------
    
    // Connect to ESD tree --------------------------------------------------
    event->IsA()==AliESDEvent::Class()
    ? fEvt = dynamic_cast<AliESDEvent*>(event)
    : fEvt = dynamic_cast<AliAODEvent*>(event);
    if (!fEvt) return;
    
    // ----------------------------------------------------------------------
    
    // Multiplicity(centrality) ---------------------------------------------
    // fCent:
    //       0-100: Selected
    //       999: Not selected
    //       -999: No MultSection
    //
    fCent = -999; // Multiplicity
    AliMultSelection *MultSelection = (AliMultSelection*) fEvt->FindListObject("MultSelection");
    if(MultSelection)
    {
        if (!(MultSelection->IsEventSelected()))
        {
            AliInfo("This event is not selected: AliMultSelection");
            fCent += MultSelection->GetMultiplicityPercentile("V0M");
        }
        fCent = 999;
    }
    else
    {
        //If this happens, re-check if AliMultSelectionTask ran before your task!
        AliInfo("Didn't find MultSelection!");
    }
    // ----------------------------------------------------------------------
    
    // Preparation for MC ---------------------------------------------------
    if (IsMC)
    {
        std::cout << "AliAnalysisTaskXi1530:: MC Check" << std::endl;
        if (fEvt->IsA()==AliESDEvent::Class()){
            //AliMCEvent  *mcEvent = static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->MCEvent();
            AliMCEvent  *mcEvent = MCEvent();
            if(!mcEvent) return;
            fMCStack = (AliStack*) mcEvent->Stack();
        }
        else{
            fMCArray = (TClonesArray*) fEvt->FindListObject("mcparticles");
        }
    }

    // Load InputHandler for each event---------------------------------------
    AliInputEventHandler* inputHandler = (AliInputEventHandler*)
    AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
    // -----------------------------------------------------------------------

    fPIDResponse = (AliPIDResponse*) inputHandler->GetPIDResponse();
    if(!fPIDResponse){
        std::cout << "AliAnalysisTaskXi1530:: No PIDd\n" << std::endl;
    }

    Bool_t IsMinimumBias = kFALSE;
    fHistos -> FillTH1("hEventNumbers","All",1);
    // In Complete DAQ Event Cut----------------------------------------------
    if (fEvt->IsIncompleteDAQ()) {
        std::cout << "Reject: IsIncompleteDAQ" << std::endl;;
        PostData(1, fHistos->GetListOfHistograms());
        return;
    }
    fHistos -> FillTH1("hEventNumbers","PS",1);
    // -----------------------------------------------------------------------
    
    // Vertex Check-----------------------------------------------------------
    //const AliVVertex* pVtx      = fEvt->GetPrimaryVertex() ;
    const AliVVertex* trackVtx  = fEvt->GetPrimaryVertexTPC() ;
    const AliVVertex* spdVtx    = fEvt->GetPrimaryVertexSPD() ;

    Bool_t IsGoodVertex = kFALSE;
    Bool_t IsGoodVertexCut = kFALSE;
    AliVTrack * track1;
    
    if (spdVtx) {
        fZ = spdVtx->GetZ();
        if (spdVtx->GetNContributors()<1) IsGoodVertex = kFALSE;
        else {
            fHistos -> FillTH1("hEventNumbers","Goodz",1);
            IsGoodVertex = kTRUE;
            fHistos->FillTH1("hZvtx",fZ);
        }
    } else IsGoodVertex = kFALSE;
    
    if ( IsGoodVertex && fabs(fZ)<10.) {
        IsGoodVertexCut = kTRUE;
        if (IsMinimumBias) {
            fHistos -> FillTH1("hEventNumbers","Goodzcut",1);
            FillTHnSparse("hCent",{fCent});
        }
    }
    zbin = binZ.FindBin(fZ) -1;
    centbin = binCent.FindBin(fCent) -1;
    
    if (IsGoodVertexCut){
        if (this -> GoodTracksSelection() && this -> GoodCascadeSelection()) this -> FillTracks();
    }
    
    PostData(1, fHistos->GetListOfHistograms());
}
//________________________________________________________________________
Bool_t AliAnalysisTaskXi1530::GoodTracksSelection(){
    // Choose Good Tracks from AliESDtracks and save the label of them, and save them for event mixing
    // it includes several cuts:
    // - TPC PID cut for pion
    // - Eta cut
    //
    const UInt_t ntracks = fEvt ->GetNumberOfTracks();
    goodtrackindices.clear();
    
    AliVTrack * track;
    
    tracklist *etl;
    eventpool *ep;
    //Event mixing pool
    if (fsetmixing){
        ep = &fEMpool[centbin][zbin];
        ep -> push_back( tracklist() );
        etl = &(ep->back());
    }
    fNTracks = 0;
    for (UInt_t it = 0; it<ntracks; it++){
        if (fEvt->IsA()==AliESDEvent::Class()){
            track = (AliESDtrack*) fEvt ->GetTrack(it);
            if (!track) continue;
            if (!fTrackCuts->AcceptTrack((AliESDtrack*) track)) continue;
            //if (!track->IsOn(AliVTrack::kITSpureSA)) continue;
            fHistos->FillTH2("hPhiEta",track->Phi(),track->Eta());
        }
        else {
            track = (AliAODTrack*) fEvt ->GetTrack(it);
            if (!track) continue;
            if( ! ((AliAODTrack*) track)->TestFilterBit(fFilterBit)) continue;
            if (track->Pt()<fptcut) continue;
            if (abs(track->Eta())>fetacut) continue;
            fHistos->FillTH2("hPhiEta",track->Phi(),track->Eta());
        }
        // PID cut for pion
        Double_t fTPCNSigPion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
        fHistos->FillTH2("hTPCPIDXi1530Pion",track->GetTPCmomentum(),fTPCNSigPion);
        
        if (abs(fTPCNSigPion) > 3.) continue;
        fHistos->FillTH2("hTPCPIDXi1530Pion_cut",track->GetTPCmomentum(),fTPCNSigPion);
        
        fNTracks++;
        
        // Eta cut
        if(abs(track->Eta())>0.8) continue;
        
        goodtrackindices.push_back(it);
        
        //Event mixing pool
        if (fsetmixing) etl->push_back( (AliVTrack*)track -> Clone() );
    }
    
    if (fsetmixing){
        if (!goodtrackindices.size()) ep->pop_back();
        if ( ep->size() > 5 ){
            for (auto it: ep->front()) delete it;
            ep->pop_front();
        }
    }
    return goodtrackindices.size();
}

Bool_t AliAnalysisTaskXi1530::GoodCascadeSelection(){
    // Choose Good Cascade from AliESDcascade and save the label of them
    // it includes several cuts:
    // - daughter particle standard track cut
    // - daughter particle PID cut
    // - DCA cuts between Lambda daughters and Xi daughters
    // - PV DCA(Impact parameter) cut for Xi/Lambda/all daughters
    // - Cosine Pointing Angle cut for Xi and Lamnbda
    // - Mass window cut for Xi
    // - Eta cut
    //
    goodcascadeindices.clear();
    const UInt_t ncascade = fEvt->GetNumberOfCascades();
    const AliVVertex* pVtx      = fEvt->GetPrimaryVertex() ;
        Double_t PVx, PVy, PVz;
        PVx = pVtx->GetX();
        PVy = pVtx->GetY();
        PVz = pVtx->GetZ();
    Double_t bField = fEvt->GetMagneticField();
    
    const AliESDcascade *Xicandidate;
    
    fNCascade = 0;
    for (UInt_t it = 0; it<ncascade; it++){
        if (fEvt->IsA()==AliESDEvent::Class()){
            Xicandidate = ((AliESDEvent*)fEvt)->GetCascade(it);
            if (!Xicandidate) continue;
            
            if (TMath::Abs( Xicandidate->GetPindex()) == TMath::Abs( Xicandidate->GetNindex())) continue;
            if (TMath::Abs( Xicandidate->GetPindex()) == TMath::Abs( Xicandidate->GetBindex())) continue;
            if (TMath::Abs( Xicandidate->GetNindex()) == TMath::Abs( Xicandidate->GetBindex())) continue;
            
            AliESDtrack *pTrackXi   = ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs( Xicandidate->GetPindex()));
            AliESDtrack *nTrackXi   = ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs( Xicandidate->GetNindex()));
            AliESDtrack *bTrackXi   = ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs( Xicandidate->GetBindex()));
            
            // Standard track QA cuts
            if (!fTrackCuts->AcceptTrack(pTrackXi)) continue;
            if (!fTrackCuts->AcceptTrack(nTrackXi)) continue;
            if (!fTrackCuts->AcceptTrack(bTrackXi)) continue;
            
            // PID cuts for Xi daughters
            Double_t fTPCNSigProton = fPIDResponse->NumberOfSigmasTPC(nTrackXi, AliPID::kProton);
            Double_t fTPCNSigPion1 = fPIDResponse->NumberOfSigmasTPC(pTrackXi, AliPID::kPion);
            Double_t fTPCNSigPion2 = fPIDResponse->NumberOfSigmasTPC(bTrackXi, AliPID::kPion);
                fHistos->FillTH2("hTPCPIDLambdaProton",pTrackXi->GetTPCmomentum(),fTPCNSigProton);
                fHistos->FillTH2("hTPCPIDLambdaPion",nTrackXi->GetTPCmomentum(),fTPCNSigProton);
                fHistos->FillTH2("hTPCPIDBachelorPion",bTrackXi->GetTPCmomentum(),fTPCNSigProton);
            
            if (abs(fTPCNSigProton) > 3.) continue; // PID for proton
            if (abs(fTPCNSigPion1) > 3.) continue; // PID for 1st pion
            if (abs(fTPCNSigPion2) > 3.) continue; // PID for 2nd pion
                fHistos->FillTH2("hTPCPIDLambdaProton_cut",pTrackXi->GetTPCmomentum(),fTPCNSigProton);
                fHistos->FillTH2("hTPCPIDLambdaPion_cut",nTrackXi->GetTPCmomentum(),fTPCNSigProton);
                fHistos->FillTH2("hTPCPIDBachelorPion_cut",bTrackXi->GetTPCmomentum(),fTPCNSigProton);
            
            // DCA cut
            // DCA between Dautgher particles
            Double_t fDCADist_Lambda = fabs(Xicandidate->GetDcaV0Daughters());
            Double_t fDCADist_Xi = fabs(Xicandidate->GetDcaXiDaughters());
                fHistos -> FillTH1("hDCADist_Lambda_BTW_Daughters",fDCADist_Lambda);
                fHistos -> FillTH1("hDCADist_Xi_BTW_Daughters",fDCADist_Xi);
            
            if( fDCADist_Lambda > 1.6) continue;// DCA proton-pion
            if( fDCADist_Xi > 1.6) continue;// DCA Lambda-pion
                fHistos -> FillTH1("hDCADist_Lambda_BTW_Daughters_cut",fDCADist_Lambda);
                fHistos -> FillTH1("hDCADist_Xi_BTW_Daughters_cut",fDCADist_Xi);
            // DCA to PV
            Double_t fDCADist_Lambda_PV       = fabs(Xicandidate->GetD(PVx, PVy, PVz));
            Double_t fDCADist_Xi_PV           = fabs(Xicandidate->GetDcascade(PVx, PVy, PVz));
            Double_t fDCADist_LambdaProton_PV = fabs(pTrackXi->GetD(PVx, PVy, bField));
            Double_t fDCADist_LambdaPion_PV   = fabs(pTrackXi->GetD(PVx, PVy, bField));
            Double_t fDCADist_BachelorPion_PV = fabs(bTrackXi->GetD(PVx, PVy, bField));
                fHistos -> FillTH1("fDCADist_lambda_to_PV",fDCADist_Lambda_PV);
                fHistos -> FillTH1("hDCADist_Xi_to_PV",fDCADist_Xi_PV2);
                fHistos -> FillTH1("hDCADist_LambdaProton_to_PV",fDCADist_LambdaProton_PV);
                fHistos -> FillTH1("hDCADist_LambdaPion_to_PV",fDCADist_LambdaPion_PV);
                fHistos -> FillTH1("hDCADist_BachelorPion_to_PV",fDCADist_BachelorPion_PV);
            
            // CPA cut
            Double_t fLambdaCPA = Xicandidate->GetV0CosineOfPointingAngle(PVx, PVy, PVz);
            Double_t fXiCPA = Xicandidate->GetCascadeCosineOfPointingAngle(PVx, PVy, PVz);
                fHistos -> FillTH1("hCosPA_lambda",fLambdaCPA);
                fHistos -> FillTH1("hCosPA_Xi",fXiCPA);
            
            if(Xicandidate->GetV0CosineOfPointingAngle(PVx, PVy, PVz) < 0.97) continue;
            if(Xicandidate->GetCascadeCosineOfPointingAngle(PVx, PVy, PVz) < 0.97) continue;
                fHistos -> FillTH1("hCosPA_lambda_cut",fLambdaCPA);
                fHistos -> FillTH1("hCosPA_Xi_cut",fXiCPA);
            
            // Mass window cut
            Double_t fMass_Xi = Xicandidate->M();
                fHistos -> FillTH1("hMass_Xi",fMass_Xi);
            if (fabs(fMass_Xi - Ximass) > 0.007) continue;
                fHistos -> FillTH1("hMass_Xi_cut",fMass_Xi);
            
            // Eta cut
            if(abs(Xicandidate->Eta())>0.8) continue;
            
            fHistos->FillTH2("hPhiEta_Xi",Xicandidate->Phi(),Xicandidate->Eta());
        }
        else {
            // !! NEED TO MODIFY !!
            /*
            track = (AliAODTrack*) fEvt ->GetTrack(it);
            if (!track) continue;
            if( ! ((AliAODTrack*) track)->TestFilterBit(fFilterBit)) continue;
            if (track->Pt()<fptcut) continue;
            if (abs(track->Eta())>fetacut) continue;
            fHistos->FillTH2("hPhiEta",track->Phi(),track->Eta());
            */
        }
        fNCascade++;
        goodcascadeindices.push_back(it);
    }
    
    return goodcascadeindices.size();
}

void AliAnalysisTaskXi1530::FillTracks(){
    
    AliVTrack *track1;
    // charged track, pion
    AliESDcascade *Xicandidate;
    // Cascade
    
    Bool_t mcXiFilled = kFALSE; // So that mctracks are never used uninitialized
    
    TLorentzVector temp1,temp2;
    TLorentzVector vecsum;
    
    const UInt_t ncascade = goodcascadeindices.size();
    const UInt_t ntracks = goodtrackindices.size();
    
    tracklist trackpool;
    if (fsetmixing){
        eventpool &ep = fEMpool[centbin][zbin];
        if (ep.size()<5 ) return;
        for (auto pool: ep){
            for (auto track: pool) trackpool.push_back((AliVTrack*)track);
        }
    }
    
    for (Int_t i = 0; i < ncascade; i++) {
        Xicandidate = ((AliESDEvent*)fEvt)->GetCascade(goodcascadeindices[i]);
        if(!Xicandidate) continue;
        
        AliESDtrack *pTrackXi   = ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs( Xicandidate->GetPindex()));
        AliESDtrack *nTrackXi   = ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs( Xicandidate->GetNindex()));
        AliESDtrack *bTrackXi   = ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs( Xicandidate->GetBindex()));
        
        temp1.SetXYZM(Xicandidate->Px(),Xicandidate->Py(), Xicandidate->Pz(), Xicandidate->M());
        
        // QA Plot for Good Xi
        FillTHnSparse("hInvMass_dXi",{fCent,Xicandidate->Pt(),Xicandidate->M()});
        
        for (UInt_t i = 0; i < ntracks; i++) {
            track1 = (AliVTrack*) fEvt->GetTrack(goodtrackindices[i]);
            if (!track1) continue;
            temp2.SetXYZM(track1->Px(),track1->Py(), track1->Pz(),pionmass);
            vecsum = temp1+temp2; // two pion vector sum
            
            // Y cut
            if (fabs(vecsum.Rapidity())>0.5) continue;
            
            auto sign = kAllType;
            if (Xicandidate->Charge()*track1->Charge()==-1) sign = kData; //Unlike sign -> Data
            else sign = kLS; //like sign bg
            
            mcXiFilled = kFALSE;
            if (IsMC) {
                if (fEvt->IsA()==AliESDEvent::Class()){
                    TParticle* MCXiD2esd = (TParticle*)fMCStack->Particle(abs(bTrackXi->GetLabel()));
                    TParticle* MCLamD1esd;
                    TParticle* MCLamD2esd;
                    TParticle* MCLamesd;
                    TParticle* MCXiesd;
                    TParticle* MCXiStaresd;
                    TParticle* MCXiStarD2esd;
                    
                    if (abs(MCXiD2esd->GetPdgCode()) == kPionCode) {
                        MCLamD1esd = (TParticle*)fMCStack->Particle(abs(pTrackXi->GetLabel()));
                        MCLamD2esd = (TParticle*)fMCStack->Particle(abs(nTrackXi->GetLabel()));
                        if (MCLamD1esd->GetMother(0) == MCLamD2esd->GetMother(0)) {
                            if ((abs(MCLamD1esd->GetPdgCode()) == kProtonCode && abs(MCLamD2esd->GetPdgCode()) == kPionCode) || (abs(MCLamD1esd->GetPdgCode()) == kPionCode && abs(MCLamD2esd->GetPdgCode()) == kProtonCode)) {
                                MCLamesd = (TParticle*)fMCStack->Particle(abs(MCLamD1esd->GetMother(0)));
                                if (abs(MCLamesd->GetPdgCode()) == kLambdaCode) {
                                    if (MCLamesd->GetMother(0) == MCXiD2esd->GetMother(0)) {
                                        MCXiesd = (TParticle*)fMCStack->Particle(abs(MCLamesd->GetMother(0)));
                                        if (abs(MCXiesd->GetPdgCode()) == kXiCode) {
                                            MCXiStarD2esd = (TParticle*)fMCStack->Particle(track1->GetLabel());
                                            if (MCXiesd->GetMother(0) == MCXiStarD2esd->GetMother(0)) {
                                                MCXiStaresd = (TParticle*)fMCStack->Particle(abs(MCXiesd->GetMother(0)));
                                                if (abs(MCXiStaresd->GetPdgCode()) == kXiStarCode) {
                                                    temp1.SetXYZM(MCXiesd->Px(),MCXiesd->Py(), MCXiesd->Pz(),Ximass);
                                                    temp2.SetXYZM(MCXiStarD2esd->Px(),MCXiStarD2esd->Py(), MCXiStarD2esd->Pz(),pionmass);
                                                    auto vecsumtrue = temp1 + temp2;
                                                    FillTHnSparse("hInvMassMCXi",{fCent,vecsumtrue.Pt(),vecsumtrue.M()});
                                                }//Xi1530 check
                                            }// Xi+pion mother check
                                        }// Xi Check
                                    }// Lambda+pion(D2esd) mother check
                                }//Lambda check
                            }//Lamda daugthers check
                        }//Same mother(lambda)
                    }//D2esd->pion
                }// ESD
                else{
                    // !! NEED TO UPDATE FOR AOD CASE !!
                    /*
                    AliAODMCParticle *par1 =
                    dynamic_cast<AliAODMCParticle*>(fMCArray->At(track1->GetLabel()));
                    AliAODMCParticle *par2 =
                    dynamic_cast<AliAODMCParticle*>(fMCArray->At(track2->GetLabel()));
                    Bool_t IsTrueOk = false;
                    if (!par1 || !par2) continue;
                    Int_t pdg1 = par1 -> PdgCode();
                    Int_t pdg2 = par2 -> PdgCode();
                    Int_t pdgtype = AliPID::ParticleCode(fParticleType);
                    Double_t pdgmass = AliPID::ParticleMass(fParticleType);
                    if (abs(pdg1) != pdgtype || abs(pdg2) != pdgtype) continue;
                    temp1.SetXYZM(par1->Px(),par1->Py(), par1->Pz(),pdgmass);
                    temp2.SetXYZM(par2->Px(),par2->Py(), par2->Pz(),pdgmass);
                    auto vecsumtrue = temp1 + temp2;
                    FillTHnSparse("hPtInvMResponse",{(double)sign,fCent
                        ,vecsum.Pt(),vecsumtrue.Pt(), vecsum.M(),vecsumtrue.M()});
                     */
                }// AOD
            }// MC
            FillTHnSparse("hInvMass",{(double)sign,fCent,vecsum.Pt(),vecsum.M()});
        }
    }
    if (fsetmixing){
        for (Int_t i = 0; i < ncascade; i++) {
            Xicandidate = ((AliESDEvent*)fEvt)->GetCascade(goodcascadeindices[i]);
            if(!Xicandidate) continue;
            temp1.SetXYZM(Xicandidate->Px(),Xicandidate->Py(), Xicandidate->Pz(), Xicandidate->M());
            
            for (UInt_t jt = 0; jt < ntracks; jt++) {
                track1 = trackpool.at(jt);
                temp2.SetXYZM(track1->Px(),track1->Py(), track1->Pz(),pionmass);
                vecsum = temp1+temp2; // two pion vector sum
                if (track1->Charge()*Xicandidate->Charge() == -1) continue;
                if (fabs(vecsum.Eta())>0.5) continue; //rapidity cut
                FillTHnSparse("hInvMass",{kMixing,fCent,vecsum.M(),vecsum.Pt()});
            }
        }
    }
}

void AliAnalysisTaskXi1530::Terminate(Option_t *)
{
}
Int_t AliAnalysisTaskXi1530::GetPID(AliPIDResponse *pid, const AliVTrack *trk){
    if (!pid) return -1; // no pid available
    
    /*Double_t sigmas[] ={-999,-999,-999,-999};
     
     Int_t ipid = kUnknown;
     Double_t lsigma = 3;
     sigmas[kPion] = pid -> NumberOfSigmasTPC(trk,AliPID::kPion);
     sigmas[kKaon] = pid -> NumberOfSigmasTPC(trk,AliPID::kKaon);
     sigmas[kProton] = pid -> NumberOfSigmasTPC(trk,AliPID::kProton);
     sigmas[kElectron] = pid -> NumberOfSigmasTPC(trk,AliPID::kElectron);
     for (int i=0; i<kUnknown; i++){
     if (fabs(sigmas[i]) < lsigma) {
     lsigma = fabs(sigmas[i]);
     ipid = i;
     }
     }
     
     // derive information, whether tof pid is available
     if (0){
     const Bool_t ka = !(trk->GetStatus() & AliESDtrack::kTOFmismatch);
     const Bool_t kb =  (trk->GetStatus() & AliESDtrack::kTOFpid);
     const Bool_t ktof = ka && kb;
     }
     
     if (lsigma>3 ) return kUnknown;
     else  return ipid;
     
     */
    Double_t prob[AliPID::kSPECIES];
    fPIDCombined->ComputeProbabilities(trk,pid,prob);
    Int_t ipid = AliPID::kUnknown;
    Double_t iprob = 0;
    for (int i=0; i<AliPID::kSPECIES; i++){
        if (prob[i]>0.6 && prob[i]>iprob) {
            iprob = prob[i];
            ipid = i;
        }
    }
    if (ipid == AliPID::kUnknown) ipid = AliPID::kPion;
    return ipid;
    
}

THnSparse* AliAnalysisTaskXi1530::CreateTHnSparse(TString name, TString title, Int_t ndim, std::vector<TAxis> bins, Option_t * opt){
    const TAxis * axises[bins.size()];
    for( UInt_t i=0;i<bins.size();i++ ) axises[i]= &bins[i];
    THnSparse * h= fHistos->CreateTHnSparse(name, title, ndim, axises,opt );
    return h;
}

Long64_t AliAnalysisTaskXi1530::FillTHnSparse( TString name, std::vector<Double_t> x, Double_t w ){
    auto hsparse = dynamic_cast<THnSparse*>( fHistos->FindObject(name) );
    if(! hsparse ){
        std::cout<<"ERROR : no "<<name<<std::endl;
        exit(1);
    }
    return FillTHnSparse( hsparse, x, w );
}

Long64_t AliAnalysisTaskXi1530::FillTHnSparse( THnSparse *h, std::vector<Double_t> x, Double_t w ){
    if( int(x.size()) != h->GetNdimensions() ){
        std::cout<<"ERROR : wrong sized of array while Fill "<<h->GetName()<<std::endl;
        exit(1);
    }
    return h->Fill( &x.front(), w );
}

TAxis AliAnalysisTaskXi1530::AxisFix
( TString name, int nbin, Double_t xmin, Double_t xmax ){
    TAxis axis(nbin, xmin, xmax);axis.SetName(name);
    return axis;
}

TAxis AliAnalysisTaskXi1530::AxisStr( TString name, std::vector<TString> bin ){
    TAxis ax = AxisFix( name, bin.size(), 0.5, bin.size()+0.5);
    UInt_t i=1;
    for( auto blabel : bin )
        ax.SetBinLabel( i++, blabel );
    return ax;
}

TAxis AliAnalysisTaskXi1530::AxisVar( TString name, std::vector<Double_t> bin ){
    TAxis axis( bin.size()-1, &bin.front() ) ;axis.SetName(name);
    return axis;
}

TAxis AliAnalysisTaskXi1530::AxisLog
( TString name, int nbin, Double_t xmin, Double_t xmax, Double_t xmin0){
    int binoffset = ( xmin0<0 || (xmin-xmin0)<1e-9) ? 0 : 1;
    std::vector<Double_t> bin(nbin+1+binoffset,0);
    double logBW3 = (log(xmax)-log(xmin))/nbin;
    for(int ij=0;ij<=nbin;ij++) bin[ij+binoffset]=xmin*exp(ij*logBW3);
    TAxis axis( nbin, &bin.front() ) ;
    axis.SetName(name);
    return axis;
}

