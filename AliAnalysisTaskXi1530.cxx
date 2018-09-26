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
#include "TList.h"
#include "TSystem.h"
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
//#include "AliAODMCHeader.h"
//#include "AliAODMCParticle.h"
#include "AliMultiplicity.h"

// Some constants
const Double_t        pi = TMath::Pi();
const Double_t  pionmass = AliPID::ParticleMass(AliPID::kPion);
const Double_t  kaonmass = AliPID::ParticleMass(AliPID::kKaon);
enum {  kPN=1, kPP, kNN, kMixing, kAllType}; //P=Positive charge, N=Negative


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
    delete fTrigger;
    delete fTrackCuts;
    delete fPIDResponse;
    delete fRunTable;
}
/*
void AliAnalysisTaskXi1530::XiStarInit()
{
    //
    //Inits cuts and analysis settings
    //
    fEventCounter = 0; // event counter initialization
    if (fDevelopeMode)std::cout << "AliAnalysisTaskXi1530 XiStarInit() call" << std::endl;
    if (fDevelopeMode)std::cout << "MC Mode?: " << fMCcase << std::endl;
    
    ///////////////////////////////////////////////
    // Track Cuts for ESD analysis
    fTrackCut = new AliESDtrackCuts();
    fTrackCut->SetPtRange(.15, 1000);
    fTrackCut->SetAcceptKinkDaughters(kFALSE);
    //fTrackCut->SetMinNClustersTPC(50);
    fTrackCut->SetRequireTPCRefit(kTRUE);
    fTrackCut->SetMaxChi2PerClusterTPC(4); //From Enrico
    
    ////////////////////////////////////////////////
    
    fZvertexBins = 20;
    fMultBins = 11;// This must also be set in AliAnalysisTaskXi1530.h
    if (fMCcase) fEventsToMix = 0;
    else fEventsToMix = 40; // original 40 jisong
    
    // multiplicity edges for event mixing bins
    fMultLimits[0] = 0, fMultLimits[1] = 5, fMultLimits[2] = 10, fMultLimits[3] = 15, fMultLimits[4] = 20, fMultLimits[5] = 25;
    fMultLimits[6] = 30, fMultLimits[7] = 35, fMultLimits[8] = 40, fMultLimits[9] = 45, fMultLimits[10] = 50, fMultLimits[11] = 150;
    
    
    fEC = new AliAnalysisTaskXi1530EventCollection **[fZvertexBins];
    for (unsigned short i = 0; i < fZvertexBins; i++) {
        
        fEC[i] = new AliAnalysisTaskXi1530EventCollection *[fMultBins];
        
        for (unsigned short j = 0; j < fMultBins; j++) {
            
            fEC[i][j] = new AliAnalysisTaskXi1530EventCollection(fEventsToMix + 1);
        }
    }
    
    fTempStruct = new AliAnalysisTaskXi1530TrackStruct[kNbinsM * 8];
    fESDTrack4 = new AliESDtrack();
    fXiTrack = new AliESDtrack();
    
    
    fMaxDecayLength = 100.;
    fMassWindow = 0.007;
    
    /////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////
    // DecayParameters Key (number represents array index)
    // NclustersTPC: 0=proton, 1=pion first, 2=pion second, 3=pion third
    // DCAVtx: 4=proton, 5=pion first, 6=pion second, 7=lambda, 8=pion third
    // 9 = DCA proton-pion
    // 10 = DCA Lambda-pion
    // 11 = Rxy Lambda
    // 12 = Rxy Xi
    // 13 = Cos PA Lambda
    // 14 = Cos PA Xi
    
    // Set Standard Reconstruction cut values
    fCutValues[0][0] = 50;
    fCutValues[0][1] = 50;
    fCutValues[0][2] = 50;
    fCutValues[0][3] = 50; // for 2010 cut (origin: 70)
    fCutValues[0][4] = 0.04;
    fCutValues[0][5] = 0.04;
    fCutValues[0][6] = 0.05;
    fCutValues[0][7] = 0.07;
    fCutValues[0][8] = 2.0;
    fCutValues[0][9] = 1.6;
    fCutValues[0][10] = 1.6;
    fCutValues[0][11] = 0.97;
    fCutValues[0][12] = 0.97;
    
    
    for (int cv = 1; cv < kNCutVariations; cv++) {
        for (int ct = 0; ct < kNCuts; ct++) {
            fCutValues[cv][ct] = fCutValues[0][ct];
        }
    }
    
    //systematic variation// Loose
    fCutValues[1][0] = 45;
    fCutValues[1][1] = 45;
    fCutValues[1][2] = 45;
    fCutValues[1][3] = 45;// 60 -> 45
    fCutValues[2][4] = 0.03;
    fCutValues[3][5] = 0.03;
    fCutValues[4][6] = 0.04;
    fCutValues[5][7] = 0.06;
    fCutValues[6][8] = 2.1;
    fCutValues[7][9] = 1.7;
    fCutValues[8][10] = 1.7;
    fCutValues[9][11] = 0.95;
    fCutValues[10][12] = 0.965;
    
    //systematic variation// tight
    fCutValues[11][0] = 55;
    fCutValues[11][1] = 55;
    fCutValues[11][2] = 55;
    fCutValues[11][3] = 55;// 70 -> 55
    fCutValues[12][4] = 0.104;
    fCutValues[13][5] = 0.104;
    fCutValues[14][6] = 0.08;
    fCutValues[15][7] = 0.1;
    fCutValues[16][8] = 1.0;
    fCutValues[17][9] = 0.94;
    fCutValues[18][10] = 1.41;
    fCutValues[19][11] = 0.99;
    fCutValues[20][12] = 0.985;
    
    
    // PDG mass values
    fTrueMassPr = .93827, fTrueMassPi = .13957, fTrueMassK = .493677, fTrueMassLam = 1.11568, fTrueMassXi = 1.32171;
    
    // The following CovMatrix is set so that PropogateToDCA() ignores track errors. Only used to propagate Xi to third pion for XiStar reconstruction
    for (Int_t i = 0; i < 21; i++) fCovMatrix[i] = 0;
    fCovMatrix[0] = 1, fCovMatrix[2] = 1, fCovMatrix[5] = 1, fCovMatrix[9] = 1, fCovMatrix[14] = 1, fCovMatrix[20] = 1;
    
    
}
 */
//________________________________________________________________________
void AliAnalysisTaskXi1530::UserCreateOutputObjects()
{
    std::cout << "User Create Object" << std::endl;
    // Histograms container
    fOutput = new TList();
    fOutput->SetOwner(kTRUE);
    
    // Offline triggers -----------------------------------------------------
    //fTrigger = new AliTriggerAnalysis; // offline trigger
    //fTrigger -> SetFMDThreshold(0.3,0.5); // FMD threshold
    //-----------------------------------------------------------------------
    
    // TrackCuts for Xi1530--------------------------------------------------
    fTrackCuts = new AliESDtrackCuts();
    fTrackCuts -> GetStandardITSTPCTrackCuts2011(kTRUE,kTRUE);
    fTrackCuts -> SetEtaRange(-0.8,0.8);
    fTrackCuts -> SetPtRange(0.15, 1e20);
    // ----------------------------------------------------------------------
    
    fHistos = new THistManager("Xi1530hists");
    
    auto binType = AxisStr("Type",{"PN","PP","NN","NP","Mixing"});
    if (IsAA) binCent = AxisFix("Cent",10,0,300);
    else binCent = AxisFix("Cent",1,0,300);
    auto binPt   = AxisFix("Pt",200,0,20);
    auto binMass = AxisFix("Mass",1000,10,20);
    
    CreateTHnSparse("hInvMass","InvMass",4,{binType,binCent,binPt,binMass},"s");
    //CreateTHnSparse("hInvMass2track","InvMass",3,{binCent,binPt,binMass},"s");
    //CreateTHnSparse("hInvMass4track","InvMass",3,{binCent,binPt,binMass},"s");
    //CreateTHnSparse("hInvMass6track","InvMass",3,{binCent,binPt,binMass},"s");
    CreateTHnSparse("hMult","Multiplicity",1,{binCent},"s");
    CreateTHnSparse("hPtInvMResponse","InvMass Res",6
                    ,{binType,binCent,binPt,binPt,binMass,binMass},"s");
    
    vector<TString> ent = {"All","PS","PSpileup","Goodz","Goodzcut"};
    auto hNofEvt = fHistos->CreateTH1("hEventNumbers","",ent.size(), 0, ent.size());
    for(auto i=0u;i<ent.size();i++) hNofEvt->GetXaxis()->SetBinLabel(i+1,ent.at(i).Data());
    
    fHistos -> CreateTH2("hPhiEta","",180,0,2*pi,40,-2,2);
    
    binZ = AxisVar("Z",{-10,-5,-3,-1,1,3,5,10});
    fHistos->CreateTH1("hMul","",200,0,200,"s");
    fHistos->CreateTH1("hZvtx","",600,-30,30,"s");
    fEMpool.resize(binCent.GetNbins(),vector<eventpool> (binZ.GetNbins()));
    
    if (IsMC)
    {
        // To get Trigger efficiency in each trk/V0M Multiplicity region
        
        vector<TString> mcent = {"All","IsINELg0","tracklet in |Eta|<1","CINT7 triggered","AliMultiSelection"};
        vector<Int_t> ntrklet = {0,5,10,15,20,25,30,35,40,50};  // # of Tracklet bin
        vector<Int_t> nmult = {0,1,5,10,15,20,30,40,50,70,100}; // V0M Multiplicity bin
        
        auto hNofEvtMC = fHistos->CreateTH1("htotalEvent_MC","",mcent.size(), 0, mcent.size());
        auto htrkINELg0 = fHistos->CreateTH1("htriggered_INELg0_tracklet","",ntrklet.size(), 0, ntrklet.size());
        auto htrkCINT7 = fHistos->CreateTH1("htriggered_CINT7_tracklet","",ntrklet.size(), 0, ntrklet.size());
        auto hmultCINT7 = fHistos->CreateTH1("htriggered_CINT7_VOM","",nmult.size(), 0, nmult.size());
        
        for(auto i=0u;i<mcent.size();i++) hNofEvtMC->GetXaxis()->SetBinLabel(i+1,mcent.at(i).Data());
        htrkINELg0->GetXaxis()->SetBinLabel(1,"0 to Inf (MB)");
        for(auto i=1u;i<ntrklet.size();i++) htrkINELg0->GetXaxis()->SetBinLabel(i+1,Form("%d to %d",ntrklet.at(i),ntrklet.at(i)+5));
        htrkCINT7->GetXaxis()->SetBinLabel(1,"0 to Inf (MB)");
        for(auto i=1u;i<ntrklet.size();i++) htrkCINT7->GetXaxis()->SetBinLabel(i+1,Form("%d to %d",ntrklet.at(i),ntrklet.at(i)+5));
        hmultCINT7->GetXaxis()->SetBinLabel(1,"0 - 100 % (MB)");
        for(auto i=1u;i<nmult.size();i++) hmultCINT7->GetXaxis()->SetBinLabel(i+1,Form("%d - %d%%",nmult.at(i),nmult.at(i+1)));
        
        CreateTHnSparse("hMult_MC","Multiplicity",1,{binCent},"s");
        CreateTHnSparse("hMult_MC_selected","Multiplicity",1,{binCent},"s");
    }

    PostData(1, fHistos->GetListOfHistograms());
}

//________________________________________________________________________
void AliAnalysisTaskXi1530::UserExec(Option_t *)
{
    std::cout << "UserExec" << std::endl;
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
    //       150-250: Not selected
    //       300: No MultSection
    //
    fCent = 150; // Multiplicity
    AliMultSelection *MultSelection = (AliMultSelection*) fEvt->FindListObject("MultSelection");
    if(MultSelection)
    {
        if (!(MultSelection->IsEventSelected()))
        {
            AliInfo("This event is not selected: AliMultSelection");
            fCent += MultSelection->GetMultiplicityPercentile("V0M");
        }
        fCent = MultSelection->GetMultiplicityPercentile("V0M");
    }
    else
    {
        //If this happens, re-check if AliMultSelectionTask ran before your task!
        AliInfo("Didn't find MultSelection!");
        fCent = 300;
    }
    // ----------------------------------------------------------------------
    
    // Preparation for MC ---------------------------------------------------
    if (IsMC)
    {
        AliMCEvent  *mcEvent        = 0x0;
        AliVEventHandler* eventHandler = AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler();
        if(eventHandler){
            AliMCEventHandler* mcEventHandler = dynamic_cast<AliMCEventHandler*>(eventHandler);
            if(mcEventHandler) mcEvent = static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->MCEvent();
        }
        if(!mcEvent) return;
        
        fMCArray = (TClonesArray*) fEvt->FindListObject("mcparticles");
        fMCStack = (AliStack*) mcEvent->Stack();
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
    
    if(fRunTable->IsAA() || fRunTable->IsPA()){
        IsMinimumBias = (inputHandler -> IsEventSelected()) & (AliVEvent::kMB);
    }
    if (IsMinimumBias) fHistos -> FillTH1("hEventNumbers","PS",1);
    
    /*
    //
    AliMCEvent  *mcEvent        = 0x0;
    AliStack    *mcstack        = 0x0;
    if (fDevelopeMode)std::cout << "TEST" << std::endl;
    //for mc study
    Bool_t IsINELg0 = kFALSE;
    Bool_t isSelectedkINT7 = kFALSE;
    
    if (fMCcase) {
        AliVEventHandler* eventHandler = AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler();
        if(eventHandler){
            AliMCEventHandler* mcEventHandler = dynamic_cast<AliMCEventHandler*>(eventHandler);
            if(mcEventHandler) mcEvent = static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->MCEvent();
        }
        if (!mcEvent) {
            if (fDevelopeMode)std::cout << "Return: No mcEvent" << std::endl;
            return;
        }
        mcstack = mcEvent->Stack();
        if (!mcstack) {
            if (fDevelopeMode)std::cout << "Return: No mcStack" << std::endl;
            return;
        }
        
        htotalEvent->Fill(0); // Total N of event
        
        // for INELg0 check.
        for (Int_t it = 0; it < mcstack->GetNprimary(); it++) {
            AliMCParticle *mcInputTrack = (AliMCParticle*)mcEvent->GetTrack(it);
            if (mcInputTrack->IsPhysicalPrimary() && TMath::Abs(mcInputTrack->Charge()) && TMath::Abs(mcInputTrack->Eta()) < 1 ) IsINELg0 = true;
        }
        if (IsINELg0) {
            if (fDevelopeMode)std::cout << "It's INELg0 event!" << std::endl;
            htotalEvent->Fill(1); // Total N of INELg0 event
            
            const AliMultiplicity* mult = fESD->GetMultiplicity();
            Int_t fSpdT_origin = mult->GetNumberOfTracklets();
            if (fDevelopeMode)std::cout << "# of total tracklet : " << fSpdT_origin << std::endl;
            
            Int_t fNSpdT = 0;
            for (Int_t i = 0; i < fSpdT_origin; ++i) if (TMath::Abs(mult->GetEta(i)) < 0.8) fNSpdT++;
            if (fDevelopeMode)std::cout << "# of passed tracklet : " << fNSpdT << std::endl;
            
            // |Eta| < 1
            if (fNSpdT > 0) {
                htotalEvent->Fill(2); // Total N of triggered event.
                
                if (fNSpdT > 0) htriggered_INELg0_tracklet->Fill(0); // INEL>0
                
                if (fNSpdT > 51) htriggered_INELg0_tracklet->Fill(10); // INEL>0
                else if (fNSpdT > 41) htriggered_INELg0_tracklet->Fill(9); // INEL>0
                else if (fNSpdT > 36) htriggered_INELg0_tracklet->Fill(8); // INEL>0
                else if (fNSpdT > 31) htriggered_INELg0_tracklet->Fill(7); // INEL>0
                else if (fNSpdT > 26) htriggered_INELg0_tracklet->Fill(6); // INEL>0
                else if (fNSpdT > 21) htriggered_INELg0_tracklet->Fill(5); // INEL>0
                else if (fNSpdT > 16) htriggered_INELg0_tracklet->Fill(4); // INEL>0
                else if (fNSpdT > 11) htriggered_INELg0_tracklet->Fill(3); // INEL>0
                else if (fNSpdT > 6) htriggered_INELg0_tracklet->Fill(2); // INEL>0
                else if (fNSpdT > 0) htriggered_INELg0_tracklet->Fill(1); // INEL>0
                
                // CINT7 Triggered event.
                isSelectedkINT7 = inputHandler->IsEventSelected() & AliVEvent::kINT7;
                if (isSelectedkINT7) {
                    
                    htotalEvent->Fill(3); // Total N of triggered event.
                    
                    if (fNSpdT > 0) htriggered_CINT7_tracklet->Fill(0); // INEL>0 in CINT7
                    
                    if (fNSpdT > 51) htriggered_CINT7_tracklet->Fill(10); // INEL>0 in CINT7
                    else if (fNSpdT > 41) htriggered_CINT7_tracklet->Fill(9); // INEL>0 in CINT7
                    else if (fNSpdT > 36) htriggered_CINT7_tracklet->Fill(8); // INEL>0 in CINT7
                    else if (fNSpdT > 31) htriggered_CINT7_tracklet->Fill(7); // INEL>0 in CINT7
                    else if (fNSpdT > 26) htriggered_CINT7_tracklet->Fill(6); // INEL>0 in CINT7
                    else if (fNSpdT > 21) htriggered_CINT7_tracklet->Fill(5); // INEL>0 in CINT7
                    else if (fNSpdT > 16) htriggered_CINT7_tracklet->Fill(4); // INEL>0 in CINT7
                    else if (fNSpdT > 11) htriggered_CINT7_tracklet->Fill(3); // INEL>0 in CINT7
                    else if (fNSpdT > 6) htriggered_CINT7_tracklet->Fill(2); // INEL>0 in CINT7
                    else if (fNSpdT > 0) htriggered_CINT7_tracklet->Fill(1); // INEL>0 in CINT7
                    
                    
                    // AliMultSelection
                    AliMultSelection *MultSelection = (AliMultSelection*) fESD->FindListObject("MultSelection");
                    Float_t lPerc = 300; //nonsense
                    
                    //Quality check                                                                          // it should be same with 1.
                    lPerc = MultSelection->GetMultiplicityPercentile("V0M");
                    fMultDist_MCpp->Fill(lPerc);
                    
                    if (lPerc < 101) htriggered_CINT7_VOM->Fill(0); // INEL>0 in CINT7
                    
                    if (lPerc > 70) htriggered_CINT7_VOM->Fill(1); // INEL>0 in CINT7
                    else if (lPerc > 50) htriggered_CINT7_VOM->Fill(2); // INEL>0 in CINT7
                    else if (lPerc > 40) htriggered_CINT7_VOM->Fill(3); // INEL>0 in CINT7
                    else if (lPerc > 30) htriggered_CINT7_VOM->Fill(4); // INEL>0 in CINT7
                    else if (lPerc > 20) htriggered_CINT7_VOM->Fill(5); // INEL>0 in CINT7
                    else if (lPerc > 15) htriggered_CINT7_VOM->Fill(6); // INEL>0 in CINT7
                    else if (lPerc > 10) htriggered_CINT7_VOM->Fill(7); // INEL>0 in CINT7
                    else if (lPerc > 5) htriggered_CINT7_VOM->Fill(8); // INEL>0 in CINT7
                    else if (lPerc > 1) htriggered_CINT7_VOM->Fill(9); // INEL>0 in CINT7
                    else if (lPerc > 0) htriggered_CINT7_VOM->Fill(10); // INEL>0 in CINT7
                    
                    if (MultSelection->IsEventSelected()) {
                        htotalEvent->Fill(4); // Total N of Multi selected event
                        fMultDist_MCpp_selected->Fill(lPerc);
                        
                        if (lPerc < 101) htriggered_AliMult_VOM->Fill(0); // INEL>0 in CINT7
                        
                        if (lPerc > 70) htriggered_AliMult_VOM->Fill(1); // INEL>0 in CINT7
                        else if (lPerc > 50) htriggered_AliMult_VOM->Fill(2); // INEL>0 in CINT7
                        else if (lPerc > 40) htriggered_AliMult_VOM->Fill(3); // INEL>0 in CINT7
                        else if (lPerc > 30) htriggered_AliMult_VOM->Fill(4); // INEL>0 in CINT7
                        else if (lPerc > 20) htriggered_AliMult_VOM->Fill(5); // INEL>0 in CINT7
                        else if (lPerc > 15) htriggered_AliMult_VOM->Fill(6); // INEL>0 in CINT7
                        else if (lPerc > 10) htriggered_AliMult_VOM->Fill(7); // INEL>0 in CINT7
                        else if (lPerc > 5) htriggered_AliMult_VOM->Fill(8); // INEL>0 in CINT7
                        else if (lPerc > 1) htriggered_AliMult_VOM->Fill(9); // INEL>0 in CINT7
                        else if (lPerc > 0) htriggered_AliMult_VOM->Fill(10); // INEL>0 in CINT7
                        
                    }// IsEventSelected in AliMultSelection
                }//eta < 1
            }//CINT7
        }//INEL>0
    }
    
    // check  : events are selected by physics selection class
    UInt_t fSelectMask = fInputHandler->IsEventSelected();
    
    isSelectedkINT7 = fSelectMask & AliVEvent::kINT7;
    if(isSelectedkINT7) ((TH1F*)hEventSelecInfo)->Fill(2);
    
    Bool_t isSelectedkHighMultV0 = fSelectMask & AliVEvent::kHighMultV0;
    if(isSelectedkHighMultV0) ((TH1F*)hEventSelecInfo)->Fill(3);
    
    
    if (fHMTrigger) {
        if (!isSelectedkHighMultV0) {
            if (fDevelopeMode)std::cout << "Event Rejected: No kHighMultV0 trigger" << std::endl;
            return;
        }
    }
    else {
        if (!isSelectedkINT7) {
            if (fDevelopeMode)std::cout << "Event Rejected: No kINT7 trigger" << std::endl;
            return;
        }
    }
    
    if (fDevelopeMode)std::cout << "TEST3" << std::endl;
    ///////////////////////////////////////////////////////////
    const AliESDVertex *PrimaryVertexESD;
    
    // ---- AliPIDResponse ---- //
    fPIDResponse = inputHandler->GetPIDResponse();
    double nSigTPCPID = 3.0;
    
    
    // TClonesArray *mcArray       = 0x0;
    TParticle   *MCLamD1esd     = 0x0;
    TParticle   *MCLamD2esd     = 0x0;
    TParticle   *MCLamesd       = 0x0;
    TParticle   *MCXiD2esd      = 0x0;
    TParticle   *MCXiesd        = 0x0;
    TParticle   *MCXiStarD2esd  = 0x0;
    TParticle   *MCXiStaresd    = 0x0;
    
    Double_t px1, py1, pz1, px2, py2, pz2;
    Double_t p1sq, p2sq, e1, e2, angle;
    Double_t dca3d;
    Float_t dca2[2];
    Double_t xiVtx[3];//, xiStarVtx[3];
    Double_t xiP[3], xiStarP[3];
    Double_t xiStarMom;
    Double_t xiMass, xiStarMass;
    Double_t xiPt, xiStarPt;
    Double_t xiY, xiStarY, MCxiY, MCxiStarY, MCxiYout, MCxiStarYout ;
    Double_t xiCharge;
    Double_t decayLengthXY;
    Double_t pDaughter1[3];
    Double_t pDaughter2[3];
    Double_t xDaughter1[3];
    Double_t xDaughter2[3];
    //
    Double_t bField = 0;
    UInt_t status = 0;
    Int_t positiveTracks = 0, negativeTracks = 0;
    Int_t myTracks = 0;
    Int_t myxiTracks = 0;
    Int_t myMCTracks = 0;
    //
    Double_t primaryVtx[3] = {0};
    Int_t mBin = 0;
    Int_t zBin = 0;
    Double_t zStep = 2 * 10 / Double_t(fZvertexBins), zStart = -10.;
    
    //
    Bool_t mcXiFilled = kFALSE; // So that mctracks are never used uninitialized
    
    fMultDist1->Fill(fESD->GetNumberOfTracks());
    PrimaryVertexESD = fESD->GetPrimaryVertex();
    if (!PrimaryVertexESD) return;
    hNumberOfEvent->Fill(0);
    fEventNumber = ((((ULong64_t)fESD->GetPeriodNumber()) << 36 ) | (((ULong64_t)fESD->GetOrbitNumber()) << 12 ) | ((ULong64_t)fESD->GetBunchCrossNumber()));
    
    primaryVtx[0] = PrimaryVertexESD->GetX();
    primaryVtx[1] = PrimaryVertexESD->GetY();
    primaryVtx[2] = PrimaryVertexESD->GetZ();
    fVertexDist1->Fill(primaryVtx[0], primaryVtx[1], primaryVtx[2]);
    
    if (fMCcase) {
        /////////////////////////////////////////////////
        // Lam mc input
        /////////////////////////////////////////////////
        for (Int_t it = 0; it < mcstack->GetNprimary(); it++) {
            TParticle *mcInputTrack = (TParticle*)mcstack->Particle(it);
            if (!mcInputTrack) {
                Error("UserExec", "Could not receive track %d", it);
                continue;
            }
            //fTreeVariableMCinputTotalSongPID = -1;
            //fTreeVariableEventNumber4 = GetMCEventNumber();
            if (mcInputTrack->GetPdgCode() != +kXiCode && mcInputTrack->GetPdgCode() != -kXiCode && mcInputTrack->GetPdgCode() != +kXiStarCode && mcInputTrack->GetPdgCode() != -kXiStarCode) continue;
            myMCTracks++;
            MCxiStarY = mcInputTrack->Y();
            fXiStarYDistMC->Fill(MCxiStarY);
            // Xi
            if (mcInputTrack->GetPdgCode() == +kXiCode) {
                hMCinputTotalXi1->Fill(mcInputTrack->Pt(), mcInputTrack->Y(), mcInputTrack->GetCalcMass());
            }
            if (mcInputTrack->GetPdgCode() == -kXiCode) {
                hMCinputTotalXibar1->Fill(mcInputTrack->Pt(), mcInputTrack->Y(), mcInputTrack->GetCalcMass());
            }
            // XiStar
            if (mcInputTrack->GetPdgCode() == +kXiStarCode) {
                hMCinputTotalXiStar1->Fill(mcInputTrack->Pt(), mcInputTrack->Y(), mcInputTrack->GetCalcMass());
            }
            if (mcInputTrack->GetPdgCode() == -kXiStarCode) {
                hMCinputTotalXiStarbar1->Fill(mcInputTrack->Pt(), mcInputTrack->Y(), mcInputTrack->GetCalcMass());
            }
        }
    }
    
    
    // Vertex systematic study default : 10 , loose : 11 , tight : 9 (cm)
    if (fabs(primaryVtx[2]) > 10.) return; // Z-Vertex Cut
    hNumberOfEvent->Fill(1);
    fMultDist2->Fill(fESD->GetNumberOfTracks());
    
    
    
    fMultDist3->Fill(fESD->GetNumberOfTracks());
    fVertexDist3->Fill(primaryVtx[0], primaryVtx[1], primaryVtx[2]);
    
    
    // Before the AliMulti
    hNumberOfEvent->Fill(2);
    
    // IncompleteDAQ Check
    if (fESD->IsIncompleteDAQ()) {
        if (fDevelopeMode)std::cout << "Reject: IsIncompleteDAQ" << std::endl;;
        return;
    }
    
    // Muliplicity Selection
    AliMultSelection *MultSelection = (AliMultSelection*) fESD->FindListObject("MultSelection");
    Float_t lPerc = 300; //nonsense
    
    if (MultSelection) {
        if (!(MultSelection->IsEventSelected())) {
            AliInfo("This event is not selected: AliMultSelection");
            lPerc += MultSelection->GetMultiplicityPercentile("V0M");
            
        }
        lPerc = MultSelection->GetMultiplicityPercentile("V0M");
    }
    else {
        //If this happens, re-check if AliMultSelectionTask ran before your task!
        AliInfo("Didn't find MultSelection!");
    }
    if (fDevelopeMode)std::cout << "Multiplicity: " << lPerc << std::endl;
    
    fMultDist_pp->Fill(lPerc);
    
    // After the AliMulti
    hNumberOfEvent->Fill(3);
    
    
    if (fDevelopeMode)std::cout << "There are " << fESD->GetNumberOfTracks() << " tracks in this event" << std::endl;;
    
    bField = fESD->GetMagneticField();
    
    // Track loop
    for (Int_t i = 0; i < fESD->GetNumberOfTracks(); i++) {
        AliESDtrack* esdtrack = fESD->GetTrack(i);
        if (!esdtrack) continue;
        status = esdtrack->GetStatus();
        
        if (!fTrackCut->AcceptTrack(esdtrack)) continue;
        
        Bool_t goodMomentum = esdtrack->GetPxPyPz(fTempStruct[myTracks].fP);
        if (!goodMomentum) continue;
        esdtrack->GetXYZ( fTempStruct[myTracks].fX);
        //=========checking PID =========//
        //// *** TPC *** ////
        Float_t fTPCPIDmom = esdtrack->GetTPCmomentum();
        Float_t sigTPC = esdtrack->GetTPCsignal();
        Float_t nsigpi = fabs(fPIDResponse->NumberOfSigmasTPC(esdtrack, AliPID::kPion));
        Float_t nsigk = fabs(fPIDResponse->NumberOfSigmasTPC(esdtrack, AliPID::kKaon));
        Float_t nsigpr = fabs(fPIDResponse->NumberOfSigmasTPC(esdtrack, AliPID::kProton));
        
        hTPCPID->Fill(fTPCPIDmom, sigTPC);
        if (nsigpi < 3.) hTPCPIDpi->Fill(fTPCPIDmom, sigTPC);
        if (nsigk < 3.) hTPCPIDk->Fill(fTPCPIDmom, sigTPC);
        if (nsigpr < 3.) hTPCPIDp->Fill(fTPCPIDmom, sigTPC);
        
        hNSig3rdPion->Fill(fTPCPIDmom, fPIDResponse->NumberOfSigmasTPC(esdtrack, AliPID::kPion));
     
        hQANSig3rdPion->Fill(fTPCPIDmom, fPIDResponse->NumberOfSigmasTPC(esdtrack, AliPID::kPion));
        
        esdtrack->GetCovarianceXYZPxPyPz( fTempStruct[myTracks].fCov);
        //esdtrack->GetImpactParameters(dca2, cov);
        dca2[0] = sqrt( pow(fTempStruct[myTracks].fX[0] - primaryVtx[0], 2) + pow(fTempStruct[myTracks].fX[1] - primaryVtx[1], 2));
        dca2[1] = sqrt( pow(fTempStruct[myTracks].fX[2] - primaryVtx[2], 2));
        dca3d = sqrt( pow(dca2[0], 2) + pow(dca2[1], 2));
        
        fDCADist->Fill(fESD->GetNumberOfTracks(), dca3d);
        // fPhiPtDist->Fill(esdtrack->Charge(), esdtrack->Phi(), esdtrack->Pt());
        // fPtEtaDist->Fill(esdtrack->Charge(), esdtrack->Pt(), esdtrack->Eta());
        
        fPtDist->Fill(esdtrack->Pt());
        fPhiDist->Fill(esdtrack->Phi());
        fEtaDist->Fill(esdtrack->Eta());
        
        fTempStruct[myTracks].fStatus = status;
        fTempStruct[myTracks].fID = esdtrack->GetID();
        fTempStruct[myTracks].fLabel = esdtrack->GetLabel();
        fTempStruct[myTracks].fPhi = atan2(fTempStruct[myTracks].fP[1], fTempStruct[myTracks].fP[0]);
        if (fTempStruct[myTracks].fPhi < 0) fTempStruct[myTracks].fPhi += 2 * PI;
        fTempStruct[myTracks].fPt = sqrt(pow(fTempStruct[myTracks].fP[0], 2) + pow(fTempStruct[myTracks].fP[1], 2));
        fTempStruct[myTracks].fMom = sqrt( pow(fTempStruct[myTracks].fPt, 2) + pow(fTempStruct[myTracks].fP[2], 2) );
        fTempStruct[myTracks].fEta = esdtrack->Eta();
        fTempStruct[myTracks].fCharge = esdtrack->Charge();
        fTempStruct[myTracks].fDCAXY = dca2[0];
        fTempStruct[myTracks].fDCAZ = dca2[1];
        fTempStruct[myTracks].fDCA = dca3d;
        fTempStruct[myTracks].fNSigmaPi = fabs(fPIDResponse->NumberOfSigmasTPC(esdtrack, AliPID::kPion));
        fTempStruct[myTracks].fNSigmaK = fabs(fPIDResponse->NumberOfSigmasTPC(esdtrack, AliPID::kKaon));
        fTempStruct[myTracks].fNSigmaPr = fabs(fPIDResponse->NumberOfSigmasTPC(esdtrack, AliPID::kProton));
        fTempStruct[myTracks].fNclusTPC = esdtrack->GetTPCNcls();
        
        if (esdtrack->Charge() > 0) positiveTracks++;
        else negativeTracks++;
        
        //if(fTempStruct[myTracks].fNclusTPC < 50) continue;  //60 to 50
        //   if(dca2[1]>3) continue;
        //   if(dca2[0]>3) continue;
        myTracks++;
    }
    hNumberOfEvent->Fill(6);
    if (myTracks >= 1) {
        fMultDist5->Fill(myTracks);
        fMultDist3d->Fill(positiveTracks + negativeTracks, positiveTracks, negativeTracks);
    }
    
    if (fDevelopeMode)std::cout << "There are " << myTracks << "  myTracks" << std::endl;
    
    // set Z Vertex bin
    for (Int_t i = 0; i < fZvertexBins; i++) {
        if (fDevelopeMode)std::cout << "fZvertexBins: " << fZvertexBins << std::endl;
        if (fDevelopeMode)std::cout << "loop: " << i << std::endl;
        if ( (primaryVtx[2] > zStart + i * zStep) && (primaryVtx[2] < zStart + (i + 1)*zStep) ) {
            zBin = i;
            break;
        }
    }
    
    // set Multiplicity bin
    for (Int_t i = 0; i < fMultBins; i++) {
        if (fDevelopeMode)std::cout << "Multi bin: " << fMultBins << std::endl;
        if (fDevelopeMode)std::cout << "loop: " << i << std::endl;
        if ( ( myTracks > fMultLimits[i]) && ( myTracks <= fMultLimits[i + 1]) ) {
            mBin = i;
            break;
        }
    }
    if (fDevelopeMode)std::cout << "zBin : " << zBin << "mBin: " << mBin << std::endl;
    
    if (fDevelopeMode)std::cout << "01" << std::endl;
    
    ////////////////////////////////////
    // Add event to buffer if > 0 tracks
    if (myTracks > 0) {
        fEC[zBin][mBin]->FIFOShift();
        (fEvt) = fEC[zBin][mBin]->fEvtStr;
        (fEvt)->fNTracks = myTracks;
        for (Int_t i = 0; i < myTracks; i++) (fEvt)->fTracks[i] = fTempStruct[i];
    }
    
    
    
    if (fMCcase) { // get Input MC information for ESD case
        
        /////////////////////////////////////////////////
        // Xi mc input
        /////////////////////////////////////////////////
        for (Int_t it = 0; it < mcstack->GetNprimary(); it++) {
            TParticle *mcInputTrackXi = (TParticle*)mcstack->Particle(it);
            if (!mcInputTrackXi) {
                Error("UserExec", "Could not receive track %d", it);
                continue;
            }
            
            
            //if(!mcstack->IsPhysicalPrimary(it)) continue;
            if (abs(mcInputTrackXi->GetPdgCode()) != kXiCode) continue;
            
            fXiYDistMC1->Fill(mcInputTrackXi->Y());
            
            
            if (mcInputTrackXi->GetPdgCode() == +kXiCode) {
                hMCinputTotalXi3->Fill(mcInputTrackXi->Pt(),  mcInputTrackXi->Y(), mcInputTrackXi->GetCalcMass());
            }
            else {
                hMCinputTotalXibar3->Fill(mcInputTrackXi->Pt(),  mcInputTrackXi->Y(), mcInputTrackXi->GetCalcMass());
                
            }
            
            
            
        }
        
        
        /////////////////////////////////////////////////
        // XiStar mc input
        /////////////////////////////////////////////////
        for (Int_t it = 0; it < mcstack->GetNprimary(); it++) {
            TParticle *mcInputTrackXiStar = (TParticle*)mcstack->Particle(it);
            if (!mcInputTrackXiStar) {
                Error("UserExec", "Could not receive track %d", it);
                continue;
            }
            if (abs(mcInputTrackXiStar->GetPdgCode()) != kXiStarCode) continue;
            
            
            fXiStarYDistMC1->Fill(mcInputTrackXiStar->Y());
            
            if (mcInputTrackXiStar->GetPdgCode() == +kXiStarCode) {
                hMCinputTotalXiStar3->Fill(mcInputTrackXiStar->Pt(), lPerc, mcInputTrackXiStar->GetCalcMass());
            }
            else {
                hMCinputTotalXiStarbar3->Fill(mcInputTrackXiStar->Pt(), lPerc, mcInputTrackXiStar->GetCalcMass());
                
            }
            
            
        }
    }
    hNumberOfEvent->Fill(7);
    ////////////////////////////////////////////////
    // Reconstruction
    if (fDevelopeMode)std::cout << "02" << std::endl;
    for (Int_t i = 0; i < fESD->GetNumberOfCascades(); i++) {
        
        AliESDcascade *Xicandidate = fESD->GetCascade(i);
        
        
        if (TMath::Abs( Xicandidate->GetPindex()) == TMath::Abs( Xicandidate->GetNindex())) continue;
        if (TMath::Abs( Xicandidate->GetPindex()) == TMath::Abs( Xicandidate->GetBindex())) continue;
        if (TMath::Abs( Xicandidate->GetNindex()) == TMath::Abs( Xicandidate->GetBindex())) continue;
        
        
        
        AliESDtrack *pTrackXi   = fESD->GetTrack(TMath::Abs( Xicandidate->GetPindex()));
        AliESDtrack *nTrackXi   = fESD->GetTrack(TMath::Abs( Xicandidate->GetNindex()));
        AliESDtrack *bTrackXi   = fESD->GetTrack(TMath::Abs( Xicandidate->GetBindex()));
        
        // Standard track QA cuts
        if (!fTrackCut->AcceptTrack(pTrackXi)) continue;
        if (!fTrackCut->AcceptTrack(nTrackXi)) continue;
        if (!fTrackCut->AcceptTrack(bTrackXi)) continue;
        
        
        
        //////////////////////
        // DecayParameters Key (number represents array index)
        // NclustersTPC: 0=proton, 1=pion first, 2=pion second, 3=pion third
        // DCAVtx: 4=proton, 5=pion first, 6=pion second, 7=lambda, 8=pion third
        // 9 = DCA proton-pion
        // 10 = DCA Lambda-pion
        // 11 = Cos PA Lambda
        // 12 = Cos PA Xi
        
        //myxiTracks++;
        
        fDecayParameters[2] = bTrackXi->GetTPCNcls();
        fTPCNcls_pi2->Fill(fDecayParameters[2]);
        
        Double_t fTPCNSigProton = 10;
        Double_t fTPCNSigPion1 = 10;
        Double_t fTPCNSigPion2 = 10;
        
        
        Double_t fTPCPIDMomXi[3] = { -10, -10, -10};
        Double_t fNSigTPCXi[3] = { -10, -10, -10};
        
        
        
        
        if (Xicandidate->Charge() == -1) {
            fDecayParameters[0] = pTrackXi->GetTPCNcls();
            fDecayParameters[1] = nTrackXi->GetTPCNcls();
            fDecayParameters[4] = fabs(pTrackXi->GetD(primaryVtx[0], primaryVtx[1], bField)); // DCA Vtx proton
            fDecayParameters[5] = fabs(nTrackXi->GetD(primaryVtx[0], primaryVtx[1], bField)); // DCA Vtx pion first
            
            fTPCNSigProton = fPIDResponse->NumberOfSigmasTPC(pTrackXi, AliPID::kProton);
            fTPCNSigPion1 = fPIDResponse->NumberOfSigmasTPC(nTrackXi, AliPID::kPion);
            fTPCNSigPion2 = fPIDResponse->NumberOfSigmasTPC(bTrackXi, AliPID::kPion);
            
            fTPCPIDMomXi[0] = pTrackXi->GetTPCmomentum();
            fNSigTPCXi[0] = pTrackXi->GetTPCsignal();
            
            fTPCPIDMomXi[1] = nTrackXi->GetTPCmomentum();
            fNSigTPCXi[1] = nTrackXi->GetTPCsignal();
            
            fTPCPIDMomXi[2] = bTrackXi->GetTPCmomentum();
            fNSigTPCXi[2] = bTrackXi->GetTPCsignal();
            
        } else {
            fDecayParameters[0] = nTrackXi->GetTPCNcls();
            fDecayParameters[1] = pTrackXi->GetTPCNcls();
            fDecayParameters[4] = fabs(nTrackXi->GetD(primaryVtx[0], primaryVtx[1], bField)); // DCA Vtx proton
            fDecayParameters[5] = fabs(pTrackXi->GetD(primaryVtx[0], primaryVtx[1], bField)); // DCA Vtx pion first
            
            fTPCNSigProton = fPIDResponse->NumberOfSigmasTPC(nTrackXi, AliPID::kProton);
            fTPCNSigPion1 = fPIDResponse->NumberOfSigmasTPC(pTrackXi, AliPID::kPion);
            fTPCNSigPion2 = fPIDResponse->NumberOfSigmasTPC(bTrackXi, AliPID::kPion);
            
            
            fTPCPIDMomXi[0] = nTrackXi->GetTPCmomentum();
            fNSigTPCXi[0] = nTrackXi->GetTPCsignal();
            
            fTPCPIDMomXi[1] = pTrackXi->GetTPCmomentum();
            fNSigTPCXi[1] = pTrackXi->GetTPCsignal();
            
            fTPCPIDMomXi[2] = bTrackXi->GetTPCmomentum();
            fNSigTPCXi[2] = bTrackXi->GetTPCsignal();
        }
        
        fTPCNcls_p->Fill(fDecayParameters[0]);
        fTPCNcls_pi1->Fill(fDecayParameters[1]);
        fDCADist_p->Fill(fDecayParameters[4]);
        fDCADist_pi1->Fill(fDecayParameters[5]);
        
        hTPCNSigProton->Fill(fTPCNSigProton);
        hTPCNSigPion1->Fill(fTPCNSigPion1);
        hTPCNSigPion2->Fill(fTPCNSigPion2);
        
        hdEdxProton->Fill(fTPCPIDMomXi[0], fNSigTPCXi[0]);
        hdEdxPion1->Fill(fTPCPIDMomXi[1], fNSigTPCXi[1]);
        hdEdxPion2->Fill(fTPCPIDMomXi[2], fNSigTPCXi[2]);
        
        
        if (fTPCNSigProton > -nSigTPCPID && fTPCNSigProton < nSigTPCPID)hdEdxProtonAfter->Fill(fTPCPIDMomXi[0], fNSigTPCXi[0]);
        if (fTPCNSigPion1 > -nSigTPCPID && fTPCNSigPion1 < nSigTPCPID) hdEdxPion1After->Fill(fTPCPIDMomXi[1], fNSigTPCXi[1]);
        if (fTPCNSigPion2 > -nSigTPCPID && fTPCNSigPion2 < nSigTPCPID)hdEdxPion2After->Fill(fTPCPIDMomXi[2], fNSigTPCXi[2]);
        
        
        // PID Cuts
        if (fPIDOption && abs(fTPCNSigProton) > nSigTPCPID) continue; // PID for proton
        if (fPIDOption && abs(fTPCNSigPion1) > nSigTPCPID) continue; // PID for 1st pion
        if (fPIDOption && abs(fTPCNSigPion2) > nSigTPCPID) continue; // PID for 2nd pion
        
        hQATPCNSigProton->Fill(fTPCNSigProton);
        hQATPCNSigPion1->Fill(fTPCNSigPion1);
        hQATPCNSigPion2->Fill(fTPCNSigPion2);
        
        
        fDecayParameters[6] = fabs(bTrackXi->GetD(primaryVtx[0], primaryVtx[1], bField)); // DCA Vtx pion second
        fDCADist_pi2->Fill(fDecayParameters[6]);
        fDecayParameters[7] = fabs(Xicandidate->GetD(primaryVtx[0], primaryVtx[1], primaryVtx[2])); // DCA Vtx Lambda
        fDCADist_lambda->Fill(fDecayParameters[7]);
        fDecayParameters[9] = fabs(Xicandidate->GetDcaV0Daughters());// DCA proton-pion
        fDCADist_pi_p->Fill(fDecayParameters[9]);
        fDecayParameters[10] = fabs(Xicandidate->GetDcaXiDaughters());// DCA Lambda-pion
        fDCADist_pi_lambda->Fill(fDecayParameters[10]);
        
        Double_t tempX[3] = {0};
        Xicandidate->GetXYZ(tempX[0], tempX[1], tempX[2]);
        
        //    fDecayParameters[11] = sqrt( pow(tempX[0],2) + pow(tempX[1],2));// Rxy Lambda
        //    fRXY_lambda->Fill(fDecayParameters[11]);
        
        fDecayParameters[11] = Xicandidate->GetV0CosineOfPointingAngle(primaryVtx[0], primaryVtx[1], primaryVtx[2]); // Cos PA Lambda
        fCosPA_lambda->Fill(fDecayParameters[11]);
        
        fDecayParameters[12] = Xicandidate->GetCascadeCosineOfPointingAngle(primaryVtx[0], primaryVtx[1], primaryVtx[2]); // Cos PA Xi
        fCosPA_Xi->Fill(fDecayParameters[12]);
        
        decayLengthXY = sqrt( pow(xiVtx[0] - primaryVtx[0], 2) + pow(xiVtx[1] - primaryVtx[1], 2) );
        //     fDecayParameters[12] = decayLengthXY;// Rxy Xi
        //     fRXY_Xi->Fill(fDecayParameters[12]);
        
        
        xiP[0] = Xicandidate->Px();
        xiP[1] = Xicandidate->Py();
        xiP[2] = Xicandidate->Pz();
        xiVtx[0] = Xicandidate->Xv();
        xiVtx[1] = Xicandidate->Yv();
        xiVtx[2] = Xicandidate->Zv();
        xiPt = Xicandidate->Pt();
        xiY = Xicandidate->RapXi();
        xiMass = Xicandidate->M();
        xiCharge = Xicandidate->Charge();
        
        myxiTracks++;
        
        
        
        if (sqrt( pow(tempX[0], 2) + pow(tempX[1], 2) ) > fMaxDecayLength) continue;
        if (decayLengthXY > fMaxDecayLength) continue;
        
        Bool_t StandardXi = kTRUE;
        if (StandardXi)fCutEvents->Fill(1, 1);
        
        if (fDecayParameters[0] < fCutValues[0][0]) StandardXi = kFALSE; // Nclus proton
        if (StandardXi)fCutEvents->Fill(2, 1);
        
        if (fDecayParameters[1] < fCutValues[0][1]) StandardXi = kFALSE; // Nclus pion first
        if (StandardXi)fCutEvents->Fill(3, 1);
        
        if (fDecayParameters[2] < fCutValues[0][2]) StandardXi = kFALSE; // Nclus pion second
        if (StandardXi)fCutEvents->Fill(4, 1);
        
        //
        if (fDecayParameters[4] < fCutValues[0][4]) StandardXi = kFALSE; // DCAVtx proton
        if (StandardXi)fCutEvents->Fill(5, 1);
        
        if (fDecayParameters[5] < fCutValues[0][5]) StandardXi = kFALSE; // DCAVtx pion first
        if (StandardXi)fCutEvents->Fill(6, 1);
        
        if (fDecayParameters[6] < fCutValues[0][6]) StandardXi = kFALSE; // DCAVtx pion second
        if (StandardXi)fCutEvents->Fill(7, 1);
        
        if (fDecayParameters[7] < fCutValues[0][7]) StandardXi = kFALSE; // DCAVtx Lambda
        if (StandardXi)fCutEvents->Fill(8, 1);
        
        //
        if (fDecayParameters[9] > fCutValues[0][9]) StandardXi = kFALSE; // DCAV proton-pion
        if (StandardXi)fCutEvents->Fill(9, 1);
        
        if (fDecayParameters[10] > fCutValues[0][10]) StandardXi = kFALSE; // DCAV Lambda-pion
        if (StandardXi)fCutEvents->Fill(10, 1);
        
        //
        // if(fDecayParameters[11] < fCutValues[0][11]) StandardXi=kFALSE;// Rxy Lambda
        // if(StandardXi)fCutEvents->Fill(11,1);
        
        // if(fDecayParameters[12] < fCutValues[0][12]) StandardXi=kFALSE;// Rxy Xi
        // if(StandardXi)fCutEvents->Fill(12,1);
        
        //
        if (fDecayParameters[11] < fCutValues[0][11]) StandardXi = kFALSE; // Cos PA Lambda
        if (StandardXi)fCutEvents->Fill(11, 1);
        
        if (fDecayParameters[12] < fCutValues[0][12]) StandardXi = kFALSE; // Cos PA Xi
        if (StandardXi)fCutEvents->Fill(12, 1);
        
        if (StandardXi) {
            if (xiCharge == -1) {
                CutVar[0].fXi->Fill(xiPt, xiY, xiMass);
                hXiInvMass->Fill(xiMass);
            }
            else {
                CutVar[0].fXibar->Fill(xiPt, xiY, xiMass);
                hXiInvMass->Fill(xiMass);
            }
        }
        
        
        
        if (fDevelopeMode)std::cout << "001" << std::endl;
        // MC associaton
        mcXiFilled = kFALSE;
        if (fMCcase ) {
            
            MCXiD2esd = (TParticle*)mcstack->Particle(abs(bTrackXi->GetLabel()));
            
            if (abs(MCXiD2esd->GetPdgCode()) == kPionCode) {
                
                MCLamD1esd = (TParticle*)mcstack->Particle(abs(pTrackXi->GetLabel()));
                MCLamD2esd = (TParticle*)mcstack->Particle(abs(nTrackXi->GetLabel()));
                
                if (MCLamD1esd->GetMother(0) == MCLamD2esd->GetMother(0)) {
                    if (abs(MCLamD1esd->GetPdgCode()) == kProtonCode || abs(MCLamD2esd->GetPdgCode()) == kProtonCode) {
                        if (abs(MCLamD1esd->GetPdgCode()) == kPionCode || abs(MCLamD2esd->GetPdgCode()) == kPionCode) {
                            
                            MCLamesd = (TParticle*)mcstack->Particle(abs(MCLamD1esd->GetMother(0)));
                            if (abs(MCLamesd->GetPdgCode()) == kLambdaCode) {
                                
                                if (MCLamesd->GetMother(0) == MCXiD2esd->GetMother(0)) {
                                    MCXiesd = (TParticle*)mcstack->Particle(abs(MCLamesd->GetMother(0)));
                                    if (abs(MCXiesd->GetPdgCode()) == kXiCode) {
                                        mcXiFilled = kTRUE;
                                        
                                        if (StandardXi) {
                                            
                                            fXiYDistMCout->Fill(xiY);
                                            
                                            if (Xicandidate->Charge() == -1) {
                                                
                                                CutVar[0].fMCrecXi->Fill(xiPt, xiY, xiMass);
                                            } else {
                                                CutVar[0].fMCrecXibar->Fill(xiPt, xiY, xiMass);
                                            }
                                            
                                            
                                            
                                        }
                                        
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }// MC association
        
        
        if (fabs(xiMass - fTrueMassXi) > fMassWindow) continue;
        
        if (StandardXi)hQAXiInvMass->Fill(xiMass);
        
        fXiTrack->Set(xiVtx, xiP, fCovMatrix, Short_t(xiCharge));
        
        if (fDevelopeMode)std::cout << "002" << std::endl;
        
        //////////////////////////////////////////////////////////
        // Reconstruct Xi(1530)
        for (Int_t EN = 0; EN < fEventsToMix + 1; EN++) { // Event buffer loop
            
            for (Int_t l = 0; l < (fEvt + EN)->fNTracks; l++) { // Present(EN=0) and Past(EN from 1 to fEventsToMix) event track loop
                
                if (EN == 0) {
                    if ((fEvt + EN)->fTracks[l].fID == pTrackXi->GetID()) continue;
                    if ((fEvt + EN)->fTracks[l].fID == nTrackXi->GetID()) continue;
                    if ((fEvt + EN)->fTracks[l].fID == bTrackXi->GetID()) continue;
                }
                
                fXiTrack->Set(xiVtx, xiP, fCovMatrix, Short_t(xiCharge));
                
                if (!fESDTrack4) continue;
                fESDTrack4->Set((fEvt + EN)->fTracks[l].fX, (fEvt + EN)->fTracks[l].fP, (fEvt + EN)->fTracks[l].fCov, (fEvt + EN)->fTracks[l].fCharge);
                
                fDecayParameters[8] = (fEvt + EN)->fTracks[l].fDCAXY; // DCA Vtx pion third
                fDCADist_3rd_pi->Fill(fDecayParameters[8]);
                
                
                
                if ((fEvt + EN)->fTracks[l].fDCAZ > 2) continue;
                if ( (((fEvt + EN)->fTracks[l].fStatus)&AliESDtrack::kITSrefit) == 0) continue; // Require itsrefit
                // no Chi^2 cut applied for ESDs.  Info not available in my track structure.
                
                
                if (fabs((fEvt + EN)->fTracks[l].fEta) > 0.8) continue;
                
                fDecayParameters[3] = (fEvt + EN)->fTracks[l].fNclusTPC;
                fTPCNcls_pi3->Fill(fDecayParameters[3]);
                
                AliVertex *XiStarVtx = new AliVertex((fEvt + EN)->fTracks[l].fX, 0, 0);
                //fESDTrack4->PropagateToDCA(fXiTrack, bField);// Propagate tracks to dca, both tracks are budged
                if (!(fXiTrack->PropagateToDCA(XiStarVtx, bField, 3))) continue; // Propagate tracks to dca, version which assumes fESDTrack4 is already primary
                /////////////
                fXiTrack->GetPxPyPz(pDaughter1);
                fXiTrack->GetXYZ(xDaughter1);
                fESDTrack4->GetPxPyPz(pDaughter2);
                fESDTrack4->GetXYZ(xDaughter2);
                //////////////////////////
                
                
                
                //xiStarVtx[0] = (xDaughter1[0]+xDaughter2[0])/2.;
                //xiStarVtx[1] = (xDaughter1[1]+xDaughter2[1])/2.;
                //xiStarVtx[2] = (xDaughter1[2]+xDaughter2[2])/2.;
                //decayLength = sqrt(pow(xiStarVtx[0]-primaryVtx[0],2)+pow(xiStarVtx[1]-primaryVtx[1],2)+pow(xiStarVtx[2]-primaryVtx[2],2));
                
                px1 = pDaughter1[0];
                py1 = pDaughter1[1];
                pz1 = pDaughter1[2];
                px2 = pDaughter2[0];
                py2 = pDaughter2[1];
                pz2 = pDaughter2[2];
                
                p1sq = px1 * px1 + py1 * py1 + pz1 * pz1;
                p2sq = px2 * px2 + py2 * py2 + pz2 * pz2;
                if (p1sq <= 0 || p2sq <= 0) continue;
                
                e1 = sqrt(p1sq + fTrueMassXi * fTrueMassXi);
                e2 = sqrt(p2sq + fTrueMassPi * fTrueMassPi);
                angle = px1 * px2 + py1 * py2 + pz1 * pz2;
                xiStarMass = fTrueMassXi * fTrueMassXi + fTrueMassPi * fTrueMassPi + 2.*e1 * e2 - 2.*angle;
                if (xiStarMass < 0.) xiStarMass = 1.e-8;
                xiStarMass = sqrt(xiStarMass);
                
                
                xiStarP[0] = px1 + px2;
                xiStarP[1] = py1 + py2;
                xiStarP[2] = pz1 + pz2;
                xiStarMom = sqrt(pow(xiStarP[0], 2) + pow(xiStarP[1], 2) + pow(xiStarP[2], 2));
                if (xiStarMom == 0) continue; // So one of the following lines doesnt break
                xiStarPt = sqrt(xiStarP[0] * xiStarP[0] + xiStarP[1] * xiStarP[1]);
                
                
                xiStarY = .5 * log( ((e1 + e2) + xiStarP[2]) / ((e1 + e2) - xiStarP[2]));
                fXiStarYDist->Fill(xiStarY);
                if (xiStarY < -0.5 || xiStarY > 0.5) continue; // here selection of rapidity for pPb analysis
                
                fQAXiStarYDist->Fill(xiStarY);
                
                for (int cv = 0; cv < kNCutVariations; cv++) {
                    if (!fSetSystematic && cv > 0) continue;
                    
                    if (fDecayParameters[0] < fCutValues[cv][0]) continue; // Nclus proton
                    if (fDecayParameters[1] < fCutValues[cv][1]) continue; // Nclus pion first
                    if (fDecayParameters[2] < fCutValues[cv][2]) continue; // Nclus pion second
                    if (fDecayParameters[3] < fCutValues[cv][3]) continue; // Nclus pion third
                    //
                    if (fDecayParameters[4] < fCutValues[cv][4]) continue; // DCAVtx proton
                    if (fDecayParameters[5] < fCutValues[cv][5]) continue; // DCAVtx pion first
                    if (fDecayParameters[6] < fCutValues[cv][6]) continue; // DCAVtx pion second
                    if (fDecayParameters[7] < fCutValues[cv][7]) continue; // DCAVtx Lambda
                    //if(fDecayParameters[8] > fCutValues[cv][8]) continue; // DCAVtx pion third
                    if (fDecayParameters[8] > (0.0105 + 0.035 / pow((fEvt + EN)->fTracks[l].fPt, 1.1))) continue; // DCAVtx pion third
                    //0.0182 + 0.035/pow((fEvt+EN)->fTracks[l].fPt,1.1 (2010 cut)
                    
                    //
                    if (fDecayParameters[9] > fCutValues[cv][9]) continue; // DCAV proton-pion
                    if (fDecayParameters[10] > fCutValues[cv][10]) continue; // DCAV Lambda-pion
                    //
                    if (fDecayParameters[11] < fCutValues[cv][11]) continue; // Cos PA Lambda
                    if (fDecayParameters[12] < fCutValues[cv][12]) continue; // Cos PA Xi
                    
                    if (EN == 0 && cv == 0) { // cut QA plot for default cut
                        fQATPCNcls_p->Fill(fDecayParameters[0]);
                        fQATPCNcls_pi1->Fill(fDecayParameters[1]);
                        fQATPCNcls_pi2->Fill(fDecayParameters[2]);
                        fQATPCNcls_pi3->Fill(fDecayParameters[3]);
                        
                        fQADCADist_p->Fill(fDecayParameters[4]);
                        fQADCADist_pi1->Fill(fDecayParameters[5]);
                        fQADCADist_pi2->Fill(fDecayParameters[6]);
                        
                        fQADCADist_lambda->Fill(fDecayParameters[7]);
                        fQADCADist_3rd_pi->Fill(fDecayParameters[8]);
                        fQADCADist_pi_p->Fill(fDecayParameters[9]);
                        fQADCADist_pi_lambda->Fill(fDecayParameters[10]);
                        fQACosPA_lambda->Fill(fDecayParameters[11]);
                        fQACosPA_Xi->Fill(fDecayParameters[12]);
                    }
                    
                    if (EN == 0 && cv == 1) {
                        fQATPCNcls_p_L->Fill(fDecayParameters[0]);
                        fQATPCNcls_pi1_L->Fill(fDecayParameters[1]);
                        fQATPCNcls_pi2_L->Fill(fDecayParameters[2]);
                        fQATPCNcls_pi3_L->Fill(fDecayParameters[3]);
                    }
                    if (EN == 0) {
                        
                        if (cv == 2)fQADCADist_p_L->Fill(fDecayParameters[4]);
                        if (cv == 3)fQADCADist_pi1_L->Fill(fDecayParameters[5]);
                        if (cv == 4)fQADCADist_pi2_L->Fill(fDecayParameters[6]);
                        if (cv == 5)fQADCADist_lambda_L->Fill(fDecayParameters[7]);
                        if (cv == 6)fQADCADist_3rd_pi_L->Fill(fDecayParameters[8]);
                        if (cv == 7)fQADCADist_pi_p_L->Fill(fDecayParameters[9]);
                        if (cv == 8)fQADCADist_pi_lambda_L->Fill(fDecayParameters[10]);
                        if (cv == 9)fQACosPA_lambda_L->Fill(fDecayParameters[11]);
                        if (cv == 10)fQACosPA_Xi_L->Fill(fDecayParameters[12]);
                        
                    }
                    
                    if (EN == 0 && cv == 11) {
                        fQATPCNcls_p_T->Fill(fDecayParameters[0]);
                        fQATPCNcls_pi1_T->Fill(fDecayParameters[1]);
                        fQATPCNcls_pi2_T->Fill(fDecayParameters[2]);
                        fQATPCNcls_pi3_T->Fill(fDecayParameters[3]);
                    }
                    if (EN == 0) {
                        
                        if (cv == 12)fQADCADist_p_T->Fill(fDecayParameters[4]);
                        if (cv == 13)fQADCADist_pi1_T->Fill(fDecayParameters[5]);
                        if (cv == 14)fQADCADist_pi2_T->Fill(fDecayParameters[6]);
                        if (cv == 15)fQADCADist_lambda_T->Fill(fDecayParameters[7]);
                        if (cv == 16)fQADCADist_3rd_pi_T->Fill(fDecayParameters[8]);
                        if (cv == 17)fQADCADist_pi_p_T->Fill(fDecayParameters[9]);
                        if (cv == 18)fQADCADist_pi_lambda_T->Fill(fDecayParameters[10]);
                        if (cv == 19)fQACosPA_lambda_T->Fill(fDecayParameters[11]);
                        if (cv == 20)fQACosPA_Xi_T->Fill(fDecayParameters[12]);
                        
                        
                    }
                    
                    
                    if (EN == 0) {
                        if (fXiTrack->Charge() == -1 &&  fESDTrack4->Charge() == -1) CutVar[cv].fXiMinusPiMinus->Fill(xiStarPt, lPerc, xiStarMass);
                        else if (fXiTrack->Charge() == -1 &&  fESDTrack4->Charge() == +1) {
                            // xiStarPt, xiStarY, xiStarMass -> xiStarPt, lPerc, xiStarMass
                            CutVar[cv].fXiMinusPiPlus->Fill(xiStarPt, lPerc, xiStarMass);
                        }
                        else if (fXiTrack->Charge() == +1 &&  fESDTrack4->Charge() == -1) {
                            CutVar[cv].fXiPlusPiMinus->Fill(xiStarPt, lPerc, xiStarMass);
                        }
                        else CutVar[cv].fXiPlusPiPlus->Fill(xiStarPt, lPerc, xiStarMass);
                    } else {
                        if (fXiTrack->Charge() == -1 &&  fESDTrack4->Charge() == -1) CutVar[cv].fXiMinusPiMinusbkg->Fill(xiStarPt, lPerc, xiStarMass);
                        else if (fXiTrack->Charge() == -1 &&  fESDTrack4->Charge() == +1) CutVar[cv].fXiMinusPiPlusbkg->Fill(xiStarPt, lPerc, xiStarMass);
                        else if (fXiTrack->Charge() == +1 &&  fESDTrack4->Charge() == -1) CutVar[cv].fXiPlusPiMinusbkg->Fill(xiStarPt, lPerc, xiStarMass);
                        else CutVar[cv].fXiPlusPiPlusbkg->Fill(xiStarPt, lPerc, xiStarMass);
                    }
                    
                    
                    
                    // MC associaton ESD
                    if (fMCcase && mcXiFilled && EN == 0) { // ESD MC's
                        MCXiStarD2esd = (TParticle*)mcstack->Particle(abs((fEvt)->fTracks[l].fLabel));
                        
                        if (abs(MCXiStarD2esd->GetPdgCode()) == kPionCode) {
                            if (MCXiesd->GetMother(0) == MCXiStarD2esd->GetMother(0)) {
                                
                                MCXiStaresd = (TParticle*)mcstack->Particle(abs(MCXiesd->GetMother(0)));
                                if (abs(MCXiStaresd->GetPdgCode()) == kXiStarCode) {
                                    
                                    fXiStarYDistMCout->Fill(xiStarY);
                                    
                                    
                                    if (fXiTrack->Charge() == -1 &&  fESDTrack4->Charge() == +1) CutVar[cv].fMCrecXiMinusPiPlus->Fill(xiStarPt, lPerc, xiStarMass);
                                    if (fXiTrack->Charge() == +1 &&  fESDTrack4->Charge() == -1) CutVar[cv].fMCrecXiPlusPiMinus->Fill(xiStarPt, lPerc, xiStarMass);
                                }
                            }
                        }
                    }
                }// Cut Variation loop
            }// 3rd pion loop
        }// Event mixing loop
    }// Xi loop
    // Post output data.
    PostData(1, fOutputList);
    */
}
//________________________________________________________________________

void AliAnalysisTaskXi1530::Terminate(Option_t *)
{
    //if(fDevelopeMode)std::cout<<"Done"<<std::endl;
}
//________________________________________________________________________
/*
 Double_t AliAnalysisTaskXi1530::LinearPropagateToDCA(AliESDtrack *v, AliESDtrack *t, Double_t b) {// Adapted from AliCascadeVertexer.cxx
 //--------------------------------------------------------------------
 // This function returns the DCA between the V0 and the track
 //--------------------------------------------------------------------
 Double_t alpha=t->GetAlpha(), cs1=TMath::Cos(alpha), sn1=TMath::Sin(alpha);
 Double_t r[3]; t->GetXYZ(r);
 Double_t x1=r[0], y1=r[1], z1=r[2];
 Double_t p[3]; t->GetPxPyPz(p);
 Double_t px1=p[0], py1=p[1], pz1=p[2];
 Double_t x2[3]={0};
 Double_t p2[3]={0};
 Double_t vx2,vy2,vz2;     // position and momentum of V0
 Double_t px2,py2,pz2;
 v->GetXYZ(x2);
 v->GetPxPyPz(p2);
 vx2=x2[0], vy2=x2[1], vz2=x2[2];
 px2=p2[0], py2=p2[1], pz2=p2[2];
 // calculation dca
 Double_t dd= Det(vx2-x1,vy2-y1,vz2-z1,px1,py1,pz1,px2,py2,pz2);
 Double_t ax= Det(py1,pz1,py2,pz2);
 Double_t ay=-Det(px1,pz1,px2,pz2);
 Double_t az= Det(px1,py1,px2,py2);
 Double_t dca=TMath::Abs(dd)/TMath::Sqrt(ax*ax + ay*ay + az*az);
 //points of the DCA
 Double_t t1 = Det(vx2-x1,vy2-y1,vz2-z1,px2,py2,pz2,ax,ay,az)/
 Det(px1,py1,pz1,px2,py2,pz2,ax,ay,az);
 x1 += px1*t1; y1 += py1*t1; //z1 += pz1*t1;
 //propagate track to the points of DCA
 x1=x1*cs1 + y1*sn1;
 if (!t->PropagateTo(x1,b)) {
 Error("PropagateToDCA","Propagation failed !");
 return 1.e+33;
 }
 return dca;
 }
 */

//________________________________________________________________________

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

THnSparse * AliAnalysisTaskXi1530::CreateTHnSparse(TString name
                                                   , TString title, Int_t ndim, std::vector<TAxis> bins, Option_t * opt){
    const TAxis * axises[bins.size()];
    for( UInt_t i=0;i<bins.size();i++ ) axises[i]= &bins[i];
    THnSparse * h= fHistos->CreateTHnSparse(name, title, ndim, axises,opt );
    return h;
}

Long64_t AliAnalysisTaskXi1530::FillTHnSparse( TString name, std::vector<Double_t> x, Double_t w ){
    auto hsparse = dynamic_cast<THnSparse*>( fHistos->FindObject(name) );
    if(! hsparse ){
        cout<<"ERROR : no "<<name<<endl;
        exit(1);
    }
    return FillTHnSparse( hsparse, x, w );
}

Long64_t AliAnalysisTaskXi1530::FillTHnSparse( THnSparse *h, std::vector<Double_t> x, Double_t w ){
    if( int(x.size()) != h->GetNdimensions() ){
        cout<<"ERROR : wrong sized of array while Fill "<<h->GetName()<<endl;
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

