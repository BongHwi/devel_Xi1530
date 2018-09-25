#ifndef ALIANALYSISTASKXi1530_H
#define ALIANALYSISTASKXi1530_H
//
// Class AliAnalysisTaskXi1530
//
// AliAnalysisTaskXi1530
// author:
//  (Original Code) Beomkyu KIM ()
//  (Modification) Bong-Hwi Lim (bong-hwi.lim@cern.ch)

#include "THnSparse.h"
#include "AliAnalysisTaskSE.h"
#include "AliTriggerAnalysis.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliStack.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "THistManager.h"
#include "AliMultiplicitySelectionCPPWA.h"
#include <deque>

class AliAnalysisTaskXi1530RunTable {
public:
    enum {kPP,kPA,kAA,kUnknownCollType};
    AliAnalysisTaskXi1530RunTable();
    AliAnalysisTaskXi1530RunTable(Int_t runnumber);
    ~AliAnalysisTaskXi1530RunTable();
    
    Bool_t IsPP(){
        return fCollisionType==kPP;
    }
    Bool_t IsPA(){
        return fCollisionType==kPA;
    }
    Bool_t IsAA(){
        return fCollisionType==kAA;
    }
private:
    Int_t  fCollisionType; //! Is proton-proton collisions?
};

class AliAnalysisTaskXi1530 : public AliAnalysisTaskSE {
    public:
        enum {kPP,kPA,kAA,kUnknownCollType};
    
        AliAnalysisTaskXi1530();
        AliAnalysisTaskXi1530(const char *name, const char *option);
        AliAnalysisTaskXi1530(const AliAnalysisTaskXi1530& ap);
        AliAnalysisTaskXi1530& operator =(const AliAnalysisTaskXi1530& ap);
        ~AliAnalysisTaskXi1530();
    
        virtual void    UserCreateOutputObjects();
        virtual void    UserExec(Option_t *);
        virtual void    FinishTaskOutput();
        virtual void    Terminate(Option_t *);
    
        void SetOption(char * option) {fOption = option;}
        void SetFilterBit(UInt_t filterbit) {fFilterBit = filterbit;}
        Int_t GetPID(AliPIDResponse *pid, const AliVTrack *trk);
    
        Bool_t  GoodTracksSelection();
        void FillTracks();
        void SetMixing (Bool_t setmixing) {fsetmixing = setmixing;}
        void SetIsAA (Bool_t isaa) {IsAA = isaa;}
        void SetIsMC (Bool_t ismc) {IsMC = ismc;}
        void SetParticleType(Int_t partype) {fParticleType = partype;}
    
        TAxis AxisFix( TString name, int nbin, Double_t xmin, Double_t xmax);
        TAxis AxisVar( TString name, std::vector<Double_t> bin );
        TAxis AxisLog( TString name, int nbin, Double_t xmin, Double_t xmax
                      , Double_t xmin0);
        TAxis AxisStr( TString name, std::vector<TString> bin );
        THnSparse * CreateTHnSparse(TString name, TString title
                                    , Int_t ndim, std::vector<TAxis> bins, Option_t * opt="");
        THnSparse * CreateTHnSparse(TString name, TString title
                                    , TString templ, Option_t * opt="");
        Long64_t FillTHnSparse( TString name, std::vector<Double_t> x, Double_t w=1.);
        Long64_t FillTHnSparse( THnSparse *h, std::vector<Double_t> x, Double_t w=1.);
    
    private:
        typedef std::vector<AliVTrack*> tracklist;
        typedef std::deque<tracklist>  eventpool;
        typedef std::vector<vector<eventpool> > mixingpool;
    
        TString                         fOption;
        TList*                          fOutput=nullptr; //!
    
        AliTriggerAnalysis*             fTrigger=nullptr; //!
        AliESDtrackCuts*                fTrackCuts=nullptr; //!
        AliVEvent*                      fEvt=nullptr; //!
        UInt_t                          fFilterBit;
        Bool_t                          IsFirstEvent=kTRUE;
        AliAnalysisTaskXi1530RunTable*   fRunTable=nullptr; //!
    
        Double_t                        fCent=-1;
        Double_t                        fZ=-30;
        std::vector < UInt_t >          goodtrackindices; //!
    
        AliPIDResponse                 *fPIDResponse=nullptr; //!
        AliPIDCombined                 *fPIDCombined=nullptr; //!
        AliMultiplicitySelectionCPPWA  *fSelec=nullptr; //!
        //Histograms below are main
        std::vector< std::vector< TH2D* > > fMass2D; //! signbins, centbins
        //Histograms for pT_pair amd pT
    
        mixingpool                      fEMpool; //!
        TAxis                           binCent; //!
        TAxis                           binZ; //!
        Int_t                           centbin = -1 ;
        Int_t                           zbin = -1 ;
        Double_t                        fptcut = 0.15;
        Double_t                        fetacut = 0.9;
        Bool_t                          fsetmixing = kFALSE;
        Bool_t                          IsDGV0=kFALSE;
        Bool_t                          IsDGV0FMD=kFALSE;
        Bool_t                          IsAA=kFALSE;
        Bool_t                          IsMC=kFALSE;
        THistManager*                   fHistos=nullptr; //!
        TClonesArray*                   fMCArray=nullptr; //!
        AliStack*                       fMCStack=nullptr; //!
        Int_t                           fParticleType;
        Int_t                           fNTracks = 0;
    
    

    ClassDef(AliAnalysisTaskXi1530, 1);
};

#endif
