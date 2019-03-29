#if !defined(__CINT__) || defined(__CLING__)
#include "AliAODInputHandler.h"
#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskXi1530temp.h"
#include "AliESDInputHandler.h"

#include "./CreateESDChain.C"
#endif

void run_KIAF(const char* dataset = "test1.list",
        const char* taskname = "Xi1530",
         const char* option =
             "LHC16k_pass2_SYS_Mix"  // when scanning AOD, add "AOD"
         ,
         const char* gridmode = "test"  // or "terminate" to merge
         ,
         UInt_t istart = 0,
         UInt_t iend = 25,
         const char* localorgrid = "local") {
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libSTEERBase. so");
    gSystem->Load("libESD.so");
    gSystem->Load("libANALYSIS.so");
    gSystem->Load("libOADB.so");
    gSystem->Load("libANALYSISalice.so");
    gSystem->Load("libpythia6_4_21.so");

    gSystem->Setenv("alien_CLOSE_SE", "working_disk_SE");
    // ROOT 6 MODE
    // add aliroot indlude path
    gInterpreter->ProcessLine(
        Form(".include %s/include", gSystem->ExpandPathName("$ROOTSYS")));
    gInterpreter->ProcessLine(
        Form(".include %s/include", gSystem->ExpandPathName("$ALICE_ROOT")));
    gInterpreter->ProcessLine(
        Form(".include %s/include", gSystem->ExpandPathName("$ALICE_PHYSICS")));

    bool isaa = kFALSE;
    bool ismc = kFALSE;
    bool setmixing = kFALSE;
    bool vertexer = false;
    int nmix = 20;
    bool highmult = kFALSE;
    TString foption = option;
    const char* suffix = "test";
    if (foption.Contains("MC"))
        ismc = kTRUE;
    if (foption.Contains("Vertex"))
        vertexer = true;

    // analysis manager
    AliAnalysisManager* mgr =
        new AliAnalysisManager(Form("%s%s", taskname, option));
    AliInputEventHandler* handler;
    if (foption.Contains("AOD"))
        handler = new AliAODInputHandler();
    else
        handler = new AliESDInputHandler();

    // handler->SetNeedField(1);
    mgr->SetInputEventHandler(handler);

    if (ismc) {
        AliMCEventHandler* mcHandler = new AliMCEventHandler();
        mgr->SetMCtruthEventHandler(mcHandler);
    }
    //
    // Physics Selection
    AliPhysicsSelectionTask* physSelTask =
        reinterpret_cast<AliPhysicsSelectionTask*>(gInterpreter->ExecuteMacro(
            Form("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C(%d)",
                 ismc)));
    if (!physSelTask) {
        Printf("no physSelTask");
        return;
    }
    // Multiplicity selection
    AliMultSelectionTask* MultSlection =
        reinterpret_cast<AliMultSelectionTask*>(gInterpreter->ExecuteMacro(
            "$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/"
            "AddTaskMultSelection.C"));
    if (!MultSlection) {
        Printf("no MultSlection");
        return;
    }
    MultSlection->SetAddInfo(kTRUE);
    MultSlection->SetSelectedTriggerClass(AliVEvent::kAny);
    // MultSlection->SetAlternateOADBforEstimators("LHC16k"); //if needed
    // PID response
    AliAnalysisTask* fPIDResponse =
        reinterpret_cast<AliAnalysisTaskXi1530temp*>(
            gInterpreter->ExecuteMacro(Form(
                "$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C(%d)", ismc)));
    if (!fPIDResponse) {
        Printf("no fPIDResponse");
        return;
    }
    // V0, Xi Super verexter by David
    if(vertexer){
        AliAnalysisTaskWeakDecayVertexer* taskWDV =
            reinterpret_cast<AliAnalysisTaskWeakDecayVertexer*>(
                gInterpreter->ExecuteMacro(
                    "$ALICE_PHYSICS/PWGLF/STRANGENESS/Cascades/Run2/macros/"
                    "AddTaskWeakDecayVertexer.C"));
        taskWDV->SetUseImprovedFinding();

        //______________________________________________________________
        // Revertexing configuration
        // WARNING: This applies only to the Cascade analysis
        taskWDV->SetUseImprovedFinding();

        // V0-Related topological selections
        taskWDV->SetV0VertexerDCAFirstToPV(0.05);
        taskWDV->SetV0VertexerDCASecondtoPV(0.05);
        taskWDV->SetV0VertexerDCAV0Daughters(1.6);
        taskWDV->SetV0VertexerCosinePA(0.97);
        taskWDV->SetV0VertexerMinRadius(0.5);
        taskWDV->SetV0VertexerMaxRadius(200);

        // Cascade-Related topological selections
        taskWDV->SetCascVertexerMinV0ImpactParameter(0.05);
        taskWDV->SetCascVertexerV0MassWindow(0.007);
        taskWDV->SetCascVertexerDCABachToPV(0.05);
        taskWDV->SetCascVertexerDCACascadeDaughters(1.6);
        taskWDV->SetCascVertexerCascadeMinRadius(.5);
        taskWDV->SetCascVertexerCascadeCosinePA(.97);
    }
    gInterpreter->LoadMacro("AliAnalysisTaskXi1530temp.cxx+g");
    // AliAnalysisTaskXi1530temp *myTask =
    // reinterpret_cast<AliAnalysisTaskXi1530temp*>(gInterpreter->ExecuteMacro(Form("AddTaskXi1530.c(\"%s\",\"%s\",%i,%d,%d,%d,%d)",taskname,option,nmix,highmult,isaa,ismc,setmixing)));

    AliAnalysisTaskXi1530temp* myTask =
        reinterpret_cast<AliAnalysisTaskXi1530temp*>(gInterpreter->ExecuteMacro(
            Form("AddTaskXi1530.C(\"%s\",\"%s\",%i,\"%s\")", taskname, option,
                 nmix, suffix)));
    // mgr->SetDebugLevel(1);
    if (!mgr->InitAnalysis())
        return;
    mgr->PrintStatus();

    // start analysis
    Printf("Starting Analysis....");

    TChain* chain = CreateESDChain(dataset);
    mgr->StartAnalysis("local",chain);
}
