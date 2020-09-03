#if !defined(__CINT__) || defined(__CLING__)
#include "AliAODInputHandler.h"
#include "AliAnalysisAlien.h"
#include "AliPhysicsSelectionTask.h"
#include "AliMultSelectionTask.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliAnalysisTaskXi1530.h"
#include "AliAnalysisTaskXi1530temp.h"
#include "AliESDInputHandler.h"
#endif

void run_KIAF(const char* dataset = "test1.list",
              const char* option = "AOD_Mix_Nano_SYS") {
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
  bool isDev = kFALSE;
  bool isNano = kFALSE;
  bool setmixing = kFALSE;
  bool isaod = kFALSE;
  int nmix = 10;
  bool highmult = kFALSE;
  TString foption = option;
  const char* taskname = "Xi1530MB";
  const char* suffix = "_";
  if (foption.Contains("MC"))
    ismc = kTRUE;
  if (foption.Contains("Dev"))
    isDev = kTRUE;
  if (foption.Contains("AOD"))
    isaod = true;
  if (foption.Contains("Nano"))
    isNano = true;

  // analysis manager
  AliAnalysisManager* mgr =
      new AliAnalysisManager(Form("%s%s", taskname, option));
  AliInputEventHandler* handler;
  if (isaod)
    handler = new AliAODInputHandler();
  else
    handler = new AliESDInputHandler();

  // handler->SetNeedField(1);
  mgr->SetInputEventHandler(handler);

  if (ismc) {
    AliMCEventHandler* mcHandler = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mcHandler);
  }
  if(!isNano){
    //
    // Physics Selection
    AliPhysicsSelectionTask* physSelTask =
        reinterpret_cast<AliPhysicsSelectionTask*>(gInterpreter->ExecuteMacro(
            Form("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C(%d)",
                ismc)));
    // Multiplicity selection
    AliMultSelectionTask* MultSlection =
        reinterpret_cast<AliMultSelectionTask*>(gInterpreter->ExecuteMacro(
            "$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/"
            "AddTaskMultSelection.C"));
    MultSlection->SetAddInfo(kTRUE);
    MultSlection->SetSelectedTriggerClass(AliVEvent::kAny);
    // MultSlection->SetAlternateOADBforEstimators("LHC16k"); //if needed
    // PID response
    AliAnalysisTask* fPIDResponse =
        reinterpret_cast<AliAnalysisTask*>(gInterpreter->ExecuteMacro(
            Form("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C(%d)", ismc)));
  }
  if(isNano)
    gInterpreter->ExecuteMacro("$ALICE_PHYSICS/PWG/DevNanoAOD/macros/AddTaskNanoAODnormalisation.C");
  if (isDev) {
    gInterpreter->LoadMacro("AliAnalysisTaskXi1530temp.cxx+g");
    AliAnalysisTaskXi1530temp* task =
        reinterpret_cast<AliAnalysisTaskXi1530temp*>(gInterpreter->ExecuteMacro(
            Form("AddTaskXi1530.C(\"%s\",\"%s\",%i)", taskname, option, nmix)));

  } else {
      AliAnalysisTaskXi1530* task =
        reinterpret_cast<AliAnalysisTaskXi1530*>(gInterpreter->ExecuteMacro(
            Form("$ALICE_PHYSICS/PWGLF/RESONANCES/PostProcessing/Xi1530/AddTaskXi1530.C(\"%s\",\"%s\",%i)", taskname, option, nmix)));
  }
  if (!mgr->InitAnalysis())
    return;
  mgr->PrintStatus();

  // start analysis
  Printf("Starting Analysis....");

  // TChain* chain = CreateESDChain(dataset,-1);
  TChain* chain = new TChain("ESDTree");
  std::stringstream esdChain;
  if(isNano){
    esdChain << ".x " << "CreateNanoAODChain.C(";
    esdChain << "\"" << dataset << "\", -1);";
  }
  else if (isaod) {
    esdChain << ".x " << gSystem->Getenv("ALICE_PHYSICS")
             << "/PWG/EMCAL/macros/CreateAODChain.C(";
    esdChain << "\"" << dataset << "\", -1);";
  } else {
    esdChain << ".x " << gSystem->Getenv("ALICE_PHYSICS")
             << "/PWG/EMCAL/macros/CreateESDChain.C(";
    esdChain << "\"" << dataset << "\", -1);";
  }
  chain = reinterpret_cast<TChain*>(gROOT->ProcessLine(esdChain.str().c_str()));
  mgr->StartAnalysis("local", chain);
}
