#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisTaskXi1530.h"
#endif


#if !defined (__CINT__) || defined (__CLING__)
// ROOT 6 MODE
//vector<Int_t> LHC16kRuns = {257605}; // for test
vector<Int_t> LHC16kRuns = {258537, 258499, 258477, 258456, 258454, 258452, 258426, 258393, 258391, 258387, 258359, 258336, 258332, 258307, 258306, 258303, 258302, 258301, 258299, 258278, 258274, 258273, 258271, 258270, 258258, 258257, 258256, 258204, 258203, 258202, 258198, 258197, 258178, 258117, 258114, 258113, 258109, 258108, 258107, 258063, 258062, 258060, 258059, 258053, 258049, 258045, 258042, 258041, 258039, 258019, 258017, 258014, 258012, 258008, 258003, 257992, 257989, 257986, 257979, 257963, 257960, 257957, 257939, 257937, 257936, 257892, 257855, 257853, 257851, 257850, 257804, 257803, 257800, 257799, 257798, 257797, 257773, 257765, 257757, 257754, 257737, 257735, 257734, 257733, 257727, 257725, 257724, 257697, 257694, 257692, 257691, 257689, 257688, 257687, 257685, 257684, 257682, 257644, 257642, 257636, 257635, 257632, 257630, 257606, 257605, 257604, 257601, 257595, 257594, 257592, 257590, 257588, 257587, 257566, 257562, 257561, 257560, 257541, 257540, 257539, 257537, 257531, 257530, 257492, 257491, 257490, 257488, 257487, 257474, 257468, 257457, 257433, 257364, 257358, 257330, 257322, 257320, 257318, 257260, 257224, 257209, 257206, 257204, 257144, 257141, 257139, 257138, 257137, 257136, 257100, 257095, 257092, 257086, 257084, 257082, 257080, 257077, 257028, 257026, 257021, 257012, 257011, 256944, 256942, 256941, 256697, 256695, 256694, 256692, 256691, 256684, 256681, 256677, 256676, 256658, 256620, 256619, 256592, 256591, 256589, 256567, 256565, 256564, 256562, 256560, 256557, 256556, 256554, 256552, 256514, 256512, 256510, 256506, 256504}; // for run
//const int LHC16kRuns[]={258537, 258499, 258477, 258456, 258454, 258452, 258426, 258393, 258391, 258387, 258359, 258336, 258332, 258307, 258306, 258303, 258302, 258301, 258299, 258278, 258274, 258273, 258271, 258270, 258258, 258257, 258256, 258204, 258203, 258202, 258198, 258197, 258178, 258117, 258114, 258113, 258109, 258108, 258107, 258063, 258062, 258060, 258059, 258053, 258049, 258045, 258042, 258041, 258039, 258019, 258017, 258014, 258012, 258008, 258003, 257992, 257989, 257986, 257979, 257963, 257960, 257957, 257939, 257937, 257936, 257892, 257855, 257853, 257851, 257850, 257804, 257803, 257800, 257799, 257798, 257797, 257773, 257765, 257757, 257754, 257737, 257735, 257734, 257733, 257727, 257725, 257724, 257697, 257694, 257692, 257691, 257689, 257688, 257687, 257685, 257684, 257682, 257644, 257642, 257636, 257635, 257632, 257630, 257606, 257605, 257604, 257601, 257595, 257594, 257592, 257590, 257588, 257587, 257566, 257562, 257561, 257560, 257541, 257540, 257539, 257537, 257531, 257530, 257492, 257491, 257490, 257488, 257487, 257474, 257468, 257457, 257433, 257364, 257358, 257330, 257322, 257320, 257318, 257260, 257224, 257209, 257206, 257204, 257144, 257141, 257139, 257138, 257137, 257136, 257100, 257095, 257092, 257086, 257084, 257082, 257080, 257077, 257028, 257026, 257021, 257012, 257011, 256944, 256942, 256941, 256697, 256695, 256694, 256692, 256691, 256684, 256681, 256677, 256676, 256658, 256620, 256619, 256592, 256591, 256589, 256567, 256565, 256564, 256562, 256560, 256557, 256556, 256554, 256552, 256514, 256512, 256510, 256506, 256504};
#else
// ROOT 5 MODE
const int LHC16kRuns[] = {257605};
#endif


class AliAnalysisGrid;
void run(
         const char *taskname = "Xi1530"
         , const char *option = "LHC16k" // when scanning AOD, add "AOD"
         , const char *gridmode = "full" // or "terminate" to merge
         , UInt_t     istart = 0
         , UInt_t     iend = 25
         , const char *localorgrid = "grid"
         )
{
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libSTEERBase. so");
    gSystem->Load("libESD.so");
    gSystem->Load("libANALYSIS.so");
    gSystem->Load("libOADB.so");
    gSystem->Load("libANALYSISalice.so");
#if !defined (__CINT__) || defined (__CLING__)
    // ROOT 6 MODE
    // add aliroot indlude path
    gInterpreter->ProcessLine(Form(".include %s/include",gSystem->ExpandPathName("$ROOTSYS")));
    gInterpreter->ProcessLine(Form(".include %s/include",gSystem->ExpandPathName("$ALICE_ROOT")));
    gInterpreter->ProcessLine(Form(".include %s/include",gSystem->ExpandPathName("$ALICE_PHYSICS")));
#else
    // ROOT 5 MODE
    // add aliroot indlude path
    gROOT->ProcessLine(Form(".include %s/include",gSystem->ExpandPathName("$ROOTSYS")));
    gROOT->ProcessLine(Form(".include %s/include",gSystem->ExpandPathName("$ALICE_ROOT")));
    gROOT->ProcessLine(Form(".include %s/include",gSystem->ExpandPathName("$ALICE_PHYSICS")));
#endif
    
    bool isaa = kFALSE;
    bool ismc = kFALSE;
    bool setmixing = kFALSE;
    TString foption = option;
    if(foption.Contains("MC")) ismc = true;
    if(foption.Contains("Mix")) setmixing = true;
    
    // analysis manager
    AliAnalysisManager* mgr = new AliAnalysisManager(Form("%s%s",taskname,option));
    AliInputEventHandler* handler;
    if (foption.Contains("AOD")){
        handler = new AliAODInputHandler();
        std::cout << "AOD Hanlder!" << std::endl;
    }
    else{
        handler = new AliESDInputHandler();
        std::cout << "ESD Hanlder!" << std::endl;
    }
    
    //handler->SetNeedField(1);
    mgr->SetInputEventHandler(handler);
    
    if(ismc) {
        AliMCEventHandler *mcHandler  = new AliMCEventHandler();
        mgr->SetMCtruthEventHandler(mcHandler);
    }
#if !defined (__CINT__) || defined (__CLING__)
    // ROOT 6 MODE
    //
    // Physics Selection
    if (!foption.Contains("AOD")){
        AliPhysicsSelectionTask *physSelTask = reinterpret_cast<AliPhysicsSelectionTask*>(gInterpreter->ExecuteMacro(Form("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C(%d)",ismc)));
        if(!physSelTask) {
            Printf("no physSelTask");
            return;
        }
    }
    // Multiplicity selection
    AliMultSelectionTask *MultSlection = reinterpret_cast<AliMultSelectionTask*>(gInterpreter->ExecuteMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"));
    if(!MultSlection) {
        Printf("no MultSlection");
        return;
    }
    // PID response
    AliAnalysisTask *fPIDResponse = reinterpret_cast<AliAnalysisTaskXi1530*>(gInterpreter->ExecuteMacro(Form("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C(%d)",ismc)));
    if(!fPIDResponse) {
        Printf("no fPIDResponse");
        return;
    }
    // V0, Xi Super verexter by David
    //gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/STRANGENESS/Cascades/Run2/macros/AddTaskWeakDecayVertexer.C");
    //AliAnalysisTaskWeakDecayVertexer *taskWDV = AddTaskWeakDecayVertexer();
    //taskWDV->SetRunV0Vertexer(kTRUE);
    //taskWDV->SetRunCascadeVertexer(kFALSE);
    //taskWDV -> SetDoImprovedCascadeVertexFinding(kTRUE);
    //taskWDV -> SetDoImprovedDCAV0DauPropagation(kTRUE);
    //taskWDV -> SetDoImprovedCascadePosition(kTRUE);
    //taskWDV->SetOnlyCombineMCTrue(kTRUE);
    //taskWDV -> SetCentralityInterval(0,0.5);
    
    gInterpreter->LoadMacro("AliAnalysisTaskXi1530.cxx+g");
    
    AliAnalysisTaskXi1530 *myTask = reinterpret_cast<AliAnalysisTaskXi1530*>(gInterpreter->ExecuteMacro(Form("AddTaskXi1530.c(\"%s\",\"%s\",%d,%d,%d)",taskname,option,isaa,ismc,setmixing)));
#else
    // ROOT 5 MODE
    //
    // Physics Selection
    if (!foption.Contains("AOD")){
        gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
        AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(ismc);
        if(!physSelTask) {
            Printf("no physSelTask");
            return;
        }
    }
    
    // Multiplicity selection
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
    AliMultSelectionTask *MultSlection = AddTaskMultSelection();
    
    // PID Response
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AliAnalysisTask *fPIDResponse = AddTaskPIDResponse(ismc); //! PID response object
    
    gROOT->LoadMacro("AliAnalysisTaskXi1530.cxx+g");
    
    gROOT->LoadMacro("AddTaskXi1530.c");
    AliAnalysisTaskXi1530 *myTask = AddTaskXi1530(taskname,option,isaa,ismc,setmixing);
#endif
    
    mgr->SetDebugLevel(0);
    if (!mgr->InitAnalysis()) return;
    mgr->PrintStatus();
    
    // start analysis
    Printf("Starting Analysis....");
    
    //----LOCAL MODE-------------------------------------------------
    if(strcmp(localorgrid,"local")==0){
        TChain* chain = new TChain("ESDTree");
#if !defined (__CINT__) || defined (__CLING__)
        // ROOT 6 MODE
        std::stringstream esdChain;
        esdChain << ".x " << gSystem->Getenv("ALICE_PHYSICS") <<  "/PWG/EMCAL/macros/CreateESDChain.C(";
        esdChain << "\"" << "data.txt" << "\", ";
        esdChain << 1 << ", ";
        esdChain << 0 << ", ";
        esdChain << std::boolalpha << kFALSE << ");";
        chain = reinterpret_cast<TChain *>(gROOT->ProcessLine(esdChain.str().c_str()));
        
        chain->Lookup();
#else
        // ROOT 5 MODE
        gROOT->LoadMacro("$ALICE_PHYSICS/PWG0/CreateESDChain.C");
        TChain* chain = CreateESDChain("data.txt");
        chain->Lookup();
#endif
        mgr->StartAnalysis(localorgrid,chain);
    }
    //----GRID  MODE-------------------------------------------------
    else{
        // create the alien handler and attach it to the manager
        AliAnalysisAlien *plugin = new AliAnalysisAlien();
        plugin->SetAPIVersion("V1.1x");
        plugin->SetAliPhysicsVersion("vAN-20180925-1");
        if(!ismc)plugin->SetRunPrefix("000");
        //plugin->SetDropToShell(0);
        plugin->SetNrunsPerMaster(0);
        
        if (foption.Contains("LHC16k")){
            if(!foption.Contains("MC")){
                plugin->SetGridDataDir("/alice/data/2016/LHC16k");
                plugin->SetDataPattern("/pass2/*/AliESDs.root");
            }
            else {
                plugin->SetGridDataDir("/alice/sim/2018/LHC18c6b2");
                plugin->SetDataPattern("*AliESDs.root");
            }
#if !defined (__CINT__) || defined (__CLING__)
            for(auto i=0u;i<LHC16kRuns.size();i++)
                plugin->AddRunNumber(LHC16kRuns.at(i));
#else
            for (int i=0; i<1; i++) plugin->AddRunNumber(LHC16kRuns[i]);
#endif
        }
        plugin->SetGridWorkingDir(Form("%s%s",taskname,option));
        plugin->SetGridOutputDir("out");
        plugin->AddIncludePath("-I$ALICE_ROOT/include  -I$ALICE_ROOT/lib -I$ALICE_PHYSICS/include -I$ALICE_PHYSICS/lib -I$ALICE_PHYSICS/OADB/macros" );
        plugin->SetAnalysisSource("AliXiStarpp13TeVDevelEventCollection.cxx AliXiStarpp13TeVDevel.cxx");
        plugin->SetAnalysisSource("AliAnalysisTaskXi1530.cxx");
        plugin->SetAdditionalLibs("AliAnalysisTaskXi1530.cxx AliAnalysisTaskXi1530.h");
        plugin->SetDefaultOutputs(kFALSE);
        //plugin->SetOutputFiles("AnalysisResults.root RecTree.root");
        plugin->SetOutputFiles("AnalysisResults.root");
        
        
        plugin->SetSplitMaxInputFileNumber(2000);
        plugin->SetExecutable("myTask.sh");
        plugin->SetTTL(10000);
        plugin->SetJDLName("myTask.jdl");
        //plugin->SetOutputToRunNo(kTRUE);
        plugin->SetKeepLogs(kTRUE);
        plugin->SetMaxMergeStages(3);
        plugin->SetMaxMergeFiles(100);
        plugin->SetMergeViaJDL(kTRUE);
        plugin->SetCheckCopy(kFALSE);
        
        plugin->SetGridWorkingDir(Form("%s%s",taskname,option));
        plugin->SetGridOutputDir("out");
        
        plugin->SetOverwriteMode(0);
        plugin->SetUser("blim");
        
        mgr->SetGridHandler(plugin);
        
        // speficy on how many files you want to run
        if(strcmp(gridmode,"test")==0)plugin->SetNtestFiles(1);
        
        plugin->SetRunMode(gridmode);
        mgr->SetGridHandler(plugin);
        mgr->StartAnalysis(localorgrid);
    }
}
