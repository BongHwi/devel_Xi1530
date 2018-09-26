#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisTaskXi1530.h"
#endif

vector<Int_t> LHC16kRuns = {257605}; // for test
//const int LHC16kRuns[]={258537, 258499, 258477, 258456, 258454, 258452, 258426, 258393, 258391, 258387, 258359, 258336, 258332, 258307, 258306, 258303, 258302, 258301, 258299, 258278, 258274, 258273, 258271, 258270, 258258, 258257, 258256, 258204, 258203, 258202, 258198, 258197, 258178, 258117, 258114, 258113, 258109, 258108, 258107, 258063, 258062, 258060, 258059, 258053, 258049, 258045, 258042, 258041, 258039, 258019, 258017, 258014, 258012, 258008, 258003, 257992, 257989, 257986, 257979, 257963, 257960, 257957, 257939, 257937, 257936, 257892, 257855, 257853, 257851, 257850, 257804, 257803, 257800, 257799, 257798, 257797, 257773, 257765, 257757, 257754, 257737, 257735, 257734, 257733, 257727, 257725, 257724, 257697, 257694, 257692, 257691, 257689, 257688, 257687, 257685, 257684, 257682, 257644, 257642, 257636, 257635, 257632, 257630, 257606, 257605, 257604, 257601, 257595, 257594, 257592, 257590, 257588, 257587, 257566, 257562, 257561, 257560, 257541, 257540, 257539, 257537, 257531, 257530, 257492, 257491, 257490, 257488, 257487, 257474, 257468, 257457, 257433, 257364, 257358, 257330, 257322, 257320, 257318, 257260, 257224, 257209, 257206, 257204, 257144, 257141, 257139, 257138, 257137, 257136, 257100, 257095, 257092, 257086, 257084, 257082, 257080, 257077, 257028, 257026, 257021, 257012, 257011, 256944, 256942, 256941, 256697, 256695, 256694, 256692, 256691, 256684, 256681, 256677, 256676, 256658, 256620, 256619, 256592, 256591, 256589, 256567, 256565, 256564, 256562, 256560, 256557, 256556, 256554, 256552, 256514, 256512, 256510, 256506, 256504};

class AliAnalysisGrid;
void run(
         const char *taskname = "Xi1530"
         , const char *option = "LHC16k" // when scanning AOD, add "AOD"
         , const char *gridmode = "full" // or "terminate" to merge
         , UInt_t     istart = 0
         , UInt_t     iend = 25
         , const char *localorgrid = "local"
         )
{
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
    TString foption = option;
    if(foption.Contains("MC")) ismc = true;
    
    // analysis manager
    AliAnalysisManager* mgr = new AliAnalysisManager(Form("%s%s",taskname,option));
    
    // create the alien handler and attach it to the manager
    AliAnalysisAlien *plugin = new AliAnalysisAlien();
    plugin->SetRunMode(gridmode);
    plugin->SetAPIVersion("V1.1x");
    plugin->SetAliPhysicsVersion("vAN-20180925-1");
    plugin->SetDropToShell(0);
    if(!ismc)plugin->SetRunPrefix("000");
    plugin->SetNrunsPerMaster(1);
    plugin->SetOutputToRunNo();
    plugin->SetMergeViaJDL(1);
    
    if (foption.Contains("LHC16k")){
        plugin->SetGridDataDir("/alice/data/2016/LHC16k/");
        for(auto i=0u;i<LHC16kRuns.size();i++)
            plugin->AddRunNumber(LHC16kRuns.at(i));
        plugin->SetDataPattern("/pass2/*/AliESDs.root");
    }
    
    
    plugin->SetGridWorkingDir(Form("%s%s",taskname,option));
    plugin->SetGridOutputDir("out");
    plugin->AddIncludePath("-I$ALICE_ROOT/include  -I$ALICE_ROOT/lib -I$ALICE_PHYSICS/include -I$ALICE_PHYSICS/lib -I$ALICE_PHYSICS/OADB/macros" );
    plugin->SetAnalysisSource("AliAnalysisTaskXi1530.cxx");
    plugin->SetAdditionalLibs("AliAnalysisTaskXi1530.cxx AliAnalysisTaskXi1530.h");
    plugin->SetDefaultOutputs(kFALSE);
    //plugin->SetOutputFiles("AnalysisResults.root RecTree.root");
    plugin->SetOutputFiles("AnalysisResults.root");
    plugin->SetSplitMaxInputFileNumber(300);
    plugin->SetMasterResubmitThreshold(90);
    //plugin->SetFileForTestMode("data.txt");
    //plugin->SetUseSubmitPolicy();
    
    // Optionally set time to live (default 30000 sec)
    plugin->SetTTL(20000);
    // Optionally set input format (default xml-single)
    plugin->SetInputFormat("xml-single");
    // Optionally modify the name of the generated JDL (default analysis.jdl)
    plugin->SetJDLName(Form("%s%s.jdl",taskname,option));
    // Optionally modify the executable name (default analysis.sh)
    plugin->SetExecutable(Form("%s%s.sh",taskname,option));
    // Optionally modify job price (default 1)
    plugin->SetPrice(1);
    // Optionally modify split mode (default 'se')
    plugin->SetSplitMode("se");
    
    
    mgr->SetGridHandler(plugin);
    AliInputEventHandler* handler;
    if (foption.Contains("AOD"))
        handler = new AliAODInputHandler();
    else
        handler = new AliESDInputHandler() ;
    
    handler->SetNeedField(1);
    mgr->SetInputEventHandler(handler);
    
    TChain* chain = 0;
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
    
    //!! Need to be added in ROOT6 mode!!
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
    
    //!! Need to be added in ROOT6 mode!!
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG0/CreateESDChain.C");
    TChain* chain = CreateESDChain("data.txt");
    chain->Lookup();
#endif
    //hybrid track : AOD 086 -> Filter bit 272
    //hybrid track : AOD 160 -> Filter bit 768
    //hybrid track : AOD 145 -> Filter bit 768
    //hybrid track : AOD 115 -> Filter bit 768
    // TPC AOD086 : 128
    std::cout << "Task Prepare" << std::endl;
    AliAnalysisTaskXi1530 *taskXi1530 = new AliAnalysisTaskXi1530(taskname, Form("%s_%s",taskname,option));
    //taskXi1530 -> SetFilterBit(768);
    taskXi1530 -> SetIsAA(isaa);
    taskXi1530 -> SetMixing(kFALSE);
    taskXi1530 -> SetIsMC(ismc);
    taskXi1530 -> SetParticleType(99999);
    std::cout << "After Task Ready" << std::endl;
    
    // Create containers for input/output
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutputXi1530 = mgr->CreateContainer("outputXi1530", TDirectory::Class(), AliAnalysisManager::kOutputContainer,"AnalysisResults.root");
    std::cout << "Container ready" << std::endl;
    
    mgr->AddTask(taskXi1530);
    std::cout << "AddTask" << std::endl;
    mgr->ConnectInput(taskXi1530, 0, cinput);
    std::cout << "Add Input" << std::endl;
    mgr->ConnectOutput(taskXi1530, 1, coutputXi1530);
    std::cout << "Add Output" << std::endl;
    
    std::cout << "After input/output connecting" << std::endl;
    // enable debug printouts
    mgr->SetDebugLevel(5);
    if (!mgr->InitAnalysis()) return;
    mgr->PrintStatus();
    
    // start analysis
    Printf("Starting Analysis....");
    mgr->StartAnalysis(localorgrid,chain);
}

