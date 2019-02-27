#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisTaskXi1530temp.h"
#endif


#if !defined (__CINT__) || defined (__CLING__)
// ROOT 6 MODE
vector<Int_t> LHC16k = {258537, 258499, 258477, 258456, 258454, 258452, 258426, 258393, 258391, 258387, 258359, 258336, 258332, 258307, 258306, 258303, 258302, 258301, 258299, 258278, 258274, 258273, 258271, 258270, 258258, 258257, 258256, 258204, 258203, 258202, 258198, 258197, 258178, 258117, 258114, 258113, 258109, 258108, 258107, 258063, 258062, 258060, 258059, 258053, 258049, 258045, 258042, 258041, 258039, 258019, 258017, 258014, 258012, 258008, 258003, 257992, 257989, 257986, 257979, 257963, 257960, 257957, 257939, 257937, 257936, 257855, 257853, 257851, 257850, 257804, 257803, 257800, 257799, 257798, 257797, 257773, 257765, 257757, 257754, 257737, 257735, 257734, 257733, 257727, 257725, 257724, 257697, 257694, 257692, 257691, 257689, 257688, 257687, 257685, 257684, 257682, 257644, 257642, 257636, 257635, 257632, 257630, 257606, 257605, 257604, 257601, 257595, 257594, 257592, 257590, 257588, 257587, 257566, 257562, 257561, 257560, 257541, 257540, 257539, 257537, 257531, 257530, 257492, 257491, 257490, 257488, 257487, 257474, 257468, 257457, 257433, 257364, 257358, 257330, 257322, 257320, 257318, 257260, 257224, 257209, 257206, 257204, 257144, 257141, 257139, 257138, 257137, 257136, 257100, 257095, 257092, 257086, 257084, 257082, 257080, 257077, 257012, 257011, 256944, 256942, 256941, 256697, 256695, 256694, 256692, 256691, 256684, 256681, 256677, 256676, 256658, 256620, 256619, 256592, 256591, 256589, 256567, 256565, 256564, 256562, 256560, 256557, 256556, 256554, 256552, 256514, 256512, 256510, 256506, 256504}; // for run
vector<Int_t> LHC16l = {259888, 259868, 259867, 259866, 259860, 259842, 259841, 259822, 259789, 259788, 259781, 259756, 259752, 259751, 259750, 259748, 259747, 259477, 259473, 259396, 259395, 259394, 259389, 259388, 259382, 259378, 259342, 259341, 259340, 259339, 259336, 259334, 259307, 259305, 259303, 259302, 259274, 259273, 259272, 259271, 259270, 259269, 259264, 259263, 259261, 259257, 259204, 259164, 259162, 259118, 259117, 259099, 259096, 259091, 259090, 259088, 258964, 258962}; // for run
vector<Int_t> LHC16o = {264035, 264033, 263985, 263984, 263981, 263978, 263977, 263923, 263920, 263917, 263916, 263905, 263866, 263863, 263810, 263803, 263793, 263792, 263790, 263787, 263786, 263785, 263784, 263744, 263743, 263741, 263739, 263738, 263737, 263691, 263690, 263682, 263663, 263662, 263657, 263654, 263652, 263647, 263529, 263497, 263496, 263490, 263487, 263332, 263331, 262858, 262855, 262853, 262849, 262847, 262844, 262842, 262841, 262778, 262777, 262776, 262768, 262760, 262727, 262725, 262723, 262719, 262717, 262713, 262708, 262706, 262705, 262428, 262426, 262425, 262424};
vector<Int_t> LHC16p = {264347, 264346, 264345, 264341, 264336, 264312, 264306, 264305, 264281, 264279, 264277, 264273, 264267, 264266, 264265, 264264, 264262, 264261, 264260, 264259, 264238, 264235, 264233, 264232, 264198, 264197, 264194, 264190, 264188, 264168, 264164, 264139, 264138, 264137, 264129, 264110, 264109, 264086, 264085, 264082, 264078, 264076};
vector<Int_t> LHC17k = {276508, 276507, 276506, 276462, 276439, 276438, 276437, 276435, 276429, 276351, 276348, 276302, 276297, 276294, 276292, 276290, 276259, 276257, 276230, 276205, 276178, 276177, 276170, 276169, 276166, 276145, 276140, 276135, 276104, 276102, 276099, 276098, 276097, 275847, 275664, 275661, 275650, 275648, 275647, 275624, 275623, 275622, 275621, 275617, 275612, 275559, 275558, 275515, 275472, 275471, 275467, 275459, 275457, 275456, 275453, 275452, 275448, 275443, 275406, 275404, 275401, 275372, 275369, 275361, 275360, 275333, 275332, 275328, 275326, 275324, 275322, 275314, 275283, 275247, 275246, 275245, 275239, 275188, 275184, 275180, 275177, 275174, 275173, 275151, 275150, 275149, 275076, 275075, 275073, 275068, 275067, 274979, 274978, 274886, 274882, 274878, 274877, 274822, 274821, 274815, 274806, 274803, 274802, 274801, 274708, 274690};
vector<Int_t> LHC17l = {278216, 278215, 278191, 278189, 278166, 278165, 278164, 278158, 278127, 278126, 278123, 278122, 277996, 277991, 277989, 277987, 277952, 277930, 277907, 277904, 277903, 277900, 277899, 277898, 277897, 277876, 277870, 277848, 277847, 277845, 277842, 277841, 277834, 277805, 277802, 277801, 277800, 277799, 277795, 277794, 277749, 277747, 277746, 277745, 277725, 277723, 277722, 277721, 277577, 277576, 277575, 277574, 277537, 277536, 277534, 277531, 277530, 277479, 277478, 277477, 277476, 277473, 277472, 277418, 277417, 277416, 277389, 277386, 277385, 277384, 277383, 277360, 277314, 277312, 277310, 277293, 277262, 277257, 277256, 277197, 277196, 277194, 277193, 277189, 277188, 277184, 277183, 277182, 277181, 277180, 277155, 277121, 277117, 277091, 277087, 277082, 277079, 277076, 277073, 277037, 277017, 277016, 277015, 276972, 276971, 276970, 276969, 276967, 276920, 276917, 276916, 276762, 276675, 276674, 276672, 276671, 276670, 276644, 276608, 276557, 276556, 276553, 276552, 276551};
vector<Int_t> LHC17m = {280140, 280135, 280134, 280131, 280126, 280118, 280114, 280111, 280108, 280107, 280066, 280052, 280051, 279879, 279855, 279854, 279853, 279830, 279827, 279826, 279773, 279749, 279747, 279719, 279718, 279715, 279689, 279688, 279687, 279684, 279683, 279682, 279679, 279677, 279676, 279642, 279641, 279632, 279630, 279559, 279550, 279491, 279488, 279487, 279483, 279441, 279439, 279435, 279410, 279391, 279355, 279354, 279349, 279348, 279344, 279342, 279312, 279310, 279309, 279274, 279273, 279270, 279268, 279267, 279265, 279264, 279242, 279238, 279235, 279234, 279232, 279208, 279207, 279201, 279199, 279157, 279155, 279130, 279123, 279122, 279118, 279117, 279107, 279106, 279075, 279074, 279073, 279069, 279068, 279044, 279043, 279041, 279036, 279035, 279008, 279007, 279005, 279000, 278999, 278964, 278963, 278960, 278959, 278941, 278939, 278936, 278915, 278914};
vector<Int_t> LHC17o = {281961, 281956, 281953, 281940, 281939, 281932, 281931, 281928, 281920, 281918, 281916, 281915, 281895, 281894, 281893, 281892, 281633, 281592, 281583, 281574, 281569, 281568, 281563, 281562, 281557, 281511, 281509, 281477, 281475, 281450, 281449, 281446, 281444, 281443, 281441, 281415, 281321, 281301, 281277, 281275, 281273, 281271, 281244, 281243, 281242, 281241, 281240, 281213, 281212, 281191, 281190, 281189, 281181, 281180, 281179, 281081, 281080, 281062, 281061, 281060, 281036, 281035, 281033, 281032, 280999, 280998, 280997, 280996, 280994, 280990, 280947, 280943, 280940, 280936, 280897, 280880, 280856, 280854, 280849, 280848, 280847, 280844, 280842, 280793, 280792, 280787, 280786, 280768, 280767, 280766, 280765, 280764, 280763, 280762, 280761, 280757, 280756, 280755, 280754, 280753, 280729, 280706, 280705, 280681, 280679, 280671, 280647, 280645, 280639, 280637, 280636, 280634, 280613, 280583, 280581, 280574, 280551, 280550, 280547, 280546, 280519, 280518, 280499, 280490, 280448, 280447, 280446, 280445, 280443, 280419, 280415, 280412, 280406, 280405, 280403, 280375, 280374, 280351, 280350, 280349, 280348, 280312, 280310, 280290, 280286, 280285, 280284, 280282};
vector<Int_t> LHC17r = {282704, 282703, 282702, 282700, 282677, 282676, 282673, 282671, 282670, 282667, 282666, 282651, 282629, 282622, 282620, 282618, 282609, 282608, 282607, 282606, 282580, 282579, 282575, 282573, 282546, 282545, 282544, 282528};
vector<Int_t> LHC10b = {114931, 115186, 115193, 115393, 115401, 116102, 116288, 116402, 116403, 116562, 116571, 116574, 116643, 116645, 117048, 117050, 117052, 117053, 117059, 117060, 117063, 117092, 117099, 117109, 117112, 117116, 117220, 117222, 117054,117065, 117077, 117086};
#else
// ROOT 5 MODE
const int LHC16k[] = {258537, 258499, 258477, 258456, 258454, 258452, 258426, 258393, 258391, 258387, 258359, 258336, 258332, 258307, 258306, 258303, 258302, 258301, 258299, 258278, 258274, 258273, 258271, 258270, 258258, 258257, 258256, 258204, 258203, 258202, 258198, 258197, 258178, 258117, 258114, 258113, 258109, 258108, 258107, 258063, 258062, 258060, 258059, 258053, 258049, 258045, 258042, 258041, 258039, 258019, 258017, 258014, 258012, 258008, 258003, 257992, 257989, 257986, 257979, 257963, 257960, 257957, 257939, 257937, 257936, 257855, 257853, 257851, 257850, 257804, 257803, 257800, 257799, 257798, 257797, 257773, 257765, 257757, 257754, 257737, 257735, 257734, 257733, 257727, 257725, 257724, 257697, 257694, 257692, 257691, 257689, 257688, 257687, 257685, 257684, 257682, 257644, 257642, 257636, 257635, 257632, 257630, 257606, 257605, 257604, 257601, 257595, 257594, 257592, 257590, 257588, 257587, 257566, 257562, 257561, 257560, 257541, 257540, 257539, 257537, 257531, 257530, 257492, 257491, 257490, 257488, 257487, 257474, 257468, 257457, 257433, 257364, 257358, 257330, 257322, 257320, 257318, 257260, 257224, 257209, 257206, 257204, 257144, 257141, 257139, 257138, 257137, 257136, 257100, 257095, 257092, 257086, 257084, 257082, 257080, 257077, 257012, 257011, 256944, 256942, 256941, 256697, 256695, 256694, 256692, 256691, 256684, 256681, 256677, 256676, 256658, 256620, 256619, 256592, 256591, 256589, 256567, 256565, 256564, 256562, 256560, 256557, 256556, 256554, 256552, 256514, 256512, 256510, 256506, 256504};
const int LHC16l[] = {259888, 259868, 259867, 259866, 259860, 259842, 259841, 259822, 259789, 259788, 259781, 259756, 259752, 259751, 259750, 259748, 259747, 259477, 259473, 259396, 259395, 259394, 259389, 259388, 259382, 259378, 259342, 259341, 259340, 259339, 259336, 259334, 259307, 259305, 259303, 259302, 259274, 259273, 259272, 259271, 259270, 259269, 259264, 259263, 259261, 259257, 259204, 259164, 259162, 259118, 259117, 259099, 259096, 259091, 259090, 259088, 258964, 258962};
#endif


class AliAnalysisGrid;
void run(
         const char *taskname = "Xi1530"
         , const char *option = "LHC16k_pass2_SYS_MC_Mix_test" // when scanning AOD, add "AOD"
         , const char *gridmode = "test" // or "terminate" to merge
         , UInt_t     istart = 0
         , UInt_t     iend = 25
         , const char *localorgrid = "local"
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
    gSystem->Load("libpythia6_4_21.so");
    
    gSystem->Setenv("alien_CLOSE_SE","working_disk_SE");
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
    int nmix = 20;
    bool highmult = kFALSE;
    TString foption = option;
    const char* suffix = "test";
    if(foption.Contains("MC"))ismc = kTRUE;
    
    // analysis manager
    AliAnalysisManager* mgr = new AliAnalysisManager(Form("%s%s",taskname,option));
    AliInputEventHandler* handler;
    if (foption.Contains("AOD")) handler = new AliAODInputHandler();
    else handler = new AliESDInputHandler();
    
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
    AliPhysicsSelectionTask *physSelTask = reinterpret_cast<AliPhysicsSelectionTask*>(gInterpreter->ExecuteMacro(Form("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C(%d)",ismc)));
    if(!physSelTask) {
        Printf("no physSelTask");
        return;
    }
    // Multiplicity selection
    AliMultSelectionTask *MultSlection = reinterpret_cast<AliMultSelectionTask*>(gInterpreter->ExecuteMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"));
    if(!MultSlection) {
        Printf("no MultSlection");
        return;
    }
    MultSlection->SetAddInfo(kTRUE);
    MultSlection->SetSelectedTriggerClass(AliVEvent::kAny);
    //MultSlection->SetAlternateOADBforEstimators("LHC16k"); //if needed
    // PID response
    AliAnalysisTask *fPIDResponse = reinterpret_cast<AliAnalysisTaskXi1530temp*>(gInterpreter->ExecuteMacro(Form("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C(%d)",ismc)));
    if(!fPIDResponse) {
        Printf("no fPIDResponse");
        return;
    }
    // V0, Xi Super verexter by David
    /*
    AliAnalysisTaskWeakDecayVertexer *taskWDV = reinterpret_cast<AliAnalysisTaskWeakDecayVertexer *>(gInterpreter->ExecuteMacro("$ALICE_PHYSICS/PWGLF/STRANGENESS/Cascades/Run2/macros/AddTaskWeakDecayVertexer.C"));
    taskWDV->SetUseImprovedFinding();

    //______________________________________________________________
    //Revertexing configuration
    //WARNING: This applies only to the Cascade analysis
    taskWDV->SetUseImprovedFinding();

    //V0-Related topological selections
    taskWDV->SetV0VertexerDCAFirstToPV(0.05);
    taskWDV->SetV0VertexerDCASecondtoPV(0.05);
    taskWDV->SetV0VertexerDCAV0Daughters(1.6);
    taskWDV->SetV0VertexerCosinePA(0.97);
    taskWDV->SetV0VertexerMinRadius(0.5);
    taskWDV->SetV0VertexerMaxRadius(200);

    //Cascade-Related topological selections
    taskWDV->SetCascVertexerMinV0ImpactParameter(0.05);
    taskWDV->SetCascVertexerV0MassWindow(0.007);
    taskWDV->SetCascVertexerDCABachToPV(0.05);
    taskWDV->SetCascVertexerDCACascadeDaughters(1.6);
    taskWDV->SetCascVertexerCascadeMinRadius(.5);
    taskWDV->SetCascVertexerCascadeCosinePA(.97);
    */
    /*  
    AliAnalysisTaskWeakDecayVertexer *taskWDV = reinterpret_cast<AliAnalysisTaskWeakDecayVertexer*>(gInterpreter->ExecuteMacro("$ALICE_PHYSICS/PWGLF/STRANGENESS/Cascades/Run2/macros/AddTaskWeakDecayVertexer.C"));
    //______________________________________________________________
    //Revertexing configuration
    //WARNING: This applies only to the Cascade analysis
    taskWDV ->SetUseImprovedFinding(); 

        //V0-Related topological selections
    taskWDV ->     SetV0VertexerDCAFirstToPV(0.05);
    taskWDV ->     SetV0VertexerDCASecondtoPV(0.05);
    taskWDV ->     SetV0VertexerDCAV0Daughters(1.6);
    taskWDV ->     SetV0VertexerCosinePA(0.97);
    taskWDV ->     SetV0VertexerMinRadius(0.5);
    taskWDV ->     SetV0VertexerMaxRadius(200);
        
        //Cascade-Related topological selections
    taskWDV ->     SetCascVertexerMinV0ImpactParameter(0.05);
    taskWDV ->     SetCascVertexerV0MassWindow(0.007);
    taskWDV ->     SetCascVertexerDCABachToPV(0.05);
    taskWDV ->     SetCascVertexerDCACascadeDaughters(1.6);
    taskWDV ->     SetCascVertexerCascadeMinRadius(.5);
    taskWDV ->     SetCascVertexerCascadeCosinePA(.97);
*/
    gInterpreter->LoadMacro("AliAnalysisTaskXi1530temp.cxx+g");
    //AliAnalysisTaskXi1530temp *myTask = reinterpret_cast<AliAnalysisTaskXi1530temp*>(gInterpreter->ExecuteMacro(Form("AddTaskXi1530.c(\"%s\",\"%s\",%i,%d,%d,%d,%d)",taskname,option,nmix,highmult,isaa,ismc,setmixing)));
    
    AliAnalysisTaskXi1530temp *myTask = reinterpret_cast<AliAnalysisTaskXi1530temp*>(gInterpreter->ExecuteMacro(Form("AddTaskXi1530.C(\"%s\",\"%s\",%i,\"%s\")",taskname,option,nmix,suffix)));
#else
    // ROOT 5 MODE
    //
    // Physics Selection
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(ismc);
    if(!physSelTask) {
        Printf("no physSelTask");
        return;
    }
    
    // Multiplicity selection
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
    AliMultSelectionTask *MultSlection = AddTaskMultSelection();

    // PID Response
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AliAnalysisTask *fPIDResponse = AddTaskPIDResponse(ismc); //! PID response object
    
    gROOT->LoadMacro("AliAnalysisTaskXi1530temp.cxx+g");
    
    gROOT->LoadMacro("AddTaskXi1530.C");
    AliAnalysisTaskXi1530temp *myTask = AddTaskXi1530(taskname,option,nmix,highmult,isaa,ismc,setmixing);
#endif
    
    //mgr->SetDebugLevel(1);
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
        if(!ismc) esdChain << "\"" << "data.txt" << "\", ";
        else esdChain << "\"" << "data_MC.txt" << "\", ";
        esdChain << 1 << ", ";
        esdChain << 0 << ", ";
        esdChain << std::boolalpha << kFALSE << ");";
        chain = reinterpret_cast<TChain *>(gROOT->ProcessLine(esdChain.str().c_str()));
        
        chain->Lookup();
#else
        // ROOT 5 MODE
        gROOT->LoadMacro("$ALICE_PHYSICS/PWG0/CreateESDChain.C");
        if(!ismc) TChain* chain = CreateESDChain("data.txt");
        else TChain* chain = CreateESDChain("data_MC.txt");
        chain->Lookup();
#endif
        mgr->StartAnalysis(localorgrid,chain);
    }
    //----GRID  MODE-------------------------------------------------
    else{
        // create the alien handler and attach it to the manager
        AliAnalysisAlien *plugin = new AliAnalysisAlien();
        plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
        //plugin->Load("libpythia6_4_21.so");
        plugin->SetAnalysisSource("AliAnalysisTaskXi1530temp.cxx");
        plugin->SetAdditionalLibs("AliAnalysisTaskXi1530temp.cxx AliAnalysisTaskXi1530temp.h libpythia6_4_21.so");
        plugin->SetAliPhysicsVersion("vAN-20190212_ROOT6-1");
        plugin->SetAPIVersion("V1.1x");
        if(!ismc)plugin->SetRunPrefix("000");
        //plugin->SetDropToShell(0);
        
        if (foption.Contains("LHC16k")){
            if(ismc){
                if (foption.Contains("Gen")){
                    plugin->SetGridDataDir("/alice/sim/2018/LHC18f1"); //general use, LHC18f1 for pass2, LHC17d20a1 for pass1
                }
                else{
                    if (foption.Contains("pass2")) plugin->SetGridDataDir("/alice/sim/2018/LHC18c6b4"); //resonance injected
                    else plugin->SetGridDataDir("/alice/sim/2018/LHC18c6b"); //resonance injected            
                }
                plugin->SetDataPattern("*AliESDs.root");
            }
            else {
                plugin->SetGridDataDir("/alice/data/2016/LHC16k");
                if (foption.Contains("pass2")) plugin->SetDataPattern("pass2/*/AliESDs.root");
                else plugin->SetDataPattern("pass1/*/AliESDs.root");
            }
            Int_t end = LHC16k.size();
            if (foption.Contains("test")) end = 1;
            for(auto i=0;i<end;i++)
                plugin->AddRunNumber(LHC16k.at(i));
        }
        if (foption.Contains("LHC16l")){
            if(ismc){
                if (foption.Contains("Gen")){
                    plugin->SetGridDataDir("/alice/sim/2018/LHC18d8"); //general use
                }
                else{
                    if (foption.Contains("pass2")) plugin->SetGridDataDir("/alice/sim/2018/LHC18c6b4"); //resonance injected
                else plugin->SetGridDataDir("/alice/sim/2018/LHC18c6b"); //resonance injected          
                }
                plugin->SetDataPattern("*AliESDs.root");
            }
            else {
                plugin->SetGridDataDir("/alice/data/2016/LHC16l");
                if (foption.Contains("pass2")) plugin->SetDataPattern("pass2/*/AliESDs.root");
                else plugin->SetDataPattern("pass1/*/AliESDs.root");
            }
            Int_t end = LHC16l.size();
            if (foption.Contains("test")) end = 1;
            for(auto i=0u;i<end;i++)
                plugin->AddRunNumber(LHC16l.at(i));
        }
        if (foption.Contains("LHC16o")){
            if(ismc){
                plugin->SetGridDataDir("/alice/sim/2018/LHC18c6b"); //resonance injected
                //plugin->SetGridDataDir("/alice/sim/2017/LHC17d20a1"); //general use
                plugin->SetDataPattern("*AliESDs.root");
            }
            else {
                plugin->SetGridDataDir("/alice/data/2016/LHC16o");
                plugin->SetDataPattern("pass1/*/AliESDs.root");
            }
            Int_t end = LHC16o.size();
            if (foption.Contains("test")) end = 1;
            for(auto i=0u;i<end;i++)
                plugin->AddRunNumber(LHC16o.at(i));
        }
        if (foption.Contains("LHC16p")){
            if(ismc){
                plugin->SetGridDataDir("/alice/sim/2018/LHC18c6b"); //resonance injected
                //plugin->SetGridDataDir("/alice/sim/2017/LHC17d20a1"); //general use
                plugin->SetDataPattern("*AliESDs.root");
            }
            else {
                plugin->SetGridDataDir("/alice/data/2016/LHC16p");
                plugin->SetDataPattern("pass1/*/AliESDs.root");
            }
            Int_t end = LHC16p.size();
            if (foption.Contains("test")) end = 1;
            for(auto i=0u;i<end;i++)
                plugin->AddRunNumber(LHC16p.at(i));
        }
        if (foption.Contains("LHC17k")){
            if(ismc){
                plugin->SetGridDataDir("/alice/sim/2018/LHC18c6a"); //resonance injected
                //plugin->SetGridDataDir("/alice/sim/2017/LHC17d20a1"); //general use
                plugin->SetDataPattern("*AliESDs.root");
            }
            else {
                plugin->SetGridDataDir("/alice/data/2017/LHC17k");
                plugin->SetDataPattern("pass1/*/AliESDs.root");
            }
            Int_t end = LHC17k.size();
            if (foption.Contains("test")) end = 1;
            for(auto i=0u;i<end;i++)
                plugin->AddRunNumber(LHC17k.at(i));
        }
        
        if (foption.Contains("LHC17l")){
            if(ismc){
                plugin->SetGridDataDir("/alice/sim/2018/LHC18c6a"); //resonance injected
                //plugin->SetGridDataDir("/alice/sim/2017/LHC17d20a1"); //general use
                plugin->SetDataPattern("*AliESDs.root");
            }
            else {
                plugin->SetGridDataDir("/alice/data/2017/LHC17l");
                plugin->SetDataPattern("pass1/*/AliESDs.root");
            }
            Int_t end = LHC17l.size();
            if (foption.Contains("test")) end = 1;
            for(auto i=0u;i<end;i++)
                plugin->AddRunNumber(LHC17l.at(i));
        }
        if (foption.Contains("LHC17m")){
            if(ismc){
                plugin->SetGridDataDir("/alice/sim/2018/LHC18c6a"); //resonance injected
                //plugin->SetGridDataDir("/alice/sim/2017/LHC17d20a1"); //general use
                plugin->SetDataPattern("*AliESDs.root");
            }
            else {
                plugin->SetGridDataDir("/alice/data/2017/LHC17m");
                plugin->SetDataPattern("pass1/*/AliESDs.root");
            }
            Int_t end = LHC17m.size();
            if (foption.Contains("test")) end = 1;
            for(auto i=0u;i<end;i++)
                plugin->AddRunNumber(LHC17m.at(i));
        }
        if (foption.Contains("LHC17o")){
            if(ismc){
                plugin->SetGridDataDir("/alice/sim/2018/LHC18c6a"); //resonance injected
                //plugin->SetGridDataDir("/alice/sim/2017/LHC17d20a1"); //general use
                plugin->SetDataPattern("*AliESDs.root");
            }
            else {
                plugin->SetGridDataDir("/alice/data/2017/LHC17o");
                plugin->SetDataPattern("pass1/*/AliESDs.root");
            }
            Int_t end = LHC17o.size();
            if (foption.Contains("test")) end = 1;
            for(auto i=0u;i<end;i++)
                plugin->AddRunNumber(LHC17o.at(i));
        }
        if (foption.Contains("LHC17r")){
            if(ismc){
                plugin->SetGridDataDir("/alice/sim/2018/LHC18c6a"); //resonance injected
                //plugin->SetGridDataDir("/alice/sim/2017/LHC17d20a1"); //general use
                plugin->SetDataPattern("*AliESDs.root");
            }
            else {
                plugin->SetGridDataDir("/alice/data/2017/LHC17r");
                plugin->SetDataPattern("pass1/*/AliESDs.root");
            }
            Int_t end = LHC17r.size();
            if (foption.Contains("test")) end = 1;
            for(auto i=0u;i<end;i++)
                plugin->AddRunNumber(LHC17r.at(i));
        }
        if (foption.Contains("LHC10b")){
            if(ismc){
                plugin->SetGridDataDir("/alice/sim/LHC10d1");
                plugin->SetDataPattern("*AliESDs.root");
            }
            else {
                plugin->SetGridDataDir("/alice/data/2010/LHC10b");
                plugin->SetDataPattern("pass4/*/AliESDs.root");
            }
            Int_t end = LHC10b.size();
            if (foption.Contains("test")) end = 1;
            for(auto i=0u;i<end;i++)
                plugin->AddRunNumber(LHC10b.at(i));
        }
        plugin->SetSplitMaxInputFileNumber(8000);
        plugin->SetExecutable(Form("%s%s.sh",taskname,option));
        plugin->SetTTL(40000);
        plugin->SetJDLName(Form("%s%s.jdl",taskname,option));
        plugin->SetKeepLogs(kTRUE);
        //plugin->SetMaxMergeStages(3);
        plugin->SetMaxMergeFiles(100);
        plugin->SetMergeViaJDL(kTRUE);
        plugin->SetOneStageMerging(kFALSE);
        plugin->SetCheckCopy(kFALSE);
        //plugin->SetNrunsPerMaster(kFALSE);
        plugin->SetUseSubmitPolicy(kTRUE);
        
        plugin->SetGridWorkingDir(Form("%s%s",taskname,option));
        plugin->SetGridOutputDir("out");
        
        //plugin->SetOutputToRunNo(kTRUE);
        plugin->SetOverwriteMode(kTRUE);
        plugin->SetUser("blim");
        mgr->SetGridHandler(plugin);
        
        plugin->SetRunMode(gridmode);
        if(strcmp(gridmode,"test")==0)plugin->SetNtestFiles(1);
        
        mgr->StartAnalysis(localorgrid);
    }
}
