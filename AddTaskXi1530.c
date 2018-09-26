#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskXi1530.h"
#include <TString.h>
#include <TList.h>
#endif

AliAnalysisTaskXi1530* AddTaskXi1530(const char *taskname = "Xi1530", const char *option = "LHC16k", bool isaa=kFALSE, bool ismc=kFALSE)
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }
    
    AliAnalysisTaskXi1530 *taskXi1530 = new AliAnalysisTaskXi1530(taskname, Form("%s_%s",taskname,option));
    //taskXi1530 -> SetFilterBit(768);
    taskXi1530 -> SetIsAA(isaa);
    taskXi1530 -> SetMixing(kFALSE);
    taskXi1530 -> SetIsMC(ismc);
    taskXi1530 -> SetParticleType(99999);
    
    if(!taskXi1530) return 0x0;
    mgr->AddTask(taskXi1530);
    
    
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutputXi1530 = mgr->CreateContainer("outputXi1530", TList::Class(), AliAnalysisManager::kOutputContainer,"AnalysisResults.root");
    
    mgr->ConnectInput(taskXi1530, 0, cinput);
    mgr->ConnectOutput(taskXi1530, 1, coutputXi1530);

    return taskXi1530;
}

