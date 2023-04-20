#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
//#include <iomanip>
//#include <vector>

#include "TPie.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLegend.h"

#include "../../BaselineFramework/Tools/Yields.C"
#include "../../BaselineFramework/Tools/Channel.C"
#include "../../BaselineFramework/Tools/EventCut.C"
#include "../../BaselineFramework/Tools/VariableDistr2D.C"
#include "../../BaselineFramework/Tools/VariableDistr.C"
#include "../../BaselineFramework/Tools/PhysicsSample.C"
#include "../../BaselineFramework/Tools/PhysicsProcess.C"

#include "../../BaselineFramework/Tools/MiniTreeAnalyzer.C"

//WW_ssEW4.tex  WWllvv.tex     WWqqlv.tex  WZlllv.tex  WZqqll.tex  ZZllll.tex     ZZqqll.tex
//WW_ssEW6.tex  WWllvv_gg.tex  WZ_EW6.tex  WZlvqq.tex  ZZ_EW6.tex  ZZllll_gg.tex



void Skimming(){
    
    

   // Channel *ch_jet1= new Channel("SRjet1","SRjet1","(@jets.size()==1)");
   
    MiniTreeAnalyzer analyzer;
    analyzer.bkgDir="/eos/user/d/dparedes/mc16d/ge0j/";
    analyzer.SetTreeName("nominal_Loose");

    analyzer.AddProcess("Zjets.list","Zjets","Zjets",14,1,"isBkg");

    analyzer.printLog=true;
    
    TString precut="nElectrons==2 && nMuons==0 && el_pt[0]>28000 && el_pt[1]>28000 && fabs(el_eta[0])<2.47 && fabs(el_eta[1])<2.47 && el_LHTight[0]>0 && el_LHTight[1]>0 && el_isoFCTight[0]>0 && el_isoFCTight[1]>0 && fabs(el_d0sig[0])<5 && fabs(el_d0sig[1])<5 && fabs(el_delta_z0_sintheta[0])<0.5 && fabs(el_delta_z0_sintheta[1])<0.5 && (fabs(el_eta[0])<1.37  || fabs(el_eta[0])>1.52) && (fabs(el_eta[1])<1.37  || fabs(el_eta[1])>1.52)";

    analyzer.MakeSkimming(precut);             


return;
}
