#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

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


void Skimming(){
    
    TString precut="(loose_SSee || OSee) && el_isTight[0] && el_isTight[1] && nMuons==0 && el_pt[0]>28000 && el_pt[1]>28000 && fabs(el_eta[0])<2.47 && fabs(el_eta[1])<2.47 && (fabs(el_eta[0])<1.37  || fabs(el_eta[0])>1.52) && (fabs(el_eta[1])<1.37  || fabs(el_eta[1])>1.52)"; 
    TString sample="Zjets";

    MiniTreeAnalyzer analyzer;

    // Main Settings
    analyzer.bkgDir="/eos/user/d/dparedes/Samples212560/LooseSSML/mc16e/ge0j/"; //"/eos/user/d/dparedes/mc16e/ge0j/";
    //analyzer.bkgDir="/eos/user/d/dparedes/LooseSSML/data/ge0j/"; //"/eos/user/d/dparedes/mc16e/ge0j/"; //Skimming/";
  
    analyzer.SetTreeName("nominal_Loose");
    analyzer.AddProcess(sample+".list",sample,sample,14,1,"isBkg");
    analyzer.printLog=true;

    // Variables to keep
    std::vector<TString> vars;
    vars.push_back("loose_SSee");
    vars.push_back("OSee");
    vars.push_back("el_isTight");
    vars.push_back("el_pt");
    vars.push_back("el_eta");
    vars.push_back("el_e");
    vars.push_back("el_phi");
    vars.push_back("nMuons");
    vars.push_back("mcChannelNumber");
    vars.push_back("el_ECIDS");
    vars.push_back("nElectrons");
    vars.push_back("el_true_type");
    vars.push_back("el_true_origin");
    vars.push_back("el_true_firstEgMotherTruthType"); 
    vars.push_back("el_true_firstEgMotherTruthOrigin");
    vars.push_back("el_true_firstEgMotherPdgId");
    vars.push_back("el_true_pdg");
    vars.push_back("el_charge");
    vars.push_back("weight_bTagSF_MV2c10_77");
    vars.push_back("weight_mc");
    vars.push_back("weight_leptonSF");
    vars.push_back("weight_pileup");
    vars.push_back("weight_jvt");
    vars.push_back("weight_normalise");
    vars.push_back("weight_indiv_SF_EL_ChargeID");
    vars.push_back("weight_indiv_SF_EL_ChargeMisID");
    vars.push_back("runNumber");
    vars.push_back("HT_all");
    vars.push_back("HT_jets");
    vars.push_back("nJets");
    vars.push_back("nBTags_MV2c10_77");
    vars.push_back("met_met");
    vars.push_back("nPrimaryVtx");
    vars.push_back("mu");
    vars.push_back("Mll01");
    vars.push_back("lep_0_phi");
    vars.push_back("lep_0_eta");
    vars.push_back("lep_0_pt");
    vars.push_back("lep_0_charge");
    vars.push_back("lep_1_phi");
    vars.push_back("lep_1_eta");
    vars.push_back("lep_1_pt");
    vars.push_back("lep_1_charge");

    // Skim the data
    analyzer.MakeSkimming(precut,vars);             

    return;
}
