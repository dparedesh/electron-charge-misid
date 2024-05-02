#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <vector>
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TString.h"
#include <map>
#include "TFile.h"
#include "TLine.h"

#include "../../BaselineFramework/Tools/Yields.C"
#include "../../BaselineFramework/Tools/Channel.C"
#include "../../BaselineFramework/Tools/EventCut.C"
#include "../../BaselineFramework/Tools/VariableDistr2D.C"
#include "../../BaselineFramework/Tools/VariableDistr.C"
#include "../../BaselineFramework/Tools/PhysicsSample.C"
#include "../../BaselineFramework/Tools/PhysicsProcess.C"

#include "../../BaselineFramework/Tools/MiniTreeAnalyzer.C"

#include "MyStyle.C"

std::vector<int> col={2,4,8,kOrange+1,kMagenta,9,13};
void GetHistos1D(TH2F *pH2,TString process,TString label);

void ComputeRatio(TString file,TString process);

TString GetLegend(TH2F *pH,int bin);

void CompareProcess(TString num,TString den,TString process);

void CompareProcesses(){

  gROOT->LoadMacro("MyStyle.C");
  SetMyStyle(); 

  ComputeRatio("Zjets_Default_nJets_ge0.root","Z+jets");  
  //ComputeRatio("Zjets_plus_nJets_ge0.root","Z+jets");   

  ComputeRatio("ttbar_Default_nJets_ge0.root","t#bar{t}");
  //ComputeRatio("ttbar_plus_nJets_ge0.root","t#bar{t}");   

  CompareProcess("Zjets_Default_nJets_ge0.root","ttbar_Default_nJets_ge0.root","wrongTrack");
  //CompareProcess("Zjets_plus_nJets_ge0.root","ttbar_plus_nJets_ge0.root","wrongTrack");

  CompareProcess("Zjets_Default_nJets_ge0.root","ttbar_Default_nJets_ge0.root","photon_conversion");
  //CompareProcess("Zjets_plus_nJets_ge0.root","ttbar_plus_nJets_ge0.root","photon_conversion");

  CompareProcess("Zjets_Default_nJets_ge0.root","ttbar_Default_nJets_ge0.root","all");
  //CompareProcess("Zjets_plus_nJets_ge0.root","ttbar_plus_nJets_ge0.root","all");

  return;
}

void CompareProcess(TString num,TString den,TString process){

  TFile *pNum= new TFile("Output/"+num);
  TFile *pDen= new TFile("Output/"+den);
 
  TString name_num=num;
  name_num.ReplaceAll(".root","");

  TString name_den=den;
  name_den.ReplaceAll(".root","");

  TH2F *pH_num= (TH2F*)pNum->Get(name_num+"_"+process+"_num");
  TH2F *pH_den= (TH2F*)pDen->Get(name_den+"_"+process+"_num");

  TH2F *pRatio=(TH2F*)pH_num->Clone();
  pRatio->SetDirectory(0);

  pRatio->Divide(pH_den);
  
  GetHistos1D(pRatio,process,"#epsilon_{Z+jets}/#epsilon_{t#bar{t}}");

  return;
}
void ComputeRatio(TString file,TString process){

  TString name=file;
  name.ReplaceAll(".root","");

  TFile *pFile = new TFile("Output/"+file);

  TH2F *pH_track= (TH2F*)pFile->Get(name+"_wrongTrack_num");
  TH2F *pH_photonConv= (TH2F*)pFile->Get(name+"_photon_conversion_num"); 

  TH2F *pRatio=(TH2F*)pH_track->Clone();
  pRatio->SetDirectory(0);

  pRatio->Divide(pH_photonConv);
  
  GetHistos1D(pRatio,process,"#epsilon_{wrongTrack}/#epsilon_{photonConversion}");
 
  //pRatio->Draw("COLZ,TEXT");
  
  //pFile->Close();

  return;
}
void GetHistos1D(TH2F *pH2,TString process,TString label){

  int nPtbins=pH2->GetNbinsY();

  std::vector<TH1D*> v_pH1D;
  v_pH1D.clear();

  std::vector<TString> v_legend;

  TCanvas *pC = new TCanvas;      
  pC->cd();

  TLegend *pLeg= new TLegend(0.6,0.6,0.9,0.9);
  TString toSave="";

  TString rename=process;

  for (unsigned int i=1; i <= pH2->GetNbinsY(); i++){

    std::cout << "-- Proyecting histogram in bin : " << i << std::endl;
    std::string bin = std::to_string(i);

    TString name="_bin"+bin;

    TH1D *pH=(TH1D*)pH2->ProjectionX(name,i,i,"e");
    pH->SetDirectory(0);
    pH->GetYaxis()->SetTitle(label);

    TString legend=GetLegend(pH2,i);

    v_legend.push_back(legend);

    std::cout <<"----- DEBUG :" << pH->GetBinContent(1) << std::endl;
    v_pH1D.push_back(pH);

    pH->GetYaxis()->SetRangeUser(0,3.5);
 
    if (process=="wrongTrack" || process=="photon_conversion" || process=="all") {

      pH->GetYaxis()->SetRangeUser(0.5,2.2);

      TString h_name=pH2->GetName();

      if (h_name.Contains("Default")) toSave="Default";
      else if (h_name.Contains("plus")) toSave="plus";

      toSave=process+"_"+toSave;
 
      std::cout << "-- In : " << toSave << std::endl;
             
    }
    else {

      toSave=pH2->GetName();
      toSave.ReplaceAll("_nJets_ge0_wrongTrack_num","");
    }


    if (process=="wrongTrack") rename="Wrong track";
    else if (process=="photon_conversion") rename="Photon conv.";
    else if (process=="all") rename="";


    pH->SetLineColor(col[i-1]);
    pH->SetLineWidth(2);
    //pH->SetMarkerSize(0);
    pH->SetLineStyle(i);
    pH->SetMarkerStyle(24);
    pH->SetMarkerColor(col[i-1]);

    pLeg->AddEntry(pH,legend,"pl");

    if (i==1) pH->Draw("pL e1");
    else pH->Draw("L e1 same");

    if (i==1) {

      TLine *pLine=new TLine(pH->GetXaxis()->GetXmin(),1.,pH->GetXaxis()->GetXmax(),1.);
      pLine->SetLineWidth(2);
      pLine->Draw("same");

    }


  }

  pLeg->Draw("same");
         
  MiniTreeAnalyzer newanalyzer;
  newanalyzer.GetATLAS("Internal",0.185,0.88,false,0.05);
  newanalyzer.GetLabel(0.186,0.83,"13 TeV",0.04);
  newanalyzer.GetLabel(0.186,0.78,rename,0.04);

  pC->SaveAs("Plots-2019-03-19/"+toSave+".pdf");

  return;
}

TString GetLegend(TH2F *pH,int bin){

  TString legend="p_{T}^{l}";

  int low=pH->GetYaxis()->GetBinLowEdge(bin);
  int up= pH->GetYaxis()->GetBinLowEdge(bin+1);

  if (low==0) legend=legend+" < "+std::to_string(up)+" GeV";
  else if (bin==pH->GetNbinsY()) legend=legend+" > "+std::to_string(low)+" GeV";
  else legend=legend+" #in ["+std::to_string(low)+","+std::to_string(up)+"] GeV";

  return legend;
}
