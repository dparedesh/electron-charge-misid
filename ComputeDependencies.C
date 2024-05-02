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

std::vector<int> col={2,4,8,kOrange+1,9,46,13};
std::vector<TString> leg_pt={"Lepton p_{T} < 60 GeV","Lepton p_{T} #in [60,90] GeV","Lepton p_{T} #in [90,130] GeV","Lepton p_{T} > 130 GeV"};

class qRatios{

   public:

       qRatios(std::vector<TString> files,TString process,TString var,TString normalize_sample,TString ylabel,TString folder);
       ~qRatios();

       void Execute();
       void GetHistos1D();
       void GetHistos2D();
       void Plot1D();
       void GetRatios();

   private:
       TString GetValue(TString str,TString out);
       TString GetRatioLegend(TString name);
       TString GetLegend(TH2F *pH,int bin);
       std::vector<TString> m_names;
       TString m_process;
       TString m_var;
       TString ylabel;
       std::map<TString,TH2F*> m_histos2D;
       std::map<TString,std::vector<TH1D*> > m_histos1D;
       TString m_sample_to_normalization;
       int nPtbins;
       std::map<TString,std::vector<TString> > m_histo_legend;
       TString m_folder;

};

void ComputeDependencies(){

  gROOT->LoadMacro("MyStyle.C");
  SetMyStyle();

  TString bin="Default"; //plus

  TString process="ttbar";

  std::vector<TString> Ht_dep;
  std::vector<TString> Bjets_dep;
  std::vector<TString> Met_dep;
  std::vector<TString> njets_dep;

  if (process=="Zjets"){

     Ht_dep={"Zjets_"+bin+"_Ht_ge0_le400","Zjets_"+bin+"_Ht_ge400_le600","Zjets_"+bin+"_Ht_ge800_le1000","Zjets_"+bin+"_Ht_ge1000","Zjets_"+bin+"_nJets_ge0"};
     Bjets_dep={"Zjets_"+bin+"_BJets_le1","Zjets_"+bin+"_BJets_ge2","Zjets_"+bin+"_nJets_ge0"};
     Met_dep={"Zjets_"+bin+"_met_ge0_le80","Zjets_"+bin+"_met_ge80","Zjets_"+bin+"_nJets_ge0"};
     njets_dep={"Zjets_"+bin+"_nJets_le4","Zjets_"+bin+"_nJets_ge5","Zjets_"+bin+"_nJets_ge0"};
  }
  else { 
     Ht_dep={"ttbar_"+bin+"_Ht_ge0_le500","ttbar_"+bin+"_Ht_ge500","ttbar_"+bin+"_nJets_ge0"};
     Bjets_dep={"ttbar_"+bin+"_BJets_le1","ttbar_"+bin+"_BJets_ge2","ttbar_"+bin+"_nJets_ge0"};
     Met_dep={"ttbar_"+bin+"_met_ge0_le80","ttbar_"+bin+"_met_ge80","ttbar_"+bin+"_nJets_ge0"};
     njets_dep={"ttbar_"+bin+"_nJets_le4","ttbar_"+bin+"_nJets_ge5","ttbar_"+bin+"_nJets_ge0"};
  }

  TString norm="#epsilon_{QmisID}("+bin+")";
  if (bin=="plus") norm="#epsilon_{QmisID}(default(+))";

  TString sample=process+"_"+bin+"_nJets_ge0";//sample used for normalization. 

  /*
  qRatios qRatio_Ht_range(Ht_dep,"ratio_Ht","ratio_Ht_range",sample,"#epsilon_{QmisID}(H_{T} dep.)/"+norm,process);
  qRatio_Ht_range.Execute();

  qRatios qRatio_met_range(Met_dep,"Ratio_met","ratio_met_range",sample,"#epsilon_{QmisID}(E_{T}^{miss} dep.)/"+norm,process);
  qRatio_met_range.Execute();

  qRatios qRatio_njets_range(njets_dep,"Ratio_nJets","ratio_njets_range",sample,"#epsilon_{QmisID}(nJets dep.)/"+norm,process);
  qRatio_njets_range.Execute();
  */

  qRatios qRatio_BJets_range(Bjets_dep,"Ratio_BJets","ratio_Bjets_range",sample,"#epsilon_{QmisID}(BJets dep.)/"+norm,process);
  qRatio_BJets_range.Execute();
  
  /*
  TString ratio="plus";
  std::vector<TString> ttbar={"ttbar_"+ratio+"_Ht_ge500","Zjets_"+ratio+"_Ht_ge500"};

  qRatios qRatio_ttbar(ttbar,"ratio_ttbar","ratio_ttbar_range","Zjets_"+ratio+"_Ht_ge500","#epsilon_{t#bar{t}}/#epsilon_{Z+jets}","Summary");
  qRatio_ttbar.Execute();
  */

  return;
}

qRatios::qRatios(std::vector<TString> files,TString process,TString var,TString normalize_sample,TString ylabel,TString folder):
  ylabel(ylabel),
  m_names(files),
  m_process(process),
  m_var(var),
  m_histo_legend(),
  m_histos1D(),
  m_histos2D(),
  nPtbins(0),
  m_folder(folder),
  m_sample_to_normalization(normalize_sample)
{}

qRatios::~qRatios()
{
  /*for (std::map<TString,TH2F*>::iterator it=m_histos2D.begin(); it!=m_histos2D.end(); ++it){
        delete it->second;
   }
 */
}

void qRatios::Execute(){

  GetHistos2D(); //Extract 2D histos

  GetHistos1D(); //Convert 2D histos in 1D (one per pt bin).

  //Plot1D();

  GetRatios();

  return;
}

void qRatios::GetRatios(){

 std::cout << "## Getting ratio plots " << std::endl;

 std::vector<TH1D*> denominator=m_histos1D[m_sample_to_normalization];
 //std::vector<TString> legend=m_histo_legend[m_sample_to_normalization];

 TCanvas *pC = 0; 
 if (nPtbins==4){
     pC=new TCanvas();
     pC->Divide(2,2);
 }
 else if (nPtbins>4){
   pC=new TCanvas(m_var,m_var,1600,800);
   pC->Divide(3,2); 
}
 
for (unsigned i=0; i<nPtbins; i++){

      pC->cd(i+1);

      TLegend *pLeg=new TLegend(0.6,0.6,0.9,0.9);

      int counter=0;
      for (std::map<TString,std::vector<TH1D*> >::iterator it=m_histos1D.begin(); it!=m_histos1D.end(); ++it){
          
            TH1D *pH=(TH1D*)it->second[i]->Clone();
            pH->Divide(denominator[i]); 
            pH->GetYaxis()->SetRangeUser(0.5,3);            
            pH->SetLineColor(col[counter]);   
            pH->SetLineWidth(2);
            pH->GetYaxis()->SetTitle(ylabel);
            //pH->SetMarkerSize(0);
            pH->SetLineStyle(counter+1);
            pH->SetMarkerStyle(24);
            pH->SetMarkerColor(col[counter]);             

            TString local=GetRatioLegend(it->first);

            if (it->first!=m_sample_to_normalization){

               if (counter==0) pH->Draw("pL e1");
               else pH->Draw("L e1 same");
               pLeg->AddEntry(pH,local,"pl");

               if (counter==0) {

                  TLine *pLine=new TLine(pH->GetXaxis()->GetXmin(),1.,pH->GetXaxis()->GetXmax(),1.);
                  pLine->SetLineWidth(2);
                  pLine->Draw("same");

               }// end counter

               counter+=1;
             } // end normalization
      }// end loop

      pLeg->Draw("same");

      MiniTreeAnalyzer newanalyzer;
      newanalyzer.GetATLAS("Internal",0.185,0.88,false,0.05);
      newanalyzer.GetLabel(0.186,0.83,"13 TeV",0.04);
      newanalyzer.GetLabel(0.186,0.78,m_histo_legend[m_sample_to_normalization][i],0.04);

      pC->SaveAs(m_folder+"/"+m_var+".pdf");
 }

 return;
}

void qRatios::GetHistos2D(){

   for (unsigned int i=0; i<m_names.size(); i++){ 

      TString file="Output/"+m_names[i]+".root";
      TString histoname=m_names[i]+"_all_num";

      std::cout << "-- Reading file: " << file << std::endl;

      TFile *pFile = new TFile(file);

      TH2F *pH=(TH2F*)pFile->Get(histoname);
      pH->SetDirectory(0);

      m_histos2D[m_names[i]]=pH;
      pFile->Close();
   }
   
   TCanvas *pC= new TCanvas;
   pC->cd();
   m_histos2D[m_names[0]]->Draw("COLZ,TEXT");

   MiniTreeAnalyzer newanalyzer;
   newanalyzer.GetATLAS("Internal",0.185,0.88,false,0.05);
   newanalyzer.GetLabel(0.186,0.81,"13 TeV",0.04);

   return;
}

void qRatios::GetHistos1D(){
    
   for (std::map<TString,TH2F*>::iterator it=m_histos2D.begin(); it!=m_histos2D.end(); ++it){

        TH2F *pH2=(TH2F*)it->second->Clone();

        nPtbins=pH2->GetNbinsY();

        std::vector<TH1D*> v_pH1D;
        v_pH1D.clear();

        std::vector<TString> v_legend;

        std::cout<<"- Will project histos from  : " << it->first << std::endl;
        
        for (unsigned int i=1; i <= pH2->GetNbinsY(); i++){
            std::cout << "-- Proyecting histogram in bin : " << i << std::endl;   
            std::string bin = std::to_string(i);

            TString name=it->first+"_bin"+bin+"_"+m_var; 

            TH1D *pH=(TH1D*)pH2->ProjectionX(name,i,i,"e");
            pH->SetDirectory(0);
            pH->GetYaxis()->SetTitle("#epsilon_{QmisID}");

            TString legend=GetLegend(pH2,i);

            v_legend.push_back(legend);

            std::cout <<"----- DEBUG :" << pH->GetBinContent(1) << std::endl;
            v_pH1D.push_back(pH);

            /*TCanvas *pC = new TCanvas;
            pC->cd();
            pH->Draw();*/
        }

        m_histos1D[it->first]=v_pH1D;
        m_histo_legend[it->first]=v_legend;
   }

   return;
}

TString qRatios::GetRatioLegend(TString name){

   TString leg=name;

   MiniTreeAnalyzer analyzer;
   std::vector<std::string> tokens;

   std::string local(name.Data());

   analyzer.tokenizeString(local,'_',tokens);

   TString str_min="";
   TString str_max="";
   TString var=TString(tokens[2]);

   for (unsigned int i=0; i<tokens.size(); i++){

      TString str=tokens[i];
      if (str.Contains("ge")) str_min=str;
      else if (str.Contains("le")) str_max=str;
   }
 
   TString min=str_min;
   TString max=str_max;

   std::cout <<"## DEBUG :" << min << " max: " << max << std::endl;

   if (str_min!="") min=GetValue(str_min,"ge");
   if (str_max!="") max=GetValue(str_max,"le");

   if (str_min !="" && str_max != "") leg=var+" #in ["+min+","+max+"]";
   else if (str_min !="") leg=var+" #geq "+min;
   else if (str_max !="") leg=var+" #leq "+max;

   if (name.Contains("met") || name.Contains("Ht")) leg=leg+" GeV";

   if (leg.Contains("met")) leg.ReplaceAll("met","E_{T}^{miss}");
   if (leg.Contains("Ht")) leg.ReplaceAll("Ht","H_{T}");

   std::cout << "-- leg :" << leg << std::endl;

   return leg;
}

TString qRatios::GetValue(TString str,TString out){

   str.ReplaceAll(out,"");

   return str;
}

TString qRatios::GetLegend(TH2F *pH,int bin){

  TString legend="p_{T}^{l}";

  int low=pH->GetYaxis()->GetBinLowEdge(bin);
  int up= pH->GetYaxis()->GetBinLowEdge(bin+1);

  if (low==0) legend=legend+" < "+std::to_string(up)+" GeV";
  else if (bin==pH->GetNbinsY()) legend=legend+" > "+std::to_string(low)+" GeV";
  else legend=legend+" #in ["+std::to_string(low)+","+std::to_string(up)+"] GeV";

  return legend;
}

void qRatios::Plot1D(){

  for (std::map<TString,std::vector<TH1D*> >::iterator it=m_histos1D.begin(); it!=m_histos1D.end(); ++it){ 

     std::cout << "- Plotting histos for : " << it->first << std::endl;

     TCanvas *pC = new TCanvas;
     pC->cd();
     pC->SetLogy();
 
     TLegend *pLeg = new TLegend(0.7,0.2,0.9,0.4);
  
     TString temp="";

     for (unsigned int j=0; j<it->second.size(); j++){

         std::cout << "-- In bin : " << j << std::endl;

         TH1D *pH=it->second[j];    
         std::cout <<"----- DEBUG :" << pH->GetBinContent(1) << std::endl;  
   
         pH->GetYaxis()->SetRangeUser(0.0001,1);;

         pH->SetMarkerStyle(20);
         pH->SetMarkerColor(col[j]);
         pH->SetLineWidth(2);
         pH->SetLineColor(col[j]);
         if (j==0) pH->Draw("L e");
         else pH->Draw("L e same");

         pLeg->AddEntry(pH,m_histo_legend[it->first][j],"l");

         pH->GetName();

     }//end for
     pLeg->Draw("same");

     TString name="mc16d: Z+jets";
     if (temp.Contains("ttbar")) name="mc16d: t#bar t"; 

     MiniTreeAnalyzer newanalyzer;
     newanalyzer.GetATLAS("Internal",0.185,0.88,false,0.05);
     newanalyzer.GetLabel(0.186,0.83,"13 TeV",0.04);
     newanalyzer.GetLabel(0.186,0.78,name,0.04);

     pC->SaveAs("Rates1D/"+it->first+".pdf");

  }

  return;
}

