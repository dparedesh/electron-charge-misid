#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <vector>
#include "TCanvas.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TString.h"
#include <map>
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TGraphAsymmErrors.h"



#include "BaselineFramework/Tools/Yields.C"
#include "BaselineFramework/Tools/Channel.C"
#include "BaselineFramework/Tools/EventCut.C"
#include "BaselineFramework/Tools/VariableDistr2D.C"
#include "BaselineFramework/Tools/VariableDistr.C"
#include "BaselineFramework/Tools/PhysicsSample.C"
#include "BaselineFramework/Tools/PhysicsProcess.C"

#include "BaselineFramework/Tools/MiniTreeAnalyzer.C"

#include "MyStyle.C"

class qHisto{

     public:
         qHisto(TH2F* pH);
         ~qHisto();
         inline float GetMin(){return min;};
         inline float GetMax(){return max;};
         inline TH2F* GetHisto(){return m_pH;};

         void SetValues();

     private:

         float GetValue(TString str,TString out);
         float min;
         float max;
         TH2F *m_pH;

};



class qWeight{

    public:
      qWeight(float w,float up=0,float d=0);
      ~qWeight();

      inline void SetUp(float w){w_up=w;};
      inline void SetDown(float w){w_down=w;};
      inline float GetUp(){return w_up;};
      inline float GetDown(){return w_down;};
      inline float GetNominal(){return w_nom;};


    private:
      float w_nom;
      float w_up;
      float w_down;


};


class AddQmisIDWeight{

  public:
   AddQmisIDWeight(TString input,TString path,TString tree,TString inputRates);
   ~AddQmisIDWeight();   

   void Execute();
   void AddDependence(TString name,std::vector<TString> rates);

  private:

    void ReadInputFile(TString input);
    qWeight* ComputeWeight(TString ch,float el_pt1,float el_pt2,float el_eta1,float el_eta2,std::vector<qHisto*> local_rates,float var=-1);
    TH2F* ChooseHisto(std::vector<qHisto*> histos,float var);
    std::vector<double> GetEpsilon(TH2F *pH,float el_pt,float el_eta);
    double GetWeight(double e1,double e2);
    double GetWeight_em(double e);
    std::vector<qHisto*> GetHistos2D(std::vector<TString> rates);

    std::map<TString,std::vector<qHisto*>> m_dependence;
    TString m_path;
    TString m_tree;
    TString m_input;
    TString m_inputRates;

};
void GetQMisIDNtuple(){


   TString path="/eos/user/d/dparedes/ChargeFlipTest/";   
   TString tree="nominal_Loose";
   TString inputRates="RatesScaled/";

   vector<TString> conf_300={"Zjets_plus_Ht_ge0_le300_all_num","Zjets_plus_Ht_ge300_le600_all_num","Zjets_plus_Ht_ge600_le900_all_num","Zjets_plus_Ht_ge900_all_num"};
   vector<TString> conf_400={"Zjets_plus_Ht_ge0_le400_all_num","Zjets_plus_Ht_ge400_le800_all_num","Zjets_plus_Ht_ge800_all_num"}; 
   vector<TString> conf_500={"Zjets_plus_Ht_ge0_le500_all_num","Zjets_plus_Ht_ge500_all_num"}; 
   vector<TString> standard={"NoFakes_plus_nJets_ge0_all_num"};



   AddQmisIDWeight qMisID("data18.root",path,tree,inputRates);

   qMisID.AddDependence("Standard",standard);
   qMisID.AddDependence("Ht_300",conf_300);
   qMisID.AddDependence("Ht_400",conf_400);
   qMisID.AddDependence("Ht_500",conf_500);

   qMisID.Execute();


 return;
}
void AddQmisIDWeight::Execute(){

  ReadInputFile(m_input); 



 return;
}
AddQmisIDWeight::AddQmisIDWeight(TString input,TString path,TString tree,TString inputRates):
 m_dependence(),
 m_path(path),
 m_inputRates(inputRates),
 m_tree(tree),
 m_input(input)
{
}


AddQmisIDWeight::~AddQmisIDWeight()
{}

void AddQmisIDWeight::ReadInputFile(TString input){

  double toGeV=0.001;

  std::cout << "-- Reading input file: " << input << std::endl;

  TFile *pFile = new TFile(m_path+input,"update");

  TTree *pTree=(TTree*)pFile->Get(m_tree);

  int nEvents=pTree->GetEntries();

  std::cout << "-- Total events: " << nEvents << std::endl;

  //Get variables
  std::vector<Float_t> *el_pt=new std::vector<Float_t>();
  std::vector<Float_t> *el_eta=new std::vector<Float_t>();
  std::vector<Char_t> *el_ECIDS=new std::vector<Char_t>();
  std::vector<Float_t> *mu_pt=new std::vector<Float_t>();
  std::vector<Float_t> *mu_eta=new std::vector<Float_t>();

  float HT_all=-1;
  int OSee=-10;
  int OSem=-10;


  float weight_qFlip_standard_nom=-1;
  float weight_qFlip_Ht300_nom=-1;
  float weight_qFlip_Ht400_nom=-1;
  float weight_qFlip_Ht500_nom=-1;

  float weight_qFlip_standard_up=-1;
  float weight_qFlip_Ht300_up=-1;
  float weight_qFlip_Ht400_up=-1;
  float weight_qFlip_Ht500_up=-1;

  float weight_qFlip_standard_down=-1;
  float weight_qFlip_Ht300_down=-1;
  float weight_qFlip_Ht400_down=-1;
  float weight_qFlip_Ht500_down=-1;


  TBranch *b_standard_nom=pTree->Branch("weight_qFlip_standard_nom",&weight_qFlip_standard_nom,"weight_qFlip_standard_nom/F");

  TBranch *b_Ht300_nom=pTree->Branch("weight_qFlip_Ht300_nom",&weight_qFlip_Ht300_nom,"weight_qFlip_Ht300_nom/F");
  TBranch *b_Ht400_nom=pTree->Branch("weight_qFlip_Ht400_nom",&weight_qFlip_Ht400_nom,"weight_qFlip_Ht400_nom/F");
  TBranch *b_Ht500_nom=pTree->Branch("weight_qFlip_Ht500_nom",&weight_qFlip_Ht500_nom,"weight_qFlip_Ht500_nom/F");

  
  TBranch *b_standard_up=pTree->Branch("weight_qFlip_standard_up",&weight_qFlip_standard_up,"weight_qFlip_standard_up/F");

  TBranch *b_Ht300_up=pTree->Branch("weight_qFlip_Ht300_up",&weight_qFlip_Ht300_up,"weight_qFlip_Ht300_up/F");
  TBranch *b_Ht400_up=pTree->Branch("weight_qFlip_Ht400_up",&weight_qFlip_Ht400_up,"weight_qFlip_Ht400_up/F");
  TBranch *b_Ht500_up=pTree->Branch("weight_qFlip_Ht500_up",&weight_qFlip_Ht500_up,"weight_qFlip_Ht500_up/F");


  TBranch *b_standard_down=pTree->Branch("weight_qFlip_standard_down",&weight_qFlip_standard_down,"weight_qFlip_standard_down/F");

  TBranch *b_Ht300_down=pTree->Branch("weight_qFlip_Ht300_down",&weight_qFlip_Ht300_down,"weight_qFlip_Ht300_down/F");
  TBranch *b_Ht400_down=pTree->Branch("weight_qFlip_Ht400_down",&weight_qFlip_Ht400_down,"weight_qFlip_Ht400_down/F");
  TBranch *b_Ht500_down=pTree->Branch("weight_qFlip_Ht500_down",&weight_qFlip_Ht500_down,"weight_qFlip_Ht500_down/F");


  pTree->SetBranchAddress("el_pt",&el_pt);
  pTree->SetBranchAddress("el_eta",&el_eta);
  pTree->SetBranchAddress("mu_pt",&mu_pt);
  pTree->SetBranchAddress("mu_eta",&mu_eta);


  pTree->SetBranchAddress("OSee",&OSee);
  pTree->SetBranchAddress("OSem",&OSem);

  pTree->SetBranchAddress("HT_all",&HT_all);
  pTree->SetBranchAddress("el_ECIDS",&el_ECIDS);

   for (unsigned int i=0; i<nEvents; i++){

        pTree->GetEntry(i);

        if (i%10000 ==0) std::cout << "-- event: " << i << std::endl;

        float el_pt1=(*el_pt)[0]*toGeV;
        float el_pt2=(*el_pt)[1]*toGeV;
        float el_eta1=fabs((*el_eta)[0]);
        float el_eta2=fabs((*el_eta)[1]);
        float Ht=HT_all*toGeV;
        char m_OSee=OSee;
        char m_OSem=OSem;
        char el_BDT1=(*el_ECIDS)[0];
        char el_BDT2=(*el_ECIDS)[1];
        //if (!el_BDT1 || !el_BDT2) continue;
     
        qWeight* chWeight_300=0;
        qWeight* chWeight_400=0;
        qWeight* chWeight_500=0;
        qWeight *chWeight_standard=0;


        if ((m_OSee==1 && el_BDT2==1 && el_BDT1==1) || (m_OSem==1 && el_BDT1==1) ){

          TString ch="ee";
          if (m_OSem==1) ch="em";

          chWeight_standard=ComputeWeight(ch,el_pt1,el_pt2,el_eta1,el_eta2,m_dependence["Standard"]);
          chWeight_300=ComputeWeight(ch,el_pt1,el_pt2,el_eta1,el_eta2,m_dependence["Ht_300"],Ht);
          chWeight_400=ComputeWeight(ch,el_pt1,el_pt2,el_eta1,el_eta2,m_dependence["Ht_400"],Ht);
          chWeight_500=ComputeWeight(ch,el_pt1,el_pt2,el_eta1,el_eta2,m_dependence["Ht_500"],Ht);


          weight_qFlip_standard_nom=chWeight_standard->GetNominal();
          weight_qFlip_Ht300_nom=chWeight_300->GetNominal();
          weight_qFlip_Ht400_nom=chWeight_400->GetNominal();
          weight_qFlip_Ht500_nom=chWeight_500->GetNominal();;


          weight_qFlip_standard_up=chWeight_standard->GetUp();
          weight_qFlip_Ht300_up=chWeight_300->GetUp();
          weight_qFlip_Ht400_up=chWeight_400->GetUp();
          weight_qFlip_Ht500_up=chWeight_500->GetUp();


          weight_qFlip_standard_down=chWeight_standard->GetDown();
          weight_qFlip_Ht300_down=chWeight_300->GetDown();
          weight_qFlip_Ht400_down=chWeight_400->GetDown();
          weight_qFlip_Ht500_down=chWeight_500->GetDown();

          


        }///endof OS reweighting
        else {
          weight_qFlip_Ht300_nom=0;
          weight_qFlip_Ht400_nom=0;
          weight_qFlip_Ht500_nom=0;
          weight_qFlip_standard_nom=0;

          weight_qFlip_Ht300_up=0;
          weight_qFlip_Ht400_up=0;
          weight_qFlip_Ht500_up=0;
          weight_qFlip_standard_up=0;

          weight_qFlip_Ht300_down=0;
          weight_qFlip_Ht400_down=0;
          weight_qFlip_Ht500_down=0;
          weight_qFlip_standard_down=0;
        }


        //if (m_OSem==1) std::cout << OSee << " " << OSem << " , weight:" << weight_qFlip_Ht300_nom << std::endl;

        b_standard_nom->Fill();
        b_standard_up->Fill();
        b_standard_down->Fill();

        b_Ht300_nom->Fill();
        b_Ht300_up->Fill();
        b_Ht300_down->Fill();

        b_Ht400_nom->Fill();
        b_Ht400_up->Fill();
        b_Ht400_down->Fill();

        b_Ht500_nom->Fill();
        b_Ht500_up->Fill();
        b_Ht500_down->Fill();      

    }// end loop over eventes


    pTree->Print(); 
    pTree->Write(); 

    pFile->Close();
    //delete pFile;



  return;
}

qWeight::qWeight(float w,float up=0,float d=0):
 w_nom(w),
 w_up(up),
 w_down(d)
{}
qWeight::~qWeight()
{}

void AddQmisIDWeight::AddDependence(TString name,std::vector<TString> rates){

  std::vector<qHisto*> histos2D=GetHistos2D(rates);

  m_dependence[name]=histos2D;

 return;
}

std::vector<qHisto*> AddQmisIDWeight::GetHistos2D(std::vector<TString> rates){

    std::vector<qHisto*> local;

    for (unsigned int i=0; i<rates.size(); i++){

      TString file=m_inputRates+"/"+rates[i]+".root";
      TString histoname=rates[i];

      std::cout << "-- Reading file with rates: " << rates[i] << std::endl;

      TFile *pFile = new TFile(file);

      TH2F *pH=(TH2F*)pFile->Get(histoname);
      pH->SetDirectory(0);

      qHisto *local_histo=new qHisto(pH);

      local.push_back(local_histo);

      pFile->Close();
   }


   return local;
}

TH2F* AddQmisIDWeight::ChooseHisto(std::vector<qHisto*> histos,float var){

  TH2F *pH=0;


  for (unsigned int i=0; i<histos.size(); i++){

     float min=histos[i]->GetMin();
     float max=histos[i]->GetMax();

     if (min >=0 && max >=0){
        if (var >= min && var <= max) pH=histos[i]->GetHisto();
     }
     else if (min>=0 && max <0){
        if (var >= min) pH=histos[i]->GetHisto();
     }
     else if (min<0 && max >0){
        if (var <= max) pH=histos[i]->GetHisto();
     }
     else std::cout << "-- Undefined values for min and max  : " << min << ", " << max << std::endl;

  }



  return pH;
}

std::vector<double> AddQmisIDWeight::GetEpsilon(TH2F *pH,float el_pt,float el_eta){

  int bin=pH->FindBin(el_eta,el_pt);

  double epsilon=pH->GetBinContent(bin);

  double up=epsilon+pH->GetBinError(bin);
  double down=epsilon-pH->GetBinError(bin);

  if (down<0) down=0;
  if (up>1) up=1;

  std::vector<double> eps={epsilon,up,down};

  return eps;
}

qWeight* AddQmisIDWeight::ComputeWeight(TString ch,float el_pt1,float el_pt2,float el_eta1,float el_eta2,std::vector<qHisto*> local_rates,float var=-1){

   double weight=-1;
   double up=-1;
   double down=-1;

   TH2F *pH=0;

   if (local_rates.size()==1) pH=local_rates[0]->GetHisto();
   else pH=ChooseHisto(local_rates,var);

   std::vector<double> epsilon_1=GetEpsilon(pH,el_pt1,el_eta1);
  

   if (ch=="ee"){

      std::vector<double> epsilon_2=GetEpsilon(pH,el_pt2,el_eta2);
   
      weight=GetWeight(epsilon_1[0],epsilon_2[0]);
      up=GetWeight(epsilon_1[1],epsilon_2[1]);
      down=GetWeight(epsilon_1[2],epsilon_2[2]);

   }
   else if (ch=="em"){

      //std::cout <<"--> Provided channeld is em: getting weight now " << std::endl;
      weight=GetWeight_em(epsilon_1[0]);
      up=GetWeight_em(epsilon_1[1]);
      down=GetWeight_em(epsilon_1[2]);
   }
   else {
     std::cout <<"--> Provided channeld does not match any of the needed ones. Something went wrong!!! Aborting!! " << std::endl;
     exit(1);
   }

   qWeight *total=new qWeight(weight,up,down);

   return total;
}
double AddQmisIDWeight::GetWeight_em(double e1){

  double num=e1;
  double den=1-e1;

  double weight=num/den;

  //std::cout << "-- weight: " << weight << std::endl;
  return weight;
}

double AddQmisIDWeight::GetWeight(double e1,double e2){

  double num=e1+e2-2*e1*e2;
  double den=1-e1-e2+2*e1*e2;

  double weight=num/den;

  return weight;
}
qHisto::qHisto(TH2F* pH):
 m_pH(pH),
 min(-1),
 max(-1)
{

  SetValues();

}
qHisto::~qHisto()
{}
void qHisto::SetValues(){


 std::string name(m_pH->GetName());
 std::cout << "---- in internal class : " << name << std::endl;

 MiniTreeAnalyzer analyzer;
 std::vector<std::string> tokens;

 analyzer.tokenizeString(name,'_',tokens);


 TString str_min="";
 TString str_max="";

 for (unsigned int i=0; i<tokens.size(); i++){

     TString l_s=tokens[i];

      TString str( l_s(0,2) );
      if (str=="ge") str_min=tokens[i];
      else if (str=="le") str_max=tokens[i];

 }

 if (str_min!="") min=GetValue(str_min,"ge");
 if (str_max!="") max=GetValue(str_max,"le");

 std::cout << "- min: " << min << ", max: " << max << std::endl;

 return;
}
float qHisto::GetValue(TString str,TString out){

   float val=-1;

   str.ReplaceAll(out,"");

   val=str.Atof();

   return val;
}



