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
#include "TTree.h"
#include "TLine.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"

#include "BaselineFramework/Tools/Yields.C"
#include "BaselineFramework/Tools/Channel.C"
#include "BaselineFramework/Tools/EventCut.C"
#include "BaselineFramework/Tools/VariableDistr2D.C"
#include "BaselineFramework/Tools/VariableDistr.C"
#include "BaselineFramework/Tools/PhysicsSample.C"
#include "BaselineFramework/Tools/PhysicsProcess.C"

#include "BaselineFramework/Tools/MiniTreeAnalyzer.C"

#include "MyStyle.C"

std::vector<int> col={13,kRed,kAzure+1,kOrange+1,kViolet+2,kGreen+3};
//std::vector<int> col={2,8,kOrange+1,6,46,13};
std::vector<int> marker{20,22,23,24,25,26,27,28};
//std::vector<int> fillStyle={0,0,0,0,0};
std::vector<int> fillStyle={3004,3005,3004,3005,3004,3005,3004};

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

class qValidate{

   public:

        qValidate(TString filelist,TString tree);
        ~qValidate();

        std::vector<TString> FillVector(TString file);
        void Execute();
        inline void SetTreeName(TString tree){m_tree=tree;};
        inline void SetPathToFile(TString path){m_path=path;};
        void AddDependence(TString name,TString latex,std::map<TString,std::vector<TString> > rates);
        void AddDependence(TString name,TString latex,std::vector<TString> rates);
        inline void AddVariable(VariableDistr *var){v_var.push_back(var);};
        inline void AddBkgComposition(std::map<TString,float> bkgComp){m_bkg_comp=bkgComp;};
        inline void SetTag(TString tag){m_tag=tag;};  
        inline void SetInputRates(TString input){File_wRates=input;};
        inline void SetMCversion(TString mcversion){m_mcversion=mcversion;};
        inline void SetApplyBDT(bool doBDT){m_doBDT=doBDT;};
        inline void SetProcess(TString process){m_process=process;};
        inline void SetOutputDir(TString output){m_outputDir=output;}; 
        inline void SetData(bool isData){m_isData=true;};

   private:

       void InitializeHistos();
       void ReadInputFile(TString input);
       void FillTGraphAsymm();
       double GetWeight(double e1,double e2);
       qWeight* ComputeWeight(float el_pt1,float el_pt2,float el_eta1,float el_eta2,std::vector<qHisto*> local_rates,float var=-1);
       qWeight* GetRateReweighting(float el_pt1,float el_pt2,float el_eta1,float el_eta2,TString dependence,float var=-1);
       std::vector<qHisto*> GetHistos2D(std::vector<TString> rates); 
       void PlotDistributions();
       TH2F* ChooseHisto(std::vector<qHisto*> histos,float var);
       void SetHistoStyle(TH1F* &pH,int col,int width,int style,int marker=1);
       void SetGraphStyle(TGraphAsymmErrors* &pG,int col,int width,int style,int marker,int fill);

       std::vector<double> GetEpsilon(TH2F *pH,float el_pt,float el_eta);
       void GetYieldsEstimation();
       TGraphAsymmErrors* GetTGraphAsym(TH1F *nom,TH1F *up,TH1F *down);
       TGraphAsymmErrors* GetRatio(TH1F *pNum,TGraphAsymmErrors *pDen);
       double GetError(double num,double num_error,double den,double den_error);

       std::vector<VariableDistr*> v_var;
       std::map<TString,TString> latex_dependence;
       std::map<TString,std::vector<qHisto*>> m_dependence;    
       std::map<TString,std::vector<qHisto*>> m_method_A;
       std::map<TString, float> m_bkg_comp;

       TString m_tag;
       TString m_list;
       TString m_tree;
       TString m_mcversion;
       bool m_isData;
       bool m_doBDT;
       TString m_path;
       TString m_outputDir;
       TString m_process;
       std::map<TString,TH1F*> m_histo_SS;
       std::map<TString,TH1F*> m_histo_OS;     
       std::map<TString,TH1F*> m_histo_OSUp;
       std::map<TString,TH1F*> m_histo_OSDown;
       std::map<TString,TGraphAsymmErrors*> m_graph_OS;
       float Zjets_comp;
       float ttbar_comp;
       TString File_wRates;
};

void MethodABatchFull(TString process,TString inputRates,TString MCversion,TString var,TString outputDir,TString TestHt){

    gROOT->LoadMacro("MyStyle.C");
    SetMyStyle();

   /*
     VariableDistr *v_pt1= new VariableDistr("el_pt[0]*0.001","PtLep","Leading lepton p_{T} [GeV]","Events /",13,0,1300,false);
    VariableDistr *v_pt2= new VariableDistr("el_pt[1]*0.001","PtSublep","Subleading lepton p_{T} [GeV]","Events /",13,0,1300,false);
    VariableDistr *v_Ht= new VariableDistr("HT_all*0.001","Ht","H_{T} [GeV]","Events /",4,0,2400,false);
    VariableDistr *v_Njets= new VariableDistr("nJets","nJets","N_{jets}","Events /",17,3.5,20.5,true);
    VariableDistr *v_met= new VariableDistr("met_met*0.001","met","E_{T}^{miss} [GeV]","Events /",5,0,600,false);
    VariableDistr *v_mu= new VariableDistr("mu","mu","mu","Events /",10,0,100,false);
    VariableDistr *v_bjets=new VariableDistr("nBTags_MV2c10_77","nbjets","nBJets","Events /",9,0.5,9.5,false);
   */
   
    VariableDistr *v_pt1= new VariableDistr("el_pt[0]*0.001","PtLep","Leading lepton p_{T} [GeV]","Events /",20,0,1000,false,false,true);
    VariableDistr *v_pt2= new VariableDistr("el_pt[1]*0.001","PtSublep","Subleading lepton p_{T} [GeV]","Events /",10,0,500,false,false,true);
    VariableDistr *v_Ht= new VariableDistr("HT_all*0.001","Ht","H_{T} [GeV]","Events /",10,0,2000,false,false,true);
    VariableDistr *v_Htjets= new VariableDistr("HT_jets*0.001","Htjets","H_{T}(jets) [GeV]","Events /",10,0,2000,false,false,true);
    VariableDistr *v_Htlep= new VariableDistr("(HT_all-HT_jets)*0.001","Htlep","H_{T}(lep) [GeV]","Events /",10,0,2000,false,false,true);
    VariableDistr *v_Njets= new VariableDistr("nJets","nJets","nJets","Events /",16,-0.5,15.5,true,false,true);
    VariableDistr *v_met= new VariableDistr("met_met*0.001","met","E_{T}^{miss} [GeV]","Events /",15,0,600,false,false,true);
    VariableDistr *v_bjets=new VariableDistr("nBTags_MV2c10_77","nbjets","BJets","Events /",10,-0.5,9.5,false,false,true);
    VariableDistr *v_mu= new VariableDistr("mu","mu","mu","Events /",10,0,100,false,false,false);
    VariableDistr *v_pv=new VariableDistr("nPrimaryVtx","pv","pv","Events /",16,0,80,false,false,false);
    VariableDistr *v_Mll= new VariableDistr("Mll","Mll","M_{ee} [GeV]","Events /",40,71,111,false,false,false);
    VariableDistr *v_eta1= new VariableDistr("el_eta[0]","EtaLep","Leading lepton #eta","Events /",10,-2.5,2.5,false,false,false);


    TString path="/eos/user/d/dparedes/Samples212560/LooseSSML/Skimming/";
    if (MCversion.Contains("mc16")) path="/eos/user/d/dparedes/Samples212560/LooseSSML/"+MCversion+"/ge0j/Skimming/";

    if (process.Contains("Data")) path="/eos/user/d/dparedes/Samples_4tops_Offline/Samples212560/Skimming/";
       
    TString bin_confi="plus";

    std::vector<TString> Ht_dep;
    std::vector<TString> njets_dep;
    std::vector<TString> Met_dep;
    std::vector<TString> Htjets_dep;
    std::vector<TString> Htlep_dep;
    std::vector<TString> PV_dep;
    std::vector<TString> mu_dep;

    std::vector<TString> standard;
    std::vector<TString> standard_new;
    std::vector<TString> standard_new_data;

    if (process=="Zjets"){

	standard={"Zjets_Default_nJets_ge0"};
	standard_new={"Zjets_plus_nJets_ge0"};
    
	Ht_dep={"Zjets_"+bin_confi+"_Ht_ge0_le300","Zjets_"+bin_confi+"_Ht_ge300_le600","Zjets_"+bin_confi+"_Ht_ge600_le900","Zjets_"+bin_confi+"_Ht_ge900"};
      	Htjets_dep={"Zjets_"+bin_confi+"_Htjets_ge0_le400","Zjets_"+bin_confi+"_Htjets_ge400"};
      	Htlep_dep= {"Zjets_"+bin_confi+"_Htlep_ge0_le400" ,"Zjets_"+bin_confi+"_Htlep_ge400"};    
      	njets_dep={"Zjets_"+bin_confi+"_nJets_le4","Zjets_"+bin_confi+"_nJets_ge5"};
      	PV_dep={"Zjets_"+bin_confi+"_pv_ge0_le20","Zjets_"+bin_confi+"_pv_ge20"};
      	mu_dep={"Zjets_"+bin_confi+"_mu_ge0_le30","Zjets_"+bin_confi+"_mu_ge30"};
    }
    else if (process=="ttbar"){
        standard_new={"ttbar_"+bin_confi+"_nJets_ge0"};
     	standard={"ttbar_Default_nJets_ge0"};
     	Ht_dep={"ttbar_"+bin_confi+"_Ht_ge0_le300","ttbar_"+bin_confi+"_Ht_ge300_le600","ttbar_"+bin_confi+"_Ht_ge600"};
     	Htjets_dep={"ttbar_"+bin_confi+"_Htjets_ge0_le400","ttbar_"+bin_confi+"_Htjets_ge400"};
     	Htlep_dep= {"ttbar_"+bin_confi+"_Htlep_ge0_le400" ,"ttbar_"+bin_confi+"_Htlep_ge400"};
     	njets_dep={"ttbar_"+bin_confi+"_nJets_le4","ttbar_"+bin_confi+"_nJets_ge5"};
     	PV_dep={"ttbar_"+bin_confi+"_pv_ge0_le20","ttbar_"+bin_confi+"_pv_ge20"};
     	mu_dep={"ttbar_"+bin_confi+"_mu_ge0_le30","ttbar_"+bin_confi+"_mu_ge30"};
   }
   else if (process.Contains("mc16") || process.Contains("Data201")) standard={process+"_"+bin_confi+"_nJets_ge0"};
   else if (process=="Data" || process=="mc"){
      standard_new={"Data_"+bin_confi+"_nJets_ge0"};
      standard={"mc_"+bin_confi+"_nJets_ge0"}; 
    }

    TString label="";
    if (bin_confi=="plus") label="(+)";
    TString tag="Z+jets";
    if (process=="ttbar") tag="t#bar{t}";

    TString process_list=process+".list";

    if (process.Contains("mc")) process_list="Zjets.list";
    process_list="InputFiles/"+process_list;

    qValidate ValidationPlots(process_list,"nominal_Loose");
    ValidationPlots.SetPathToFile(path);
    ValidationPlots.SetTag(tag);
    ValidationPlots.SetInputRates(inputRates);
    ValidationPlots.SetMCversion(MCversion);
    ValidationPlots.SetProcess(process);
    ValidationPlots.SetOutputDir(outputDir);

    if (process.Contains("Data")) ValidationPlots.SetData(true);

    if (inputRates.Contains("wBDT")) ValidationPlots.SetApplyBDT(true);

    //lines below are just for closure in both data and MC.

    if (process=="Data" || process=="mc"){

         ValidationPlots.AddDependence("Standard","MC Lik.",standard);

         ValidationPlots.AddDependence("Def_New","Data Lik.",standard_new);

    }
    else ValidationPlots.AddDependence("Standard","Likelihood",standard);
 

    /* Below definitons are for closure test only for MC.
    ValidationPlots.AddDependence("Def_New","Default (+)",standard_new);
    ValidationPlots.AddDependence("Standard","Default",standard);

    ValidationPlots.AddDependence("Ht_300","H_{T} dep."+label,Ht_dep);
  
    if (TestHt==false){
       ValidationPlots.AddDependence("nJets_dep","nJets dep."+label,njets_dep);
       ValidationPlots.AddDependence("pv_dep","nPrimaryVtx dep."+label,PV_dep);
       ValidationPlots.AddDependence("mu_dep","mu dep."+label,mu_dep);
    }
    else if (TestHt==true){

        ValidationPlots.AddDependence("Htlep_300","H_{T}(lep) dep."+label,Htlep_dep);
        ValidationPlots.AddDependence("Htjets_300","H_{T}(jets) dep."+label,Htjets_dep);
    }*/

    if (var=="pt1") ValidationPlots.AddVariable(v_pt1);
    if (var=="pt2") ValidationPlots.AddVariable(v_pt2);
    if (var=="Mll") ValidationPlots.AddVariable(v_Mll);     
    if (var=="eta1") ValidationPlots.AddVariable(v_eta1);

    if (var=="Ht") ValidationPlots.AddVariable(v_Ht);
    if (var=="Htlep") ValidationPlots.AddVariable(v_Htlep);
    if (var=="Htjets") ValidationPlots.AddVariable(v_Htjets);
    if (var=="nJets") ValidationPlots.AddVariable(v_Njets);
    if (var=="met") ValidationPlots.AddVariable(v_met);
    if (var=="BJets") ValidationPlots.AddVariable(v_bjets);
    if (var=="mu") ValidationPlots.AddVariable(v_mu);
    if (var=="pv") ValidationPlots.AddVariable(v_pv);      

    ValidationPlots.Execute();

   /*    
    std::vector<TString> standard_new={"Zjets_plus_nJets_ge0"};
    std::vector<TString> standard={"Zjets_Default_nJets_ge0"};

    std::vector<TString> Ht_dep_200_new={"Zjets_plus_Ht_ge0_le300","Zjets_plus_Ht_ge300_le600","Zjets_plus_Ht_ge400_le600","Zjets_plus_Ht_ge600_le800","Zjets_plus_Ht_ge800_le1000","Zjets_plus_Ht_ge1000"};
    //std::vector<TString> Ht_dep_200={"Zjets_Default_Ht_ge0_le200","Zjets_Default_Ht_ge200_le400","Zjets_Default_Ht_ge400_le600","Zjets_Default_Ht_ge600_le800","Zjets_Default_Ht_ge800_le1000","Zjets_Default_Ht_ge1000"};

    std::vector<TString> ttbar={"ttbar_Default_Ht_ge0_le500","ttbar_Default_Ht_ge500"};

    std::map<TString,std::vector<TString> > Input_MethodA={{"Zjets",Ht_dep_200_new},{"ttbar",ttbar}};
    std::map<TString,float> bkgPer={{"Zjets",0.097},{"ttbar",0.903}}; //{{"Zjets",1},{"ttbar",0}}; //{{"Zjets",0.097},{"ttbar",0.903}};//{{"Zjets",0.213},{"ttbar",0.787}};

    qValidate ValidationPlots("total.list","nominal_Loose");
    ValidationPlots.SetPathToFile(path);
    ValidationPlots.SetTag("Z+jets,t#bar{t}");

    ValidationPlots.AddDependence("Def_New","Default (+)",standard_new);
    ValidationPlots.AddDependence("Standard","Default",standard);
    ValidationPlots.AddDependence("Ht_200","H_{T} dep.",Ht_dep_200);
    ValidationPlots.AddDependence("Ht_200_new","H_{T} dep. (+)",Ht_dep_200_new);
    //ValidationPlots.AddDependence("methodA_Ht","Method A",Input_MethodA);
    //ValidationPlots.AddBkgComposition(bkgPer);    
     
    ValidationPlots.AddVariable(v_pt1);
    ValidationPlots.AddVariable(v_Ht);
    //ValidationPlots.AddVariable(v_Njets);
    //ValidationPlots.AddVariable(v_met);
    //ValidationPlots.AddVariable(v_bjets);
    //ValidationPlots.AddVariable(v_mu);

    ValidationPlots.Execute();
*/ 


    return;
}
qValidate::qValidate(TString filelist,TString tree):
  m_list(filelist),
  m_tree(tree),
  m_path(""),
  v_var(),
  m_tag(""),
  m_bkg_comp(),
  latex_dependence(),
  File_wRates(""),
  m_method_A(),
  m_mcversion(""),
  m_process(""),
  m_outputDir(""),
  m_doBDT(false),
  m_isData(false),
  m_dependence(),
  m_histo_SS(),
  m_histo_OS(),
  m_histo_OSUp(),
  m_histo_OSDown(),
  m_graph_OS()
{
}

qValidate::~qValidate()
{}

void qValidate::ReadInputFile(TString input){

  double toGeV=0.001;

  std::cout << "-- Reading input file: " << m_path+input << std::endl;

  TFile *pFile = new TFile(m_path+input);

  TTree *pTree=(TTree*)pFile->Get(m_tree);

  int nEvents=pTree->GetEntries();

  std::cout << "-- Total events: " << nEvents << std::endl;

  //setting vars
  std::vector<Float_t> *el_pt=new std::vector<Float_t>();
  std::vector<Float_t> *el_eta=new std::vector<Float_t>();
  std::vector<Int_t> *el_charge=new std::vector<Int_t>();
  std::vector<Char_t> *el_ECIDS=new std::vector<Char_t>();
  std::vector<Float_t> *el_phi=new std::vector<Float_t>();
  std::vector<Float_t> *el_e=new std::vector<Float_t>();

  int nJets=-1;
  int nbjets=-1;
  float HT_all=-1;
  float HT_jets=-1;
  float met_met=-1;
  float mu=-1;
  int runNumber=-999;
  float Mll=-1;

  int OSee=0;
  int loose_SSee=0;
  int nMuons=0;

  float weight_bTagSF_MV2c10_77=1;
  float weight_mc=1;
  float weight_leptonSF=1;
  float weight_pileup=1;
  float weight_jvt=1;
  float weight_normalise=1;
  float weight_indiv_SF_EL_ChargeID=1;
  float weight_indiv_SF_EL_ChargeMisID=1;

  int nPrimaryVtx=-1;

  //getting branches
  pTree->SetBranchAddress("el_pt",&el_pt);
  pTree->SetBranchAddress("el_eta",&el_eta);
  pTree->SetBranchAddress("el_phi",&el_phi);
  pTree->SetBranchAddress("el_e",&el_e);
  pTree->SetBranchAddress("el_charge",&el_charge);
  pTree->SetBranchAddress("HT_all",&HT_all);
  pTree->SetBranchAddress("HT_jets",&HT_jets);
  pTree->SetBranchAddress("nJets",&nJets);
  pTree->SetBranchAddress("met_met",&met_met);
  pTree->SetBranchAddress("mu",&mu);
  pTree->SetBranchAddress("nPrimaryVtx",&nPrimaryVtx);
  pTree->SetBranchAddress("nBTags_MV2c10_77",&nbjets);
  pTree->SetBranchAddress("runNumber",&runNumber);
  pTree->SetBranchAddress("el_ECIDS",&el_ECIDS);
  pTree->SetBranchAddress("OSee",&OSee);
  pTree->SetBranchAddress("loose_SSee",&loose_SSee);
  pTree->SetBranchAddress("nMuons",&nMuons);

  //getting weights 
  if (!m_isData){
    pTree->SetBranchAddress("weight_bTagSF_MV2c10_77",&weight_bTagSF_MV2c10_77);
    pTree->SetBranchAddress("weight_mc",&weight_mc);
    pTree->SetBranchAddress("weight_leptonSF",&weight_leptonSF);
    pTree->SetBranchAddress("weight_pileup",&weight_pileup);
    pTree->SetBranchAddress("weight_jvt",&weight_jvt);
    pTree->SetBranchAddress("weight_normalise",&weight_normalise);
    pTree->SetBranchAddress("weight_indiv_SF_EL_ChargeID",&weight_indiv_SF_EL_ChargeID); 
    pTree->SetBranchAddress("weight_indiv_SF_EL_ChargeMisID",&weight_indiv_SF_EL_ChargeMisID);
   }

  //weight_bTagSF_MV2c10_77*weight_mc*weight_leptonSF*weight_pileup*weight_jvt*weight_normalise

  for (unsigned int i=0; i<nEvents; i++){

        pTree->GetEntry(i);
      
        //std::cout << "-- print pt: " << (*el_charge)[0] << ", ... "<< (*el_charge)[1] << std::endl;    
        if (i%100000 ==0) std::cout << "-- event: " << i << std::endl;

        float el_pt1=(*el_pt)[0]*toGeV;
        float el_pt2=(*el_pt)[1]*toGeV;
        float el_eta1=fabs((*el_eta)[0]);
        float el_eta2=fabs((*el_eta)[1]);
        int el_charge1=(*el_charge)[0];
        int el_charge2=(*el_charge)[1];
        float el_phi1=(*el_phi)[0];
        float el_phi2=(*el_phi)[1];
        float el_e1=(*el_e)[0]*toGeV;
        float el_e2=(*el_e)[1]*toGeV;

        float Ht=HT_all*toGeV;
        float Htjets=HT_jets*toGeV;
        float Htlep=Ht-Htjets;
        float met=met_met*toGeV;
        char el_BDT1=(*el_ECIDS)[0];
        char el_BDT2=(*el_ECIDS)[1];
        float npv=nPrimaryVtx;

       //std::cout << "-- print pt: " << OSee << ", ... "<< loose_SSee << std::endl;
       if (!OSee && !loose_SSee) continue;
       if (nMuons>0) continue;     

       if (m_doBDT){
          if (!el_BDT1 || !el_BDT2) continue;
       }

       TLorentzVector v1, v2;
       v1.SetPtEtaPhiE(el_pt1,(*el_eta)[0],el_phi1,el_e1);
       v2.SetPtEtaPhiE(el_pt2,(*el_eta)[1],el_phi2,el_e2);
       Mll=(v1+v2).M();

       //if (nJets<1) continue; //request at least 1 jet for closure. 

       if (Mll<81 || Mll>101) continue; //DO CLOSURE ONLY FOR Z-events (remember that other data events are fully blinded!).

       if (Ht>0){  //Standard
       //if (Ht>500 && met>40 && nbjets==2 && nJets>=6){  //SR2b
       //if (Ht>500 && met>40 && nbjets==3 && nJets>=5){     //SR3b
       //if (Ht>500 && met>40 && nbjets>=4 && nJets>=6) {//SR4b
  
       float event_weight=1;

       if (!m_isData) {
               event_weight=weight_bTagSF_MV2c10_77*weight_mc*weight_leptonSF*weight_pileup*weight_jvt*weight_normalise*(36184.86*(runNumber == 284500) + 43587.3*(runNumber == 300000) + 45691.0*(runNumber == 310000));

               if (m_doBDT) event_weight*=weight_indiv_SF_EL_ChargeID*weight_indiv_SF_EL_ChargeMisID;
       } 
 
       //std::cout << "-- print pt: " << el_pt1 << ", ... "<< el_pt2 << std::endl;
       //std::cout << "-- print eta: " << el_eta1 << ", ... "<< el_eta2 << std::endl; 
       //if (el_pt2>200) std::cout << "-- print pt: " << el_pt1 << ", ... "<< el_pt2 << std::endl;
       //std::cout << "-- Print mu: " << mu << std::endl;
       //if (Ht>500) std::cout <<"######################################## DEBUG ############## " << std::endl;

       for (unsigned int v=0; v<v_var.size(); v++){

                float var_to_fill=0;
                //std::cout << "-- in var :" << v_var[v]->GetName() << std::endl;

                if (v_var[v]->GetTitle().Contains("PtLep")) var_to_fill=el_pt1;
                else if (v_var[v]->GetTitle().Contains("PtSublep")) var_to_fill=el_pt2;
                else if (v_var[v]->GetTitle().Contains("EtaLep")) var_to_fill=(*el_eta)[0];
                else if (v_var[v]->GetTitle().Contains("Htlep"))  var_to_fill=Htlep;
                else if (v_var[v]->GetTitle().Contains("Htjets")) var_to_fill=Htjets;
                else if (v_var[v]->GetTitle().Contains("Ht"))     var_to_fill=Ht;
                else if (v_var[v]->GetTitle().Contains("nJets"))  var_to_fill=nJets;
                else if (v_var[v]->GetTitle().Contains("pv"))     var_to_fill=nPrimaryVtx;
                else if (v_var[v]->GetTitle().Contains("mu"))     var_to_fill=mu;
                else if (v_var[v]->GetTitle().Contains("met"))    var_to_fill=met;
                else if (v_var[v]->GetTitle().Contains("nbjets")) var_to_fill=nbjets;
                else if (v_var[v]->GetTitle().Contains("Mll"))    var_to_fill=Mll;
                else {
                    std::cout << "## Not sure which variable has to be filled... the code is not automatic: the variables to fill have to be setted by hand! -> Aborting." << std::endl;
                    exit(1);
                }
                //std::cout << "-- var to fill: " << var_to_fill << std::endl;        

                if ( el_charge1*el_charge2 > 0){
                     //std::cout << "SS event -> Mll :" << Mll << std::endl; 

                     m_histo_SS[v_var[v]->GetTitle()+"_SS"]->Fill(var_to_fill,event_weight);
                }
                else if ( el_charge1*el_charge2 < 0 ){
                 
                    for (std::map<TString,std::vector<qHisto*> >::iterator it=m_dependence.begin(); it!=m_dependence.end(); ++it){      
                                        
                        float var_to_pass=0;

                        if (it->first.Contains("Standard")) var_to_pass=var_to_fill;
                        else if (it->first.Contains("Def_New")) var_to_pass=var_to_fill;
                        else if (it->first.Contains("Htjets")) var_to_pass=Htjets;
                        else if (it->first.Contains("Htlep")) var_to_pass=Htlep;
                        else if (it->first.Contains("Ht")) var_to_pass=Ht;
                        else if (it->first.Contains("nJets")) var_to_pass=nJets;
                        else if (it->first.Contains("mu")) var_to_pass=mu;
                        else if (it->first.Contains("pv")) var_to_pass=nPrimaryVtx;
                        else if (it->first.Contains("met")) var_to_pass=met;
                        else if (it->first.Contains("BJets")) var_to_pass=nbjets;
                        //else if (it->first.Contains("methodA")) var_to_pass=var_to_fill;   --> methodA_Ht
                        else {
                            std::cout << "## Not sure which variable has to be passed: Dependence is probably unknown or not added... the code is not automatic: the variables to fill have to be setted by hand! -> Aborting." << std::endl;
                            exit(1);
                        }
                        qWeight* chWeight=0;

                        if (it->first.Contains("methodA")) chWeight=GetRateReweighting(el_pt1,el_pt2,el_eta1,el_eta2,it->first,var_to_pass);
                        else chWeight=ComputeWeight(el_pt1,el_pt2,el_eta1,el_eta2,m_dependence[it->first],var_to_pass);

                        double nom=chWeight->GetNominal();
                        double up=chWeight->GetUp();
                        double down=chWeight->GetDown();

                        m_histo_OS[v_var[v]->GetTitle()+"_OS_"+it->first]->Fill(var_to_fill,event_weight*nom);
                        m_histo_OSUp[v_var[v]->GetTitle()+"_OSUp_"+it->first]->Fill(var_to_fill,event_weight*up);                         
                        m_histo_OSDown[v_var[v]->GetTitle()+"_OSDown_"+it->first]->Fill(var_to_fill,event_weight*down);

                     }//filling OS weighted for different dependencies
                }
                else {
                     std::cout << "-- Electron charge has not be filled: " << el_charge1 << " : " << el_charge2 << std::endl;

                } 

        }//filling histos for different vars  

    }//events cut!

  }//end loop over events

  return;
}

void qValidate::Execute(){

   std::vector<TString> samples=FillVector(m_list);

   if (m_dependence.count("methodA")==1 && (Zjets_comp <0 || ttbar_comp <0)){
     std::cout << " Method A habilitated but background composition not given!! --> Aborting !"<< std::endl;
     exit(1);
    } 

   InitializeHistos();  

   for (unsigned int i=0; i<samples.size(); i++) ReadInputFile(samples[i]);   

   FillTGraphAsymm(); 
  
   PlotDistributions();
    
   GetYieldsEstimation();

   return;
}

void qValidate::FillTGraphAsymm(){

   for (unsigned int v=0; v<v_var.size(); v++){ 

       for (std::map<TString,std::vector<qHisto*> >::iterator it=m_dependence.begin(); it!=m_dependence.end(); ++it){


           m_graph_OS[v_var[v]->GetTitle()+"_OS_"+it->first]=GetTGraphAsym(m_histo_OS[v_var[v]->GetTitle()+"_OS_"+it->first],
                                                                           m_histo_OSUp[v_var[v]->GetTitle()+"_OSUp_"+it->first],
                                                                           m_histo_OSDown[v_var[v]->GetTitle()+"_OSDown_"+it->first]);

       } //method
   }//end var

   return;
}
void qValidate::SetGraphStyle(TGraphAsymmErrors* &pG,int col,int width,int style,int marker,int fill){

   pG->SetLineColor(col);
   pG->SetLineWidth(width);
   pG->SetMarkerStyle(marker);
   pG->SetMarkerColor(col);
   pG->SetFillStyle(fill);
   pG->SetLineStyle(style);
   pG->SetFillColor(col);
 
   return;
}
void qValidate::SetHistoStyle(TH1F* &pH,int col,int width,int style,int marker=1){

   pH->SetLineColor(col);
   pH->SetLineStyle(style);
   pH->SetLineWidth(width);
   pH->SetMarkerColor(col);
   pH->SetMarkerStyle(marker);

   return;
}

void qValidate::GetYieldsEstimation(){

  std::cout<<"-## Getting estimated yieds: " << std::endl;

  for (unsigned int v=0; v<v_var.size(); v++){

	if  (v_var[v]->GetYield()){

	  std::cout << " SS events: " << m_histo_SS[v_var[v]->GetTitle()+"_SS"]->Integral() << std::endl;

                   for (std::map<TString,std::vector<qHisto*> >::iterator it=m_dependence.begin(); it!=m_dependence.end(); ++it){

                     double nom=m_histo_OS[v_var[v]->GetTitle()+"_OS_"+it->first]->Integral();
                     double down=nom-m_histo_OSDown[v_var[v]->GetTitle()+"_OSDown_"+it->first]->Integral();
                     double up=m_histo_OSUp[v_var[v]->GetTitle()+"_OSUp_"+it->first]->Integral()-nom;

		     std::cout<< "OS weighted with '"<< it->first << "' : " << nom << " +/-" << up << "/" << down << std::endl;                            

                   }//loop over dependencies
 
	}// get yield for this distribution

  }// loop over var

  return;
}

TGraphAsymmErrors* qValidate::GetTGraphAsym(TH1F *nom,TH1F *up,TH1F *down){

  TGraphAsymmErrors *pG= new TGraphAsymmErrors();

  int nbins=nom->GetNbinsX();

  for (unsigned int i=1; i<=nbins; i++){ 

     float x=nom->GetBinCenter(i);
     float y=nom->GetBinContent(i);
     float ex=0.5*nom->GetBinWidth(i);
     float y_up=up->GetBinContent(i)-y;
     float y_down=y-down->GetBinContent(i);

     pG->SetPoint(i-1,x,y);
     pG->SetPointError(i-1,ex,ex,y_down,y_up);
  } 

  return pG;
}

void qValidate::PlotDistributions(){

   for (unsigned int v=0; v<v_var.size(); v++){

       TCanvas *pC = new TCanvas("h_"+v_var[v]->GetTitle(),"h_"+v_var[v]->GetTitle(),700,600);

       TLegend *pLeg= new TLegend(0.6,0.65,0.9,0.9); //TLegend(0.5,0.5,0.9,0.9);

       TPad* c_plot = new TPad("plot","",0.,0.3,1.0,1.0);
       c_plot->SetRightMargin(0.08);
       c_plot->SetLeftMargin(0.15);
       c_plot->SetBottomMargin(0.0);
       c_plot->SetTopMargin(0.05);
       c_plot->Draw();

       TPad* c_ratio= new TPad("ratio","",0.,0,1,0.3);
       c_ratio->SetRightMargin(0.08);
       c_ratio->SetLeftMargin(0.15);
       c_ratio->SetBottomMargin(0.35);
       c_ratio->SetTopMargin(0.0);
       c_ratio->Draw();
       c_ratio->SetGridy(kTRUE);

       c_plot->cd();
       float plot_min=0;
       float plot_max=0;

       int maxBin=m_histo_SS[v_var[v]->GetTitle()+"_SS"]->GetMaximumBin();
       plot_max=m_histo_SS[v_var[v]->GetTitle()+"_SS"]->GetBinContent(maxBin);

       if (v_var[v]->plotLogY()){
            c_plot->SetLogy();  //forNormal
            plot_min=0.1;
            plot_max=plot_max*10000;

       }
       else plot_max=1.5*plot_max;
      
       m_histo_SS[v_var[v]->GetTitle()+"_SS"]->GetYaxis()->SetRangeUser(plot_min,plot_max);

       SetHistoStyle(m_histo_SS[v_var[v]->GetTitle()+"_SS"],1,2,1);

       m_histo_SS[v_var[v]->GetTitle()+"_SS"]->Draw("hist");      
       pLeg->AddEntry(m_histo_SS[v_var[v]->GetTitle()+"_SS"],"SS","l");
  
       if (m_dependence.count("Standard")==1){

          SetGraphStyle(m_graph_OS[v_var[v]->GetTitle()+"_OS_Standard"],4,2,1,1,3003);
          SetHistoStyle(m_histo_OS[v_var[v]->GetTitle()+"_OS_Standard"],4,2,1);
         
          m_graph_OS[v_var[v]->GetTitle()+"_OS_Standard"]->Draw("2,p");
          m_histo_OS[v_var[v]->GetTitle()+"_OS_Standard"]->Draw("hist,same");
          pLeg->AddEntry(m_graph_OS[v_var[v]->GetTitle()+"_OS_Standard"],"OS weighted("+latex_dependence["Standard"]+")","fl");
       }

       std::vector<TH1F*> ratios;
       std::map<TString,TGraphAsymmErrors*> gRatios;      

       int counter=0;
       for (std::map<TString,std::vector<qHisto*> >::iterator it=m_dependence.begin(); it!=m_dependence.end(); ++it){

             if (!it->first.Contains("Standard")){
                SetGraphStyle(m_graph_OS[v_var[v]->GetTitle()+"_OS_"+it->first],col[counter],2,7,marker[counter],fillStyle[counter]);
	        SetHistoStyle(m_histo_OS[v_var[v]->GetTitle()+"_OS_"+it->first],col[counter],2,7,marker[counter]);

                m_graph_OS[v_var[v]->GetTitle()+"_OS_"+it->first]->Draw("2,p");
                m_histo_OS[v_var[v]->GetTitle()+"_OS_"+it->first]->Draw("hist,same");                
                pLeg->AddEntry(m_graph_OS[v_var[v]->GetTitle()+"_OS_"+it->first],"OS weighted ("+latex_dependence[it->first]+")","lp");
             }

             TH1F *pratio=(TH1F*) m_histo_SS[v_var[v]->GetTitle()+"_SS"]->Clone();
             pratio->Divide(m_histo_OS[v_var[v]->GetTitle()+"_OS_"+it->first]);
             pratio->SetName("Ratio_"+v_var[v]->GetTitle()+"_OS_"+it->first);

             TGraphAsymmErrors *g_pratio=GetRatio(m_histo_SS[v_var[v]->GetTitle()+"_SS"],m_graph_OS[v_var[v]->GetTitle()+"_OS_"+it->first]);
             g_pratio->SetName("gRatio_"+v_var[v]->GetTitle()+"_OS_"+it->first);

             if (!it->first.Contains("Standard")) {
                    SetHistoStyle(pratio,col[counter],2,7,marker[counter]);
                    SetGraphStyle(g_pratio,col[counter],2,7,marker[counter],fillStyle[counter]);       
             }
             else{
                    SetHistoStyle(pratio,4,2,1);
                    SetGraphStyle(g_pratio,4,2,1,1,3003);
             }
             ratios.push_back(pratio);
             gRatios[it->first]=g_pratio;

             if (!it->first.Contains("Standard")) counter+=1;
        }//end of depen
        pLeg->Draw("same");

        MiniTreeAnalyzer newanalyzer;
        newanalyzer.GetATLAS("Internal",0.185,0.88,false,0.055);

        //36184.86*(runNumber == 284500) + 43587.3*(runNumber == 300000) + 45691.0*(runNumber == 310000)
        TString lumi="125 fb^{-1}";
        if (m_mcversion=="mc16a") lumi="36.1 fb^{-1}";
        else if (m_mcversion=="mc16d") lumi="43.6 fb^{-1}"; 
        else if (m_mcversion=="mc16e") lumi="46.7 fb^{-1}";
        else if (m_mcversion.Contains("201")) lumi=m_mcversion;

        newanalyzer.GetLabel(0.186,0.82,"13 TeV, "+lumi,0.045);
        TString local_tag=m_tag+" ";
        if (m_mcversion.Contains("Data")) local_tag="";

        //if (m_doBDT) local_tag+="w/ ECIDS";

        if (m_mcversion.Contains("mc")){ 
             TString mc_tag=m_mcversion;
             if (m_mcversion=="mc") mc_tag="mc16a/d/e: "; 
             newanalyzer.GetLabel(0.186,0.77,mc_tag+":"+local_tag,0.045); 
             if (m_doBDT) newanalyzer.GetLabel(0.186,0.72,"w/ ECIDS",0.045);
        }
        else {
          if (m_doBDT) newanalyzer.GetLabel(0.186,0.77,"w/ ECIDS",0.045);
        }
        c_ratio->cd();       

        int kk=0;

        if (gRatios.count("Standard")==1){

            gRatios["Standard"]->GetYaxis()->SetNdivisions(505);
            gRatios["Standard"]->GetYaxis()->SetRangeUser(0,2.001);
            gRatios["Standard"]->GetXaxis()->SetLabelSize(0.1);
            gRatios["Standard"]->GetXaxis()->SetTitleSize(0.1);
            gRatios["Standard"]->GetYaxis()->SetLabelSize(0.1);
            gRatios["Standard"]->GetYaxis()->SetTitleSize(0.1);
            gRatios["Standard"]->GetYaxis()->SetTitleOffset(0.5);
            gRatios["Standard"]->GetYaxis()->SetTitle("SS/OS");
            gRatios["Standard"]->GetYaxis()->CenterTitle();

            gRatios["Standard"]->Draw("Ap");
            //gRatios["Standard"]->Draw("p");

            kk+=1;
        }

        for (std::map<TString,TGraphAsymmErrors*>::iterator it=gRatios.begin(); it!=gRatios.end(); ++it){

            gRatios[it->first]->GetYaxis()->SetNdivisions(505); //404
            gRatios[it->first]->GetYaxis()->SetRangeUser(0,2.0);
            gRatios[it->first]->GetXaxis()->SetLabelSize(0.1);
            gRatios[it->first]->GetXaxis()->SetTitleSize(0.1);
            gRatios[it->first]->GetYaxis()->SetLabelSize(0.1);
            gRatios[it->first]->GetYaxis()->SetTitleSize(0.1);
            gRatios[it->first]->GetYaxis()->SetTitleOffset(0.5);
            gRatios[it->first]->GetYaxis()->SetTitle("SS/OS");       
            gRatios[it->first]->GetYaxis()->CenterTitle();

            if (!it->first.Contains("Standard")){

                  if (kk==0) gRatios[it->first]->Draw("Ap");
                  else  gRatios[it->first]->Draw("p");

                  kk+=1;
            }

        }

        /*
        for (unsigned int k=0; k<ratios.size(); k++){

            ratios[k]->GetYaxis()->SetRangeUser(0.5,3);
       
            ratios[k]->GetXaxis()->SetLabelSize(0.1);
            ratios[k]->GetXaxis()->SetTitleSize(0.1);
            ratios[k]->GetYaxis()->SetLabelSize(0.1);
            ratios[k]->GetYaxis()->SetTitleSize(0.1);
            ratios[k]->GetYaxis()->SetTitleOffset(0.5);
            ratios[k]->GetYaxis()->SetTitle("Ratio");

            if (k==0) ratios[k]->Draw("l");
            else ratios[k]->Draw("l,same");

            gRatios[k]->Draw("2,p");
        }
        */

        TLine* line1 = new TLine(ratios[0]->GetXaxis()->GetXmin(),1.,ratios[0]->GetXaxis()->GetXmax(),1.);
        line1->SetLineWidth(2);
        line1->Draw("same");

        system("mkdir -p "+m_outputDir);
        TString local_output=m_outputDir+"/h_"+m_process+"_"+v_var[v]->GetTitle();

        if (m_mcversion.Contains("mc16") && !m_process.Contains("mc") && !m_process.Contains("Data")) {
            system("mkdir -p "+m_outputDir+"/"+m_mcversion);
            local_output=m_outputDir+"/"+m_mcversion+"/h_"+m_process+"_"+v_var[v]->GetTitle();
        }

        if (m_doBDT) local_output+="wBDT";

        pC->SaveAs(local_output+".pdf");
        pC->SaveAs(local_output+".root");

   }//end of var

 return;
}

TGraphAsymmErrors* qValidate::GetRatio(TH1F *pNum,TGraphAsymmErrors *pDen){

    TGraphAsymmErrors *pG=new TGraphAsymmErrors();

    if (pDen->GetN()!=pNum->GetNbinsX()){
 
       std::cout << "--> Will try to divide different things for the ratio plot.. aborting!" << std::endl;
       exit(1);
    }

    Double_t* x_d=pDen->GetX();
    Double_t* y_d=pDen->GetY();
    Double_t* y_u=pDen->GetEYhigh();
    Double_t* y_l=pDen->GetEYlow();
    Double_t* x_e=pDen->GetEXlow();

    for (unsigned int i=1; i<=pNum->GetNbinsX(); i++){
  
        //std::cout << "DEBUG: " << pNum->GetBinCenter(i) << ", " << x_d[i-1] << std::endl;

        if (pNum->GetBinCenter(i)==x_d[i-1]){

            //std::cout << "DEBUG: " << pNum->GetBinCenter(i) << ", " << x_d[i-1] << std::endl;

            double num_error=pNum->GetBinError(i);
            double den_error_up=y_u[i-1];
            double den_error_down=y_l[i-1];
            double num=pNum->GetBinContent(i);
            double den=y_d[i-1];         

            double error_up=GetError(num,num_error,den,den_error_up);
            double error_down=GetError(num,num_error,den,den_error_down);
            pG->SetPoint(i-1,x_d[i-1],num/den);
            pG->SetPointError(i-1,x_e[i-1],x_e[i-1],error_down,error_up);
       }
       else {
	    std::cout <<"--> Something wrong happended... aborting!" << std::endl;
       	    exit(1); 
       }
    }

    pG->GetXaxis()->SetTitle(pNum->GetXaxis()->GetTitle());
    pG->GetXaxis()->SetRangeUser(pNum->GetBinLowEdge(1),pNum->GetBinLowEdge(pNum->GetNbinsX()+1));

    return pG;
}

double qValidate::GetError(double num,double num_error,double den,double den_error){

    double term1=pow(num_error/num,2);
    double term2=pow(den_error/den,2);

    double val=(num/den)*TMath::Sqrt(term1+term2);
 
    return val;
}

TH2F* qValidate::ChooseHisto(std::vector<qHisto*> histos,float var){

    TH2F *pH=0;  

    //std::cout << "var: " << var << std::endl;

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

    //std::cout << "-- min: " << min << ", max=" << max << std::endl; 

    //std::cout << "-- Will use histo : " << pH->GetName() << std::endl;

    return pH;
}

std::vector<double> qValidate::GetEpsilon(TH2F *pH,float el_pt,float el_eta){

    int bin=pH->FindBin(el_eta,el_pt);

    double epsilon=pH->GetBinContent(bin);

    double up=epsilon+pH->GetBinError(bin);
    double down=epsilon-pH->GetBinError(bin);

    if (down<0) down=0;
    if (up>1) up=1;

    std::vector<double> eps={epsilon,up,down};

    return eps;
}

qWeight* qValidate::ComputeWeight(float el_pt1,float el_pt2,float el_eta1,float el_eta2,std::vector<qHisto*> local_rates,float var=-1){

    double weight=-1;

    TH2F *pH=0;

    if (local_rates.size()==1) pH=local_rates[0]->GetHisto();
    else pH=ChooseHisto(local_rates,var);
  
    std::vector<double> epsilon_1=GetEpsilon(pH,el_pt1,el_eta1);
    std::vector<double> epsilon_2=GetEpsilon(pH,el_pt2,el_eta2);

    weight=GetWeight(epsilon_1[0],epsilon_2[0]);

    double up=GetWeight(epsilon_1[1],epsilon_2[1]);
    double down=GetWeight(epsilon_1[2],epsilon_2[2]);
      
    qWeight *total=new qWeight(weight,up,down);

    return total;
}

qWeight* qValidate::GetRateReweighting(float el_pt1,float el_pt2,float el_eta1,float el_eta2,TString dependence,float var=-1){
   
    double nom=0;
    double up2=0;
    double down2=0;

    for (std::map<TString,std::vector<qHisto*> >::iterator it=m_method_A.begin(); it!=m_method_A.end(); ++it){

         qWeight *qW=ComputeWeight(el_pt1,el_pt2,el_eta1,el_eta2,it->second,var);         

         float f=m_bkg_comp[it->first]; 
         float local_nom=qW->GetNominal();
         float sigma_up=qW->GetUp()-local_nom;
         float sigma_down=local_nom-qW->GetDown();

         nom=nom+f*qW->GetNominal();
         up2=up2+pow(f*sigma_up,2);
         down2=down2+pow(sigma_down,2);
    }
   
    qWeight *weight=new qWeight(nom,nom+TMath::Sqrt(up2),nom-TMath::Sqrt(down2));

    return weight;
}

double qValidate::GetWeight(double e1,double e2){

    double num=e1+e2-2*e1*e2;
    double den=1-e1-e2+2*e1*e2;

    double weight=num/den;

    return weight;
}

void qValidate::InitializeHistos(){

    std::cout << "-- Initliazing histos..." << std::endl;

    MiniTreeAnalyzer analyzer;  
    analyzer.AddWeight("weight_normalise"); //dummy just to avoid that the y-label show "Unweighted events"
 
    for (unsigned int i=0; i< v_var.size(); i++){

       TH1F *pH_SS=analyzer.CreateHisto(v_var[i],"SS");

       m_histo_SS[v_var[i]->GetTitle()+"_SS"]=pH_SS;

            for (std::map<TString,std::vector<qHisto*> >::iterator it=m_dependence.begin(); it!=m_dependence.end(); ++it){

                   TH1F *pH_OS=analyzer.CreateHisto(v_var[i],"OS_"+it->first);
                   m_histo_OS[v_var[i]->GetTitle()+"_OS_"+it->first]=pH_OS;

                   TH1F *pH_OSUp=analyzer.CreateHisto(v_var[i],"OSUp_"+it->first);
                   m_histo_OSUp[v_var[i]->GetTitle()+"_OSUp_"+it->first]=pH_OSUp;

                   TH1F *pH_OSDown=analyzer.CreateHisto(v_var[i],"OSDown_"+it->first);
                   m_histo_OSDown[v_var[i]->GetTitle()+"_OSDown_"+it->first]=pH_OSDown;


            }//end depen
    }//end var

    return;
}

std::vector<TString> qValidate::FillVector(TString file){

    std::vector<TString> samples;

    std::ifstream infile(file.Data());

    if (!infile){
        std::cout <<"-- File with the list of samples: '"<< file <<"' does not exist --- ABORTING " << std::endl;
        exit(1);
    }

    std::string line;

    while (std::getline(infile,line) ){

        if (line.empty()) continue;

        TString ss = line;

        samples.push_back(ss);
    }

    return samples;
}

void qValidate::AddDependence(TString name,TString latex,std::vector<TString> rates){

    std::vector<qHisto*> histos2D=GetHistos2D(rates);

    m_dependence[name]=histos2D;

    latex_dependence[name]=latex;

    return;
}

void qValidate::AddDependence(TString name,TString latex,std::map<TString,std::vector<TString> > rates){

    std::vector<qHisto*> histos2D=GetHistos2D(rates.begin()->second);

    for (std::map<TString,std::vector<TString> >::iterator it=rates.begin(); it!=rates.end(); ++it){ 
    
        m_method_A[it->first]=GetHistos2D(it->second);
    }

    m_dependence[name]=histos2D; //dummy: it will not be used->it has to be filled by consistency

    latex_dependence[name]=latex;

    return;
}

std::vector<qHisto*> qValidate::GetHistos2D(std::vector<TString> rates){

    std::vector<qHisto*> local;

    for (unsigned int i=0; i<rates.size(); i++){

        TString file=File_wRates+"/"+rates[i]+".root";
        TString histoname=rates[i]+"_all_num";

        std::cout << "-- Reading file with rates: " << file << std::endl;

        TFile *pFile = new TFile(file);

        TH2F *pH=(TH2F*)pFile->Get(histoname);
        pH->SetDirectory(0);

        qHisto *local_histo=new qHisto(pH);

        local.push_back(local_histo);   

        pFile->Close();
    }

    if (local.size()==0){

        std::cout << "-- Could not read file with rates from : " << File_wRates  << std::endl;
        exit(1);
    }

    return local;
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
        
        /*  TString str=tokens[i];
        if (str.Contains("ge")) str_min=str;
        else if (str.Contains("le")) str_max=str;      
        */
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

qWeight::qWeight(float w,float up=0,float d=0):
 w_nom(w),
 w_up(up),
 w_down(d)
{}
qWeight::~qWeight()
{}

