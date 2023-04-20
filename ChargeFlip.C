#include "TString.h"
#include <vector>
#include "TH1F.h"
#include "TH2F.h"

class qMisID{

   public:
      qMisID(TString name,TString latex,TString file,TString tree,TString extra="");
      ~qMisID();

      inline void PrintLog(bool val){printLog=val;};

      void PrintSummary();

      void InitializeSelections();
      void InitializeHistos();
      void Execute();
      TH2F *GetRatio(TH2F *num,TH2F* den); 

      inline TString GetName(){return m_name;};
      inline TString GetLatex(){return m_latex;};

      inline void SetPath(TString path){m_path=path;};
      inline void SetBaselineSelection(TString base){m_baseline_sel="("+base+")";};     
      inline void SetExtraSelection(TString extra){m_extra_sel=extra;};
      inline void SetSelectionLepton1(TString wrong_track,TString photon_conversion){m_sel_wrongTrack_1="("+wrong_track+")"; m_sel_photonConversion_1="("+photon_conversion+")";};
      inline void SetSelectionLepton2(TString wrong_track,TString photon_conversion){m_sel_wrongTrack_2="("+wrong_track+")"; m_sel_photonConversion_2="("+photon_conversion+")";};
      inline void SetWeight(TString weight){m_weight=weight;};

      inline TString GetFullSelectionWrongTrack1(){return m_full_sel_wrongTrack_1;};
      inline TString GetFullSelectionWrongTrack2(){return m_full_sel_wrongTrack_2;};     
      inline TString GetFullSelectionPhotonConversion1(){return m_full_sel_photonConversion_1;};
      inline TString GetFullSelectionPhotonConversion2(){return m_full_sel_photonConversion_2;};
      inline TString GetTotalSelectionLep1(){return m_full_sel_1;};
      inline TString GetTotalSelectionLep2(){return m_full_sel_2;};

      inline TString GetFullSelectionWrongTrack(){return m_full_sel_wrongTrack_all;};
      inline TString GetFullSelectionPhotonConversion(){return m_full_sel_photonConversion_all;};
      inline TString GetTotalSelectionLep(){return m_full_sel_all;}
 
      inline TString GetTotalBaselineSelection(){return m_full_baseline_all;};


      inline void SetVarEta(TString var1,TString var2,Float_t *binning,int n){m_var_eta1=var1; m_var_eta2=var2; m_eta_binning=binning; m_nbins_eta=n;};
      inline void SetVarPt(TString var1,TString var2,Float_t *binning,int n){m_var_pt1=var1; m_var_pt2=var2;m_pt_binning=binning; m_nbins_pt=n;};
             
      std::vector<TString> FillVector(TString list);
   
      std::vector<TH2F*> pH_Ratios;
  
      void ComputeRatios();
      void GetPlots();

   private:

      TH2F* CreateEtaPtHisto(TString label);
      void ReadInputFile(TString file);
      bool printLog;

      TString m_path;
      TString m_name;
      TString m_latex; 
      TString m_file;
      TString m_tree;
      TString m_baseline_sel;
      TString m_extra_sel;
      TString m_sel_wrongTrack_1;
      TString m_sel_photonConversion_1;
      TString m_sel_wrongTrack_2;
      TString m_sel_photonConversion_2;
      TString m_var_eta1;
      TString m_var_pt1;
      TString m_var_eta2;
      TString m_var_pt2;

      TString m_weight;
      TString m_full_sel_wrongTrack_1;
      TString m_full_sel_wrongTrack_2; 
      TString m_full_sel_photonConversion_1;
      TString m_full_sel_photonConversion_2;
      TString m_full_sel_1;
      TString m_full_sel_2;
      TString m_full_sel_wrongTrack_all;
      TString m_full_sel_photonConversion_all;
      TString m_full_sel_all;
      TString m_full_baseline_all;

      Float_t *m_eta_binning;
      Float_t *m_pt_binning;
      int m_nbins_eta;
      int m_nbins_pt;

      TH2F *m_histo_all_den;
      TH2F *m_histo_all_den1;
      TH2F *m_histo_all_den2;

      TH2F *m_histo_wrongTrack_num;
      TH2F *m_histo_photonConversion_num;
      TH2F *m_histo_all_num;

      TH2F *m_histo_wrongTrack_num1;
      TH2F *m_histo_wrongTrack_num2;
      TH2F *m_histo_photonConversion_num1;
      TH2F *m_histo_photonConversion_num2;     
      TH2F *m_histo_total_num1;
      TH2F *m_histo_total_num2; 

      TH2F *m_histo_ratio_wrongTrack_1;
      TH2F *m_histo_ratio_wrongTrack_2;
      TH2F *m_histo_ratio_photonConversion_1;
      TH2F *m_histo_ratio_photonConversion_2;
      TH2F *m_histo_ratio_total_1;
      TH2F *m_histo_ratio_total_2;

      TH2F *m_histo_ratio_wrongTrack;
      TH2F *m_histo_ratio_photonConversion;
      TH2F *m_histo_all_ratio;

};

qMisID MeasureChargeFlip(TString name,TString latex,TString file,TString extra_cut="");

void ChargeFlipProcess(TString process,TString process_latex,TString tag);

void ChargeFlip(){

   //TString tag="_NewBinning";
   TString tag="_Default";


   ChargeFlipProcess("ttbar","t#bar{t}",tag);

   ChargeFlipProcess("Zjets","Z+jets",tag);


 return;
}
void ChargeFlipProcess(TString process,TString process_latex,TString tag){
 
   qMisID process_ge0j= MeasureChargeFlip(process+tag+"_nJets_ge0",process_latex,process+".list");

 
   //#########  priority Ht ####### 
 /* qMisID process_ht_ge0_le200=MeasureChargeFlip(process+tag+"_Ht_ge0_le200",process_latex,process+".list","HT_all>0 && HT_all<=200000");
  qMisID process_ht_ge200_le400=MeasureChargeFlip(process+tag+"_Ht_ge200_le400",process_latex,process+".list","HT_all>200000 && HT_all<=400000");
  qMisID process_ht_ge400_le600=MeasureChargeFlip(process+tag+"_Ht_ge400_le600",process_latex,process+".list","HT_all>400000 && HT_all<=600000");
  qMisID process_ht_ge600_le800=MeasureChargeFlip(process+tag+"_Ht_ge600_le800",process_latex,process+".list","HT_all>600000 && HT_all<=800000");
  qMisID process_ht_ge800_le1000=MeasureChargeFlip(process+tag+"_Ht_ge800_le1000",process_latex,process+".list","HT_all>800000 && HT_all<=1000000");

  qMisID process_ht_ge1000=MeasureChargeFlip(process+tag+"_Ht_ge1000",process_latex,process+".list","HT_all>1000000");

  qMisID process_ht_ge0_le500=MeasureChargeFlip(process+tag+"_Ht_ge0_le500",process_latex,process+".list","HT_all>0 && HT_all<=500000");
  qMisID process_ht_ge500_le1000=MeasureChargeFlip(process+tag+"_Ht_ge500_le1000",process_latex,process+".list","HT_all>500000 && HT_all<=1000000");


  qMisID process_ge0j= MeasureChargeFlip(process+tag+"_nJets_ge0",process_latex,process+".list");
  qMisID process_le4j= MeasureChargeFlip(process+tag+"_nJets_le4",process_latex,process+".list","nJets<=4");
  qMisID process_ge5j= MeasureChargeFlip(process+tag+"_nJets_ge5",process_latex,process+".list","nJets>=5");

  qMisID process_le2bj= MeasureChargeFlip(process+tag+"_BJets_le2",process_latex,process+".list","nBTags_MV2c10_77<=2");
  qMisID process_ge3bj= MeasureChargeFlip(process+tag+"_BJets_ge3",process_latex,process+".list","nBTags_MV2c10_77>=3");
  */
  //#############################


   // qMisID process_met_ge80=MeasureChargeFlip(process+tag+"_met_ge80",process_latex,process+".list","met_met>=80000");
   // qMisID process_met_ge0_le80=MeasureChargeFlip(process+tag+"_met_ge0_le80",process_latex,process+".list","met_met<=80000");


  // Rates related to Etmiss
  /*qMisID process_met_ge0_le40=MeasureChargeFlip(process+tag+"_met_ge0_le40",process_latex,process+".list","met_met<=40000");
  qMisID process_met_ge40_le80=MeasureChargeFlip(process+tag+"_met_ge40_le80",process_latex,process+".list","met_met>=40000 && met_met<=80000");
  qMisID process_met_ge80_le120=MeasureChargeFlip(process+tag+"_met_ge80_le120",process_latex,process+".list","met_met>=80000 && met_met<=120000");
  qMisID process_met_ge120_le160=MeasureChargeFlip(process+tag+"_met_ge120_le160",process_latex,process+".list","met_met>=120000 && met_met<=160000");
  qMisID process_met_ge160_le200=MeasureChargeFlip(process+tag+"_met_ge160_le200",process_latex,process+".list","met_met>=160000 && met_met<=200000");
  qMisID process_met_ge200=MeasureChargeFlip(process+tag+"_met_ge200",process_latex,process+".list","met_met>=200000");

  qMisID process_met_ge40=MeasureChargeFlip(process+tag+"_met_ge40",process_latex,process+".list","met_met>=40000");
  qMisID process_met_ge80=MeasureChargeFlip(process+tag+"_met_ge80",process_latex,process+".list","met_met>=80000");
  qMisID process_met_ge120=MeasureChargeFlip(process+tag+"_met_ge120",process_latex,process+".list","met_met>=120000");
  qMisID process_met_ge160=MeasureChargeFlip(process+tag+"_met_ge160",process_latex,process+".list","met_met>=160000");
  qMisID process_met_ge240=MeasureChargeFlip(process+tag+"_met_ge240",process_latex,process+".list","met_met>=240000");
  qMisID process_met_ge280=MeasureChargeFlip(process+tag+"_met_ge280",process_latex,process+".list","met_met>=280000");
  */

  

  // Rates for primary vertex
  /*qMisID process_pv_ge10=MeasureChargeFlip(process+tag+"_pv_ge10",process_latex,process+".list","nPrimaryVtx>=10");
  qMisID process_pv_ge20=MeasureChargeFlip(process+tag+"_pv_ge20",process_latex,process+".list","nPrimaryVtx>=20");
  qMisID process_pv_ge30=MeasureChargeFlip(process+tag+"_pv_ge30",process_latex,process+".list","nPrimaryVtx>=30");

  qMisID process_pv_ge0_le10=MeasureChargeFlip(process+tag+"_pv_ge0_le10",process_latex,process+".list","nPrimaryVtx>=0 && nPrimaryVtx<=10");
  qMisID process_pv_ge10_le20=MeasureChargeFlip(process+tag+"_pv_ge10_le20",process_latex,process+".list","nPrimaryVtx>=10 && nPrimaryVtx<=20");
  qMisID process_pv_ge20_le30=MeasureChargeFlip(process+tag+"_pv_ge20_le30",process_latex,process+".list","nPrimaryVtx>=20 && nPrimaryVtx<=30");

  qMisID process_pv_ge0_le20=MeasureChargeFlip(process+tag+"_pv_ge0_le20",process_latex,process+".list","nPrimaryVtx>=0 && nPrimaryVtx<=20");
  */
 
  // Rates for HT_jets
  /*qMisID process_htjets_ge0_le100=MeasureChargeFlip(process+tag+"_Htjets_ge0_le100",process_latex,process+".list","HT_jets>0 && HT_jets<=100000");
  qMisID process_htjets_ge100_le200=MeasureChargeFlip(process+tag+"_Htjets_ge100_le200",process_latex,process+".list","HT_jets>100000 && HT_jets<=200000");
  qMisID process_htjets_ge200_le300=MeasureChargeFlip(process+tag+"_Htjets_ge200_le300",process_latex,process+".list","HT_jets>200000 && HT_jets<=300000");
  qMisID process_htjets_ge300_le400=MeasureChargeFlip(process+tag+"_Htjets_ge300_le400",process_latex,process+".list","HT_jets>300000 && HT_jets<=400000");
  qMisID process_htjets_ge400_le500=MeasureChargeFlip(process+tag+"_Htjets_ge400_le500",process_latex,process+".list","HT_jets>400000 && HT_jets<=500000");
  qMisID process_htjets_ge500_le600=MeasureChargeFlip(process+tag+"_Htjets_ge500_le600",process_latex,process+".list","HT_jets>500000 && HT_jets<=600000");
  qMisID process_htjets_ge600_le700=MeasureChargeFlip(process+tag+"_Htjets_ge600_le700",process_latex,process+".list","HT_jets>600000 && HT_jets<=700000");
  qMisID process_htjets_ge700_le800=MeasureChargeFlip(process+tag+"_Htjets_ge700_le800",process_latex,process+".list","HT_jets>700000 && HT_jets<=800000");

  qMisID process_htjets_ge0_le300=MeasureChargeFlip(process+tag+"_Htjets_ge0_le300",process_latex,process+".list","HT_jets>0 && HT_jets<=300000");
  qMisID process_htjets_ge300_le600=MeasureChargeFlip(process+tag+"_Htjets_ge300_le600",process_latex,process+".list","HT_jets>300000 && HT_jets<=600000");

  qMisID process_htjets_ge0_le400=MeasureChargeFlip(process+tag+"_Htjets_ge0_le400",process_latex,process+".list","HT_jets>0 && HT_jets<=400000");
  qMisID process_htjets_ge400_le800=MeasureChargeFlip(process+tag+"_Htjets_ge400_le800",process_latex,process+".list","HT_jets>400000 && HT_jets<=800000");

  qMisID process_htjets_ge100=MeasureChargeFlip(process+tag+"_Htjets_ge100",process_latex,process+".list","HT_jets>100000");
  qMisID process_htjets_ge200=MeasureChargeFlip(process+tag+"_Htjets_ge200",process_latex,process+".list","HT_jets>200000");
  qMisID process_htjets_ge300=MeasureChargeFlip(process+tag+"_Htjets_ge300",process_latex,process+".list","HT_jets>300000");
  qMisID process_htjets_ge400=MeasureChargeFlip(process+tag+"_Htjets_ge400",process_latex,process+".list","HT_jets>400000");
  qMisID process_htjets_ge500=MeasureChargeFlip(process+tag+"_Htjets_ge500",process_latex,process+".list","HT_jets>500000");
  qMisID process_htjets_ge600=MeasureChargeFlip(process+tag+"_Htjets_ge600",process_latex,process+".list","HT_jets>600000");
  qMisID process_htjets_ge700=MeasureChargeFlip(process+tag+"_Htjets_ge700",process_latex,process+".list","HT_jets>700000");
  qMisID process_htjets_ge800=MeasureChargeFlip(process+tag+"_Htjets_ge800",process_latex,process+".list","HT_jets>800000");
 */
  

 
  //#-----Rates related to Ht

 /* qMisID process_ht_ge0_le100=MeasureChargeFlip(process+tag+"_Ht_ge0_le100",process_latex,process+".list","HT_all>0 && HT_all<=100000");
  qMisID process_ht_ge100_le200=MeasureChargeFlip(process+tag+"_Ht_ge100_le200",process_latex,process+".list","HT_all>100000 && HT_all<=200000");
  qMisID process_ht_ge200_le300=MeasureChargeFlip(process+tag+"_Ht_ge200_le300",process_latex,process+".list","HT_all>200000 && HT_all<=300000");
  qMisID process_ht_ge300_le400=MeasureChargeFlip(process+tag+"_Ht_ge300_le400",process_latex,process+".list","HT_all>300000 && HT_all<=400000");
  qMisID process_ht_ge400_le500=MeasureChargeFlip(process+tag+"_Ht_ge400_le500",process_latex,process+".list","HT_all>400000 && HT_all<=500000");
  qMisID process_ht_ge500_le600=MeasureChargeFlip(process+tag+"_Ht_ge500_le600",process_latex,process+".list","HT_all>500000 && HT_all<=600000");
  qMisID process_ht_ge600_le700=MeasureChargeFlip(process+tag+"_Ht_ge600_le700",process_latex,process+".list","HT_all>600000 && HT_all<=700000");
  qMisID process_ht_ge700_le800=MeasureChargeFlip(process+tag+"_Ht_ge700_le800",process_latex,process+".list","HT_all>700000 && HT_all<=800000");
  qMisID process_ht_ge800_le900=MeasureChargeFlip(process+tag+"_Ht_ge800_le900",process_latex,process+".list","HT_all>800000 && HT_all<=900000");
  qMisID process_ht_ge900_le1000=MeasureChargeFlip(process+tag+"_Ht_ge900_le1000",process_latex,process+".list","HT_all>900000 && HT_all<=1000000");
  */

 /* qMisID process_ht_ge0_le300=MeasureChargeFlip(process+tag+"_Ht_ge0_le300",process_latex,process+".list","HT_all>0 && HT_all<=300000");
  qMisID process_ht_ge300_le600=MeasureChargeFlip(process+tag+"_Ht_ge300_le600",process_latex,process+".list","HT_all>300000 && HT_all<=600000");    
  qMisID process_ht_ge600_le900=MeasureChargeFlip(process+tag+"_Ht_ge600_le900",process_latex,process+".list","HT_all>600000 && HT_all<=900000");
 
  qMisID process_ht_ge0_le400=MeasureChargeFlip(process+tag+"_Ht_ge0_le400",process_latex,process+".list","HT_all>0 && HT_all<=400000");
  qMisID process_ht_ge400_le800=MeasureChargeFlip(process+tag+"_Ht_ge400_le800",process_latex,process+".list","HT_all>400000 && HT_all<=800000");    
  
  qMisID process_ht_ge100=MeasureChargeFlip(process+tag+"_Ht_ge100",process_latex,process+".list","HT_all>100000");
  qMisID process_ht_ge200=MeasureChargeFlip(process+tag+"_Ht_ge200",process_latex,process+".list","HT_all>200000");
  qMisID process_ht_ge300=MeasureChargeFlip(process+tag+"_Ht_ge300",process_latex,process+".list","HT_all>300000");
  qMisID process_ht_ge400=MeasureChargeFlip(process+tag+"_Ht_ge400",process_latex,process+".list","HT_all>400000");
  qMisID process_ht_ge500=MeasureChargeFlip(process+tag+"_Ht_ge500",process_latex,process+".list","HT_all>500000");
  qMisID process_ht_ge600=MeasureChargeFlip(process+tag+"_Ht_ge600",process_latex,process+".list","HT_all>600000");
  qMisID process_ht_ge700=MeasureChargeFlip(process+tag+"_Ht_ge700",process_latex,process+".list","HT_all>700000");
  qMisID process_ht_ge800=MeasureChargeFlip(process+tag+"_Ht_ge800",process_latex,process+".list","HT_all>800000");
  qMisID process_ht_ge900=MeasureChargeFlip(process+tag+"_Ht_ge900",process_latex,process+".list","HT_all>900000");
  */ 


 //##--- Rates related to nJets
 
 //qMisID process_ge0j= MeasureChargeFlip(process+tag+"_nJets_ge0",process_latex,process+".list");
 /*qMisID process_ge1j= MeasureChargeFlip(process+tag+"_nJets_ge1",process_latex,process+".list","nJets>=1");
 qMisID process_ge2j= MeasureChargeFlip(process+tag+"_nJets_ge2",process_latex,process+".list","nJets>=2");  
 qMisID process_ge3j= MeasureChargeFlip(process+tag+"_nJets_ge3",process_latex,process+".list","nJets>=3");
 qMisID process_ge4j= MeasureChargeFlip(process+tag+"_nJets_ge4",process_latex,process+".list","nJets>=4");
 //qMisID process_ge5j= MeasureChargeFlip(process+tag+"_nJets_ge5",process_latex,process+".list","nJets>=5"); //included in priority
 qMisID process_ge6j= MeasureChargeFlip(process+tag+"_nJets_ge6",process_latex,process+".list","nJets>=6");
 qMisID process_ge7j= MeasureChargeFlip(process+tag+"_nJets_ge7",process_latex,process+".list","nJets>=7");
 
 qMisID process_le0j= MeasureChargeFlip(process+tag+"_nJets_le0",process_latex,process+".list","nJets<=0");
 qMisID process_le1j= MeasureChargeFlip(process+tag+"_nJets_le1",process_latex,process+".list","nJets<=1");
 qMisID process_le2j= MeasureChargeFlip(process+tag+"_nJets_le2",process_latex,process+".list","nJets<=2");  
 qMisID process_le3j= MeasureChargeFlip(process+tag+"_nJets_le3",process_latex,process+".list","nJets<=3");    #### from HERE!!!!
 //qMisID process_le4j= MeasureChargeFlip(process+tag+"_nJets_le4",process_latex,process+".list","nJets<=4"); //included in priority
 qMisID process_le5j= MeasureChargeFlip(process+tag+"_nJets_le5",process_latex,process+".list","nJets<=5");
 qMisID process_le6j= MeasureChargeFlip(process+tag+"_nJets_le6",process_latex,process+".list","nJets<=6");
 qMisID process_le7j= MeasureChargeFlip(process+tag+"_nJets_le7",process_latex,process+".list","nJets<=7");

 qMisID process_ge2_le4= MeasureChargeFlip(process+tag+"_nJets_ge2_le4",process_latex,process+".list","nJets >=2 && nJets<=4"); 
 qMisID process_ge5_le6= MeasureChargeFlip(process+tag+"_nJets_ge5_le6",process_latex,process+".list","nJets >=5 && nJets<=6");
 */

 // qMisID process_le1bj= MeasureChargeFlip(process+tag+"_BJets_le1",process_latex,process+".list","nBTags_MV2c10_77<=1");
 // qMisID process_ge2bj= MeasureChargeFlip(process+tag+"_BJets_ge2",process_latex,process+".list","nBTags_MV2c10_77>=2");  

  //----- Rates for b-jets 
  
/* qMisID process_ge1bj= MeasureChargeFlip(process+tag+"_BJets_ge1",process_latex,process+".list","nBTags_MV2c10_77>=1");
 qMisID process_ge2bj= MeasureChargeFlip(process+tag+"_BJets_ge2",process_latex,process+".list","nBTags_MV2c10_77>=2");  
 //qMisID process_ge3bj= MeasureChargeFlip(process+tag+"_BJets_ge3",process_latex,process+".list","nBTags_MV2c10_77>=3"); //included in priority
 qMisID process_ge4bj= MeasureChargeFlip(process+tag+"_BJets_ge4",process_latex,process+".list","nBTags_MV2c10_77>=4");

 qMisID process_le1bj= MeasureChargeFlip(process+tag+"_BJets_le1",process_latex,process+".list","nBTags_MV2c10_77<=1");
 //qMisID process_le2bj= MeasureChargeFlip(process+tag+"_BJets_le2",process_latex,process+".list","nBTags_MV2c10_77<=2"); //included in priority
 qMisID process_le3bj= MeasureChargeFlip(process+tag+"_BJets_le3",process_latex,process+".list","nBTags_MV2c10_77<=3");
 qMisID process_le4bj= MeasureChargeFlip(process+tag+"_BJets_le4",process_latex,process+".list","nBTags_MV2c10_77<=4");

 qMisID process_BJets_ge2_le3= MeasureChargeFlip(process+tag+"_BJets_ge2_le3",process_latex,process+".list","nBTags_MV2c10_77 >=2 && nBTags_MV2c10_77<=3");
 */


  //Ranges for mu
 /* qMisID process_mu_ge10=MeasureChargeFlip(process+tag+"_mu_ge10",process_latex,process+".list","mu>=10");
  qMisID process_mu_ge20=MeasureChargeFlip(process+tag+"_mu_ge20",process_latex,process+".list","mu>=20");
  qMisID process_mu_ge30=MeasureChargeFlip(process+tag+"_mu_ge30",process_latex,process+".list","mu>=30");
  qMisID process_mu_ge40=MeasureChargeFlip(process+tag+"_mu_ge40",process_latex,process+".list","mu>=40");
  qMisID process_mu_ge50=MeasureChargeFlip(process+tag+"_mu_ge50",process_latex,process+".list","mu>=50");
  qMisID process_mu_ge60=MeasureChargeFlip(process+tag+"_mu_ge60",process_latex,process+".list","mu>=60");
  qMisID process_mu_ge70=MeasureChargeFlip(process+tag+"_mu_ge70",process_latex,process+".list","mu>=70"); 
  
  qMisID process_mu_ge0_le30=MeasureChargeFlip(process+tag+"_mu_ge0_le30",process_latex,process+".list","mu>=0 && mu<=30");
  qMisID process_mu_ge30_le40=MeasureChargeFlip(process+tag+"_mu_ge30_le40",process_latex,process+".list","mu>=30 && mu<=40"); 
  qMisID process_mu_ge40_le55=MeasureChargeFlip(process+tag+"_mu_ge40_le55",process_latex,process+".list","mu>=40 && mu<=55");
  qMisID process_mu_ge55=MeasureChargeFlip(process+tag+"_mu_ge55",process_latex,process+".list","mu>=55");

  qMisID process_mu_ge0_le20=MeasureChargeFlip(process+tag+"_mu_ge0_le20",process_latex,process+".list","mu>=0 && mu<=20");
  qMisID process_mu_ge20_le40=MeasureChargeFlip(process+tag+"_mu_ge20_le40",process_latex,process+".list","mu>=20 && mu<=40"); 
  qMisID process_mu_ge40_le60=MeasureChargeFlip(process+tag+"_mu_ge40_le60",process_latex,process+".list","mu>=40 && mu<=60");
  */
 
  
  


  
 //#################### Old staff....
  /*
  qMisID process_0j= MeasureChargeFlip(process+"_0j",process_latex,process+".list","nJets==0");
  qMisID process_1j= MeasureChargeFlip(process+"_1j",process_latex,process+".list","nJets==1");
  qMisID process_2j= MeasureChargeFlip(process+"_2j",process_latex,process+".list","nJets==2");  
  qMisID process_3j= MeasureChargeFlip(process+"_3j",process_latex,process+".list","nJets==3");
  qMisID process_4j= MeasureChargeFlip(process+"_4j",process_latex,process+".list","nJets==4");
  qMisID process_5j= MeasureChargeFlip(process+"_5j",process_latex,process+".list","nJets==5");
  qMisID process_6j= MeasureChargeFlip(process+"_6j",process_latex,process+".list","nJets==6");
  qMisID process_7j= MeasureChargeFlip(process+"_7j",process_latex,process+".list","nJets==7");
 

  qMisID process_0bj= MeasureChargeFlip(process+"_0bj",process_latex,process+".list","nBTags_MV2c10_77==0");
  qMisID process_1bj= MeasureChargeFlip(process+"_1bj",process_latex,process+".list","nBTags_MV2c10_77==1");
  qMisID process_2bj= MeasureChargeFlip(process+"_2bj",process_latex,process+".list","nBTags_MV2c10_77==2");
  qMisID process_3bj= MeasureChargeFlip(process+"_3bj",process_latex,process+".list","nBTags_MV2c10_77==3");
  qMisID process_4bj= MeasureChargeFlip(process+"_4bj",process_latex,process+".list","nBTags_MV2c10_77==4");
   */  

 // std::cout<< "-- Name " << process.GetName() << std::endl;


 return;
}

qMisID MeasureChargeFlip(TString name,TString latex,TString file,TString extra_cut=""){
 //Preselection: nElectrons==2 && nMuons==0 && el_pt[0]>28000 && el_pt[1]>28000 && fabs(el_eta[0])<2.47 && fabs(el_eta[1])<2.47 && el_LHTight[0]>0 && el_LHTight[1]>0 && el_isoFixedCutTight[0]>0 && el_isoFixedCutTight[1]>0 && fabs(el_d0sig[0])<5 && fabs(el_d0sig[1])<5 && fabs(el_delta_z0_sintheta[0])<0.5 && fabs(el_delta_z0_sintheta[1])<0.5 && (fabs(el_eta[0])<1.37  || fabs(el_eta[0])>1.52) && (fabs(el_eta[1])<1.37  || fabs(el_eta[1])>1.52)
 TString ori="-1"; 

 if (name.Contains("Zjets")) ori="13";
 else if (name.Contains("ttbar")) ori="10";

 TString preselection="nElectrons==2 && ((el_true_type[0]==2 && el_true_origin[0]=="+ori+") || (el_true_type[0]==4 && el_true_origin[0]==5 && fabs(el_firstEgMotherPdgId[0])==11) ) && ((el_true_type[1]==2 && el_true_origin[1]=="+ori+") || (el_true_type[1]==4 && el_true_origin[1]==5 && fabs(el_firstEgMotherPdgId[1])==11) )";
 //TString extra_cut="nJets>=1";

 TString selection_wrongTrack_lep1="el_true_type[0]==2 && el_true_origin[0]=="+ori+" && el_charge[0]*el_true_pdg[0]>0";  
 TString selection_wrongTrack_lep2="el_true_type[1]==2 && el_true_origin[1]=="+ori+" && el_charge[1]*el_true_pdg[1]>0";
 TString selection_photonConversion_lep1="el_true_type[0]==4 && el_true_origin[0]==5 && el_charge[0]*el_firstEgMotherPdgId[0]>0"; 
 TString selection_photonConversion_lep2="el_true_type[1]==4 && el_true_origin[1]==5 && el_charge[1]*el_firstEgMotherPdgId[1]>0";

 TString TreeName="nominal_Loose";
 TString varEta0="fabs(el_eta[0])";
 TString varEta1="fabs(el_eta[1])";
 TString varPt0="el_pt[0]*0.001"; //converted to GeV
 TString varPt1="el_pt[1]*0.001"; //converted to GeV
 TString path="/eos/user/d/dparedes/mc16d/ge0j/Skimming/";
 TString weight="weight_bTagSF_MV2c10_77*weight_mc*weight_leptonSF*weight_pileup*weight_jvt*weight_normalise*46.9";

 Float_t EtaBinning[]={0,0.6,1.1,1.52,1.7,2.3,2.5};
 Float_t PtBinning_Temp[]={0,60,90,130,200,2500};

 int nPt=5;
 
 Float_t PtBinning_Default[]={0,60,90,130,2500};

 if (name.Contains("Default")) nPt=4;
 

 Float_t PtBinning[]={};

 if (name.Contains("Default")) std::copy(PtBinning_Default, PtBinning_Default+sizeof(PtBinning_Default)/sizeof(PtBinning_Default[0]), PtBinning);
 else std::copy(PtBinning_Temp,PtBinning_Temp+sizeof(PtBinning_Temp)/sizeof(PtBinning_Temp[0]), PtBinning);


 qMisID qRates(name,latex,file,TreeName);
 
 qRates.SetPath(path);

 qRates.SetBaselineSelection(preselection);
 qRates.SetSelectionLepton1(selection_wrongTrack_lep1,selection_photonConversion_lep1);
 qRates.SetSelectionLepton2(selection_wrongTrack_lep2,selection_photonConversion_lep2);
 qRates.SetExtraSelection(extra_cut);
 
 qRates.SetVarEta(varEta0,varEta1,EtaBinning,6);
 qRates.SetVarPt(varPt0,varPt1,PtBinning,nPt);
 qRates.PrintLog(true);
 qRates.SetWeight(weight);

 qRates.Execute();



 return qRates;  
}
qMisID::qMisID(TString name,TString latex,TString file,TString tree,TString extra=""):

 m_name(name),
 m_latex(latex),
 m_file(file),
 m_tree(tree),
 m_path(""),
 m_weight(""),
 m_baseline_sel(""),
 m_extra_sel(extra),
 m_sel_wrongTrack_1(""),
 m_sel_photonConversion_1(""),
 m_sel_wrongTrack_2(""),
 m_sel_photonConversion_2(""),
 m_full_sel_wrongTrack_1(""),
 m_full_sel_wrongTrack_2(""), 
 m_full_sel_photonConversion_1(""),
 m_full_sel_photonConversion_2(""),
 m_full_sel_1(""),
 m_full_sel_2(""),
 m_full_sel_wrongTrack_all(""),
 m_full_sel_photonConversion_all(""),
 m_full_sel_all(""),
 m_full_baseline_all(""),
 m_var_eta1(""),
 m_var_eta2(""),
 m_var_pt1(""),
 m_var_pt2(""),
 m_nbins_eta(0),
 m_nbins_pt(0),
 printLog(false),
 m_histo_all_den(0),
 m_histo_all_den1(0),
 m_histo_all_den2(0),
 m_histo_wrongTrack_num(0),
 m_histo_photonConversion_num(0),
 m_histo_all_num(0),
 m_histo_ratio_wrongTrack(0),
 m_histo_ratio_photonConversion(0),
 m_histo_all_ratio(0),
 m_histo_wrongTrack_num1(0),
 m_histo_wrongTrack_num2(0),
 m_histo_photonConversion_num1(0),
 m_histo_photonConversion_num2(0),
 m_histo_total_num1(0),
 m_histo_total_num2(0),
 m_histo_ratio_wrongTrack_1(0),
 m_histo_ratio_wrongTrack_2(0),
 m_histo_ratio_photonConversion_1(0),
 m_histo_ratio_photonConversion_2(0),
 m_histo_ratio_total_1(0),
 m_histo_ratio_total_2(0),
 pH_Ratios(0)
{}

qMisID::~qMisID()
{

 //delete m_histo_num;
 //delete m_histo_den;
 //delete m_histo_ratio;

};

void qMisID::ReadInputFile(TString file){

  TString local_file=m_path+file;
 
  TFile *pInput = new TFile(local_file);
  if (!pInput){
        printf("-- File %s is missing\n", pInput->GetName());
        exit(1);
  }
  else {
      if (printLog) cout << "-- Opening file :" << file <<  endl;
  }


 TTree* local_Tree = (TTree*)pInput->Get(m_tree);
    
 if (printLog) std::cout<< "-- Tree entries:" << local_Tree->GetEntries() << std::endl;


 //Creating local histos
 TH2F *local_wrongTrack_num=CreateEtaPtHisto("local_wrongTrack_num");
 TH2F *local_photonConversion_num=CreateEtaPtHisto("local_photonConversion_num");
 TH2F *local_all_num=CreateEtaPtHisto("local_all_num");

 TH2F *local_wrongTrack_num1=CreateEtaPtHisto("local_wrongTrack_num1");
 TH2F *local_wrongTrack_num2=CreateEtaPtHisto("local_wrongTrack_num2");
 TH2F *local_photonConversion_num1=CreateEtaPtHisto("local_photonConversion_num1");
 TH2F *local_photonConversion_num2=CreateEtaPtHisto("local_photonConversion_num2");
 TH2F *local_total_num1=CreateEtaPtHisto("local_total_num1");
 TH2F *local_total_num2=CreateEtaPtHisto("local_total_num2");

 TH2F *local_den1=CreateEtaPtHisto("local_all_den1");
 TH2F *local_den2=CreateEtaPtHisto("local_all_den2");

 TString wrongTrack1=GetFullSelectionWrongTrack1();
 TString wrongTrack2=GetFullSelectionWrongTrack2();
 
 TString photon1=GetFullSelectionPhotonConversion1();
 TString photon2=GetFullSelectionPhotonConversion2();

 TString total1=GetTotalSelectionLep1();
 TString total2=GetTotalSelectionLep2();

 TString base=GetTotalBaselineSelection();

 if (m_weight!=""){
    wrongTrack1=m_weight+"*("+wrongTrack1+")";
    wrongTrack2=m_weight+"*("+wrongTrack2+")";
    photon1=m_weight+"*("+photon1+")";
    photon2=m_weight+"*("+photon2+")";
    total1=m_weight+"*("+total1+")";
    total2=m_weight+"*("+total2+")";
    base=m_weight+"*("+base+")";
 }

 //Filling histos
 local_Tree->Draw(m_var_pt1+":"+m_var_eta1+" >> "+local_wrongTrack_num1->GetName(),wrongTrack1,"goff"); //Fill charge flipped electrons just for wrong track
 local_Tree->Draw(m_var_pt2+":"+m_var_eta2+" >> "+local_wrongTrack_num2->GetName(),wrongTrack2,"goff"); //Fill charge flipped electrons just for wrong track

 local_Tree->Draw(m_var_pt1+":"+m_var_eta1+" >> "+local_photonConversion_num1->GetName(),photon1,"goff"); //Fill charge flipped electrons just for photon conversion
 local_Tree->Draw(m_var_pt2+":"+m_var_eta2+" >> "+local_photonConversion_num2->GetName(),photon2,"goff"); //Fill charge flipped electrons just for photon conversion

 local_Tree->Draw(m_var_pt1+":"+m_var_eta1+" >> "+local_total_num1->GetName(),total1,"goff"); 
 local_Tree->Draw(m_var_pt2+":"+m_var_eta2+" >> "+local_total_num2->GetName(),total2,"goff");     

 local_wrongTrack_num->Add(local_wrongTrack_num1); //total wrong track 1+2
 local_wrongTrack_num->Add(local_wrongTrack_num2); //total wrong track

 local_photonConversion_num->Add(local_photonConversion_num1); //total photon conversion 1+2
 local_photonConversion_num->Add(local_photonConversion_num2);

 local_all_num->Add(local_wrongTrack_num);  //adding wrong track + photon conversion
 local_all_num->Add(local_photonConversion_num);


 local_Tree->Draw(m_var_pt1+":"+m_var_eta1+" >> "+local_den1->GetName(),base,"goff"); //baseline selection
 local_Tree->Draw(m_var_pt2+":"+m_var_eta2+" >> "+local_den2->GetName(),base,"goff");
 


 m_histo_wrongTrack_num1->Add(local_wrongTrack_num1);
 m_histo_wrongTrack_num2->Add(local_wrongTrack_num2);
 m_histo_photonConversion_num1->Add(local_photonConversion_num1); 
 m_histo_photonConversion_num2->Add(local_photonConversion_num2);

 m_histo_total_num1->Add(local_total_num1); 
 m_histo_total_num2->Add(local_total_num2);

 m_histo_wrongTrack_num->Add(local_wrongTrack_num);
 m_histo_photonConversion_num->Add(local_photonConversion_num);

 m_histo_all_num->Add(local_all_num);

 m_histo_all_den1->Add(local_den1);
 m_histo_all_den2->Add(local_den2);
 m_histo_all_den->Add(local_den1);
 m_histo_all_den->Add(local_den2);



 std::cout << "### DEBUG : track 1: " << local_wrongTrack_num1->Integral() << std::endl;
 std::cout << "### DEBUG : track 2: " << local_wrongTrack_num2->Integral() << std::endl;
 std::cout << "### DEBUG : photon 1: " << local_photonConversion_num1->Integral() << std::endl;
 std::cout << "### DEBUG : photon 2: " << local_photonConversion_num2->Integral() << std::endl;
 std::cout << "### DEBUG : total 1: " << local_total_num1->Integral() << std::endl;
 std::cout << "### DEBUG : total 2: " << local_total_num2->Integral() << std::endl;
 std::cout << "### DEBUG : total track: " << local_wrongTrack_num->Integral() << std::endl;
 std::cout << "### DEBUG : total photon: " << local_photonConversion_num->Integral() << std::endl;
 std::cout << "### DEBUG : total all num: " << local_all_num->Integral() << std::endl;
 std::cout << "### DEBUG : total all den1: " << local_den1->Integral() << std::endl;
 std::cout << "### DEBUG : total all den2: " << local_den2->Integral() << std::endl;
 std::cout << "### DEBUG : total all den: " << local_den1->Integral()+local_den2->Integral() << std::endl;

 //TCanvas *pC= new TCanvas;
 //pC->cd(1);
 //m_histo_num->Draw("COLZ,TEXT");


  return;
}




TH2F* qMisID::CreateEtaPtHisto(TString label){
 
  TH2F *hist=new TH2F(m_name+"_"+label,m_name+"_"+label,m_nbins_eta,m_eta_binning,m_nbins_pt,m_pt_binning);
  hist->GetYaxis()->SetTitle("Lepton p_{T} [GeV]");
  hist->GetXaxis()->SetTitle("Lepton #eta");
  hist->Reset();

  hist->Sumw2(); 

 return hist;
}

void qMisID::Execute(){


  std::vector<TString> samples=FillVector(m_file);


  InitializeHistos();  

  InitializeSelections();

  if (printLog) PrintSummary();
  for (unsigned int i=0; i<samples.size(); i++) ReadInputFile(samples[i]);

  ComputeRatios();
  
  GetPlots();   
  



 return;
}

void qMisID::GetPlots(){


 //  m_histo_ratio->Draw();

 //for (unsigned int i=0; i<m_nbins_pt; i++){ 


 //   TH1D *pH=m_histo_ratio->ProjectionX("bin_pt",1,)


// }


 return;
}
TH2F *qMisID::GetRatio(TH2F *num,TH2F* den){

 TH2F *pH=(TH2F*)num->Clone();

 pH->Divide(den);

 return pH;
}
void qMisID::ComputeRatios(){

  m_histo_ratio_wrongTrack_1=GetRatio(m_histo_wrongTrack_num1,m_histo_all_den1);
  m_histo_ratio_wrongTrack_2=GetRatio(m_histo_wrongTrack_num2,m_histo_all_den2);
 
  m_histo_ratio_photonConversion_1=GetRatio(m_histo_photonConversion_num1,m_histo_all_den1);
  m_histo_ratio_photonConversion_2=GetRatio(m_histo_photonConversion_num2,m_histo_all_den2);

  m_histo_ratio_total_1=GetRatio(m_histo_total_num1,m_histo_all_den1);
  m_histo_ratio_total_2=GetRatio(m_histo_total_num2,m_histo_all_den2);

  m_histo_ratio_wrongTrack=GetRatio(m_histo_wrongTrack_num,m_histo_all_den);
  m_histo_ratio_photonConversion=GetRatio(m_histo_photonConversion_num,m_histo_all_den); 

  m_histo_all_ratio=GetRatio(m_histo_all_num,m_histo_all_den);


  pH_Ratios.push_back(m_histo_ratio_wrongTrack_1);
  pH_Ratios.push_back(m_histo_ratio_wrongTrack_2);
  pH_Ratios.push_back(m_histo_ratio_photonConversion_1);
  pH_Ratios.push_back(m_histo_ratio_photonConversion_2);
  pH_Ratios.push_back(m_histo_ratio_total_1);
  pH_Ratios.push_back(m_histo_ratio_total_2);
  pH_Ratios.push_back(m_histo_ratio_wrongTrack);
  pH_Ratios.push_back(m_histo_ratio_photonConversion);
  pH_Ratios.push_back(m_histo_all_ratio);

  TString fileOutput="Output/"+GetName()+".root";


  TFile *pOutput = new TFile(fileOutput,"RECREATE");


  for (unsigned int i=0; i<pH_Ratios.size(); i++){

    TCanvas *pC = new TCanvas;
    pC->cd();
    pH_Ratios[i]->Draw("COLZ,TEXT");
    pH_Ratios[i]->Write();

  }

  pOutput->Close();


 /* 
  m_histo_all_ratio;
*/



 /*
 
 TCanvas *pC = new TCanvas;
 pC->Divide(2,2,0,0);
 pC->cd(1);
 m_histo_num->Draw("COLZ,TEXT");

 pC->cd(2);
 m_histo_den->Draw("COLZ,TEXT");

 pC->cd(3);

 m_histo_ratio->Draw("COLZ,TEXT");

 std::cout << "-- Integral rates: " << m_histo_ratio->Integral() << std::endl;

 */
 return;
}
void qMisID::InitializeHistos(){

  m_histo_all_den=CreateEtaPtHisto("all_den"); 
  m_histo_all_den1=CreateEtaPtHisto("all_den1");
  m_histo_all_den2=CreateEtaPtHisto("all_den2");

  m_histo_wrongTrack_num=CreateEtaPtHisto("wrongTrack_num");
  m_histo_photonConversion_num=CreateEtaPtHisto("photon_conversion_num");
  m_histo_all_num=CreateEtaPtHisto("all_num");

  m_histo_ratio_wrongTrack=CreateEtaPtHisto("ratio_wrongTrack");
  m_histo_ratio_photonConversion=CreateEtaPtHisto("ratio_photonConversion");
  m_histo_all_ratio=CreateEtaPtHisto("all_ratio");

  m_histo_wrongTrack_num1=CreateEtaPtHisto("wrongTrack_num1");
  m_histo_wrongTrack_num2=CreateEtaPtHisto("wrongTrack_num2");
  m_histo_photonConversion_num1=CreateEtaPtHisto("photonConversion_num1");
  m_histo_photonConversion_num2=CreateEtaPtHisto("photonConversion_num2");
  m_histo_total_num1=CreateEtaPtHisto("total_num1");
  m_histo_total_num2=CreateEtaPtHisto("total_num2");

  m_histo_ratio_wrongTrack_1=CreateEtaPtHisto("ratio_wrongTrack_1");
  m_histo_ratio_wrongTrack_2=CreateEtaPtHisto("ratio_wrongTrack_2");
  m_histo_ratio_photonConversion_1=CreateEtaPtHisto("ratio_photonConversion_1");
  m_histo_ratio_photonConversion_2=CreateEtaPtHisto("ratio_photonConversion_2");
  m_histo_ratio_total_1=CreateEtaPtHisto("ratio_total_1");
  m_histo_ratio_total_2=CreateEtaPtHisto("tatio_total_2");

 return;
}
void qMisID::InitializeSelections(){


 m_full_baseline_all=m_baseline_sel;
 if (m_extra_sel != "") m_full_baseline_all="("+m_baseline_sel+" && ("+m_extra_sel+"))";

     
 m_full_sel_wrongTrack_1=m_full_baseline_all+" && "+m_sel_wrongTrack_1;
 m_full_sel_wrongTrack_2=m_full_baseline_all+" && "+m_sel_wrongTrack_2;

 m_full_sel_photonConversion_1=m_full_baseline_all+" && "+m_sel_photonConversion_1;
 m_full_sel_photonConversion_2=m_full_baseline_all+" && "+m_sel_photonConversion_2;

 m_full_sel_1=m_full_baseline_all+" && ("+m_sel_wrongTrack_1+" || "+m_sel_photonConversion_1+")";
 m_full_sel_2=m_full_baseline_all+" && ("+m_sel_wrongTrack_2+" || "+m_sel_photonConversion_2+")";

 //m_full_sel_wrongTrack_all=m_full_baseline_all+" && ("+m_sel_wrongTrack_1+" || "+m_sel_wrongTrack_2+")";
 //m_full_sel_photonConversion_all=m_full_baseline_all+" && ("+m_sel_photonConversion_1+" || "+m_sel_photonConversion_2+")";

 //m_full_sel_all=m_full_baseline_all+" && ("+m_sel_photonConversion_1+" || "+m_sel_photonConversion_2+" ||"+m_sel_wrongTrack_1+" || "+m_sel_wrongTrack_2+")";


  return;
}
void qMisID::PrintSummary(){

 std::cout << "### Printing settings: " << std::endl;
 std::cout << "- Preselection: " << GetTotalBaselineSelection() << std::endl;
 std::cout << "- Extra cut: " << m_extra_sel << std::endl;
 std::cout << "- Sel. lep 1: " << GetTotalSelectionLep1() << std::endl;
 std::cout << "- Sel. lep 2: " << GetTotalSelectionLep2() << std::endl;



 return;
}
std::vector<TString> qMisID::FillVector(TString file){

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
