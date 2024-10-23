#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

void Mass_in_bdt()
{
   Float_t pt0,pt1,pt2,pt3,pt4,pt5;
   Float_t ptB1,ptB2,etaB1,etaB2,phiB1,phiB2,BtgB1_DeepCSV,BtgB1_DeepJet,BtgB2_DeepCSV,BtgB2_DeepJet,OrdB1,OrdB2;
   Float_t ptQ1,ptQ2,etaQ1,etaQ2,phiQ1,phiQ2,BtgQ1_DeepCSV,BtgQ1_DeepJet,BtgQ2_DeepCSV,BtgQ2_DeepJet,QtgQ1,QtgQ2;
   Float_t MinvQQ,dETAQQ,dPHIQQ,NJ_InGap,ptAll,pzAll,dETA_bb_qq,alphabb,alphaqq;
   Float_t NJ_20,E_rest_20,HTT_rest_20,NJ_30,E_rest_30,HTT_rest_30;
   Float_t zep,phiA_bb_qq,minR,dR_leadqH,dR_subleadqH,norm_pt,prod_pt,ptHtoM;
   Float_t MinvBB;
   Float_t MinvBB_Reg;

   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

   reader->AddVariable( "T_mqq",                    &MinvQQ);
   reader->AddVariable( "T_dETAqq",                 &dETAQQ);
   reader->AddVariable( "T_dPHIqq",                 &dPHIQQ);
   reader->AddVariable( "T_btgb1_DeepCSV",          &BtgB1_DeepCSV);
   reader->AddVariable( "T_btgb1_DeepJet",          &BtgB1_DeepJet);
   reader->AddVariable( "T_btgb2_DeepCSV",          &BtgB2_DeepCSV);
   reader->AddVariable( "T_btgb2_DeepJet",          &BtgB2_DeepJet);
   reader->AddVariable( "T_qtgq1",                  &QtgQ1);
   reader->AddVariable( "T_qtgq2",                  &QtgQ2);
   reader->AddVariable( "T_NJ_ingap",               &NJ_InGap);
//   reader->AddVariable( "T_ptAll",                  &ptAll);
   reader->AddVariable( "T_pzAll",                  &pzAll);

   reader->AddVariable( "T_NJ_20",                  &NJ_20);
   reader->AddVariable( "T_E_rest_20",              &E_rest_20);
   reader->AddVariable( "T_HTT_rest_20",            &HTT_rest_20);

//   reader->AddVariable( "T_NJ_30",                  &NJ_30);
//   reader->AddVariable( "T_E_rest_30",              &E_rest_30);
//   reader->AddVariable( "T_HTT_rest_30",            &HTT_rest_30);

//   reader->AddVariable( "T_dETA_bb_qq",             &dETA_bb_qq);
   reader->AddVariable( "T_alphabb",                &alphabb);
   reader->AddVariable( "T_alphaqq",                &alphaqq);


   reader->AddVariable( "T_pt0",                    &pt0);
//   reader->AddVariable( "T_pt1",                    &pt1);
//   reader->AddVariable( "T_pt2",                    &pt2);
//   reader->AddVariable( "T_pt3",                    &pt3);
   reader->AddVariable( "T_pt4",                    &pt4);
//   reader->AddVariable( "T_pt5",                    &pt5);

//   reader->AddVariable( "T_btgq1_DeepCSV",          &BtgQ1_DeepCSV);
//   reader->AddVariable( "T_btgq1_DeepJet",          &BtgQ1_DeepJet);
//   reader->AddVariable( "T_btgq2_DeepCSV",          &BtgQ2_DeepCSV);
//   reader->AddVariable( "T_btgq2_DeepJet",          &BtgQ2_DeepJet);
//   reader->AddVariable( "T_ptb1",                   &ptB1);
//   reader->AddVariable( "T_ptb2",                   &ptB2);
   reader->AddVariable( "T_ptq1",                   &ptQ1);
   reader->AddVariable( "T_ptq2",                   &ptQ2);
   reader->AddVariable( "T_etaq1",                  &etaQ1);
   reader->AddVariable( "T_etaq2",                  &etaQ2);
//   reader->AddVariable( "T_phiq1",                  &phiQ1);
//   reader->AddVariable( "T_phiq2",                  &phiQ2);
   reader->AddVariable( "T_Ordb1",                  &OrdB1);
   reader->AddVariable( "T_Ordb2",                  &OrdB2);

   reader->AddVariable("T_zep",                     &zep);
   reader->AddVariable("T_phiA_bb_qq",              &phiA_bb_qq);
   reader->AddVariable("T_minR",                    &minR);
   reader->AddVariable("T_dR_leadqH",               &dR_leadqH);
   reader->AddVariable("T_dR_subleadqH",            &dR_subleadqH);
   reader->AddVariable("T_norm_pt",                 &norm_pt);
//   reader->AddVariable("T_prod_pt",                 &prod_pt);
//   reader->AddVariable("T_ptHtoM",                  &ptHtoM);

   reader->AddSpectator( "T_mbb",                   &MinvBB);

   reader->BookMVA("TMVAClassification_BDTG", "dataset/weights/TMVAClassification_BDTG.weights.xml" );

// ***************************************************************************************************************
// Signal
//   TFile *f00 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_VBFHToBB_M_125_.root");
   TFile *f00 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_VBFHToBB_M_125_dipoleRecoilOn_.root");
   TFile *f01 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_GluGluHToBB_M125_.root");
   TTree * ntps[2] = {(TTree*)f00->Get("Sig"), (TTree*)f01->Get("Sig")};
   int nSProc = sizeof(ntps)/sizeof(ntps[0]);

   Float_t S_weight[nSProc];
   Float_t S_pt0[nSProc],S_pt1[nSProc],S_pt2[nSProc],S_pt3[nSProc],S_pt4[nSProc],S_pt5[nSProc];
   Float_t S_ptB1[nSProc],S_ptB2[nSProc],S_etaB1[nSProc],S_etaB2[nSProc],S_phiB1[nSProc],S_phiB2[nSProc];
   Float_t S_BtgB1_DeepCSV[nSProc],S_BtgB1_DeepJet[nSProc],S_BtgB2_DeepCSV[nSProc],S_BtgB2_DeepJet[nSProc],S_OrdB1[nSProc],S_OrdB2[nSProc];
   Float_t S_ptQ1[nSProc],S_ptQ2[nSProc],S_etaQ1[nSProc],S_etaQ2[nSProc],S_phiQ1[nSProc],S_phiQ2[nSProc];
   Float_t S_BtgQ1_DeepCSV[nSProc],S_BtgQ1_DeepJet[nSProc],S_BtgQ2_DeepCSV[nSProc],S_BtgQ2_DeepJet[nSProc],S_QtgQ1[nSProc],S_QtgQ2[nSProc];
   Float_t S_MinvQQ[nSProc],S_dETAQQ[nSProc],S_dPHIQQ[nSProc],S_NJ_InGap[nSProc],S_ptAll[nSProc],S_pzAll[nSProc],S_dETA_bb_qq[nSProc],S_alphabb[nSProc],S_alphaqq[nSProc];
   Float_t S_NJ_20[nSProc],S_E_rest_20[nSProc],S_HTT_rest_20[nSProc],S_NJ_30[nSProc],S_E_rest_30[nSProc],S_HTT_rest_30[nSProc];
   Float_t S_zep[nSProc],S_phiA_bb_qq[nSProc],S_minR[nSProc],S_dR_leadqH[nSProc],S_dR_subleadqH[nSProc],S_norm_pt[nSProc],S_prod_pt[nSProc],S_ptHtoM[nSProc];
   Float_t S_MinvBB[nSProc];
   Float_t S_MinvBB_Reg[nSProc];

   for(int k=0;k<nSProc;k++)
   {

    ntps[k]->SetBranchAddress("T_weight",        &S_weight[k]);

    ntps[k]->SetBranchAddress("T_pt0",           &S_pt0[k]);
    ntps[k]->SetBranchAddress("T_pt1",           &S_pt1[k]);
    ntps[k]->SetBranchAddress("T_pt2",           &S_pt2[k]);
    ntps[k]->SetBranchAddress("T_pt3",           &S_pt3[k]);
    ntps[k]->SetBranchAddress("T_pt4",           &S_pt4[k]);
    ntps[k]->SetBranchAddress("T_pt5",           &S_pt5[k]);

    ntps[k]->SetBranchAddress("T_ptb1",          &S_ptB1[k]);
    ntps[k]->SetBranchAddress("T_ptb2",          &S_ptB2[k]);
    ntps[k]->SetBranchAddress("T_etab1",         &S_etaB1[k]);
    ntps[k]->SetBranchAddress("T_etab2",         &S_etaB2[k]);
    ntps[k]->SetBranchAddress("T_phib1",         &S_phiB1[k]);
    ntps[k]->SetBranchAddress("T_phib2",         &S_phiB2[k]);
    ntps[k]->SetBranchAddress("T_btgb1_DeepCSV", &S_BtgB1_DeepCSV[k]);
    ntps[k]->SetBranchAddress("T_btgb1_DeepJet", &S_BtgB1_DeepJet[k]);
    ntps[k]->SetBranchAddress("T_btgb2_DeepCSV", &S_BtgB2_DeepCSV[k]);
    ntps[k]->SetBranchAddress("T_btgb2_DeepJet", &S_BtgB2_DeepJet[k]);
    ntps[k]->SetBranchAddress("T_Ordb1",         &S_OrdB1[k]);
    ntps[k]->SetBranchAddress("T_Ordb2",         &S_OrdB2[k]);

    ntps[k]->SetBranchAddress("T_ptq1",          &S_ptQ1[k]);
    ntps[k]->SetBranchAddress("T_ptq2",          &S_ptQ2[k]);
    ntps[k]->SetBranchAddress("T_etaq1",         &S_etaQ1[k]);
    ntps[k]->SetBranchAddress("T_etaq2",         &S_etaQ2[k]);
    ntps[k]->SetBranchAddress("T_phiq1",         &S_phiQ1[k]);
    ntps[k]->SetBranchAddress("T_phiq2",         &S_phiQ2[k]);
    ntps[k]->SetBranchAddress("T_btgq1_DeepCSV", &S_BtgQ1_DeepCSV[k]);
    ntps[k]->SetBranchAddress("T_btgq1_DeepJet", &S_BtgQ1_DeepJet[k]);
    ntps[k]->SetBranchAddress("T_btgq2_DeepCSV", &S_BtgQ2_DeepCSV[k]);
    ntps[k]->SetBranchAddress("T_btgq2_DeepJet", &S_BtgQ2_DeepJet[k]);
    ntps[k]->SetBranchAddress("T_qtgq1",         &S_QtgQ1[k]);
    ntps[k]->SetBranchAddress("T_qtgq2",         &S_QtgQ2[k]);

    ntps[k]->SetBranchAddress("T_mqq",           &S_MinvQQ[k]);
    ntps[k]->SetBranchAddress("T_dETAqq",        &S_dETAQQ[k]);
    ntps[k]->SetBranchAddress("T_dPHIqq",        &S_dPHIQQ[k]);
    ntps[k]->SetBranchAddress("T_NJ_ingap",      &S_NJ_InGap[k]);
    ntps[k]->SetBranchAddress("T_ptAll",         &S_ptAll[k]);
    ntps[k]->SetBranchAddress("T_pzAll",         &S_pzAll[k]);
    ntps[k]->SetBranchAddress("T_NJ_20",         &S_NJ_20[k]);
    ntps[k]->SetBranchAddress("T_E_rest_20",     &S_E_rest_20[k]);
    ntps[k]->SetBranchAddress("T_HTT_rest_20",   &S_HTT_rest_20[k]);
    ntps[k]->SetBranchAddress("T_NJ_30",         &S_NJ_30[k]);
    ntps[k]->SetBranchAddress("T_E_rest_30",     &S_E_rest_30[k]);
    ntps[k]->SetBranchAddress("T_HTT_rest_30",   &S_HTT_rest_30[k]);
    ntps[k]->SetBranchAddress("T_HTT_rest_30",   &S_HTT_rest_30[k]);
    ntps[k]->SetBranchAddress("T_dETA_bb_qq",    &S_dETA_bb_qq[k]);
    ntps[k]->SetBranchAddress("T_alphabb",       &S_alphabb[k]);
    ntps[k]->SetBranchAddress("T_alphaqq",       &S_alphaqq[k]);

    ntps[k]->SetBranchAddress("T_zep",           &S_zep[k]);
    ntps[k]->SetBranchAddress("T_phiA_bb_qq",    &S_phiA_bb_qq[k]);
    ntps[k]->SetBranchAddress("T_minR",          &S_minR[k]);
    ntps[k]->SetBranchAddress("T_dR_leadqH",     &S_dR_leadqH[k]);
    ntps[k]->SetBranchAddress("T_dR_subleadqH",  &S_dR_subleadqH[k]);
    ntps[k]->SetBranchAddress("T_norm_pt",       &S_norm_pt[k]);
    ntps[k]->SetBranchAddress("T_prod_pt",       &S_prod_pt[k]);
    ntps[k]->SetBranchAddress("T_ptHtoM",        &S_ptHtoM[k]);

    ntps[k]->SetBranchAddress("T_mbb",           &S_MinvBB[k]);
    ntps[k]->SetBranchAddress("T_reg_mbb",       &S_MinvBB_Reg[k]);

   }

// QCD
   TFile *f0 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_QCD_HT100to200_TuneCP5_.root");
   TFile *f1 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_QCD_HT200to300_TuneCP5_.root");
   TFile *f2 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_QCD_HT300to500_TuneCP5_.root");
   TFile *f3 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_QCD_HT500to700_TuneCP5_.root");
   TFile *f4 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_QCD_HT700to1000_TuneCP5_.root");
   TFile *f5 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_QCD_HT1000to1500_TuneCP5_.root");
   TFile *f6 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_QCD_HT1500to2000_TuneCP5_.root");
   TFile *f7 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_QCD_HT2000toInf_TuneCP5_.root");
//Single_Top
   TFile *f8  = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_ST_tW_top_5f_inclusiveDecays_TuneCP5_.root");
   TFile *f9  = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_ST_tW_antitop_5f_inclusiveDecays_TuneCP5_.root");
   TFile *f10  = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_ST_t_channel_top_4f_InclusiveDecays_TuneCP5_.root");
   TFile *f11 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_ST_t_channel_antitop_4f_InclusiveDecays_TuneCP5_.root");
//ttbar
   TFile *f12 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_TTTo2L2Nu_TuneCP5_.root");
   TFile *f13 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_TTToHadronic_TuneCP5_.root");
   TFile *f14 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_TTToSemiLeptonic_TuneCP5_.root");
//WJets
   TFile *f15 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_WJetsToQQ_HT200to400_qc19_3j_TuneCP5_.root");
   TFile *f16 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_WJetsToQQ_HT400to600_qc19_3j_TuneCP5_.root");
   TFile *f17 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_WJetsToQQ_HT600to800_qc19_3j_TuneCP5_.root");
   TFile *f18 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_WJetsToQQ_HT800toInf_qc19_3j_TuneCP5_.root");
//ZJets
   TFile *f19 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_ZJetsToQQ_HT200to400_qc19_4j_TuneCP5_.root");
   TFile *f20 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_ZJetsToQQ_HT400to600_qc19_4j_TuneCP5_.root");
   TFile *f21 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_ZJetsToQQ_HT600to800_qc19_4j_TuneCP5_.root");
   TFile *f22 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/MC_Derived/NTuples_ZJetsToQQ_HT800toInf_qc19_4j_TuneCP5_.root");
// Full Data
   TFile *f23 = TFile::Open("../../../Making_Ntuples_For_SigBkgDiscrimination/Data_Derived/NTuples_Full_Data_2016_.root");

   TTree * ntpb[] =  { (TTree*)f0->Get("Bkg"),(TTree*)f1->Get("Bkg"),(TTree*)f2->Get("Bkg"),(TTree*)f3->Get("Bkg"),(TTree*)f4->Get("Bkg"),(TTree*)f5->Get("Bkg"),(TTree*)f6->Get("Bkg"), (TTree*)f7->Get("Bkg"),
                        (TTree*)f8->Get("Bkg"),(TTree*)f9->Get("Bkg"),(TTree*)f10->Get("Bkg"),(TTree*)f11->Get("Bkg"),
                        (TTree*)f12->Get("Bkg"),(TTree*)f13->Get("Bkg"),(TTree*)f14->Get("Bkg"),
                        (TTree*)f15->Get("Bkg"),(TTree*)f16->Get("Bkg"),(TTree*)f17->Get("Bkg"),(TTree*)f18->Get("Bkg"),
                        (TTree*)f19->Get("Bkg"),(TTree*)f20->Get("Bkg"),(TTree*)f21->Get("Bkg"),(TTree*)f22->Get("Bkg"),
                        (TTree*)f23->Get("Bkg")
                      };

   int nBProc = sizeof(ntpb)/sizeof(ntpb[0]);

   Float_t B_weight[nBProc];
   Float_t B_pt0[nBProc],B_pt1[nBProc],B_pt2[nBProc],B_pt3[nBProc],B_pt4[nBProc],B_pt5[nBProc];
   Float_t B_ptB1[nBProc],B_ptB2[nBProc],B_etaB1[nBProc],B_etaB2[nBProc],B_phiB1[nBProc],B_phiB2[nBProc];
   Float_t B_BtgB1_DeepCSV[nBProc],B_BtgB1_DeepJet[nBProc],B_BtgB2_DeepCSV[nBProc],B_BtgB2_DeepJet[nBProc],B_OrdB1[nBProc],B_OrdB2[nBProc];
   Float_t B_ptQ1[nBProc],B_ptQ2[nBProc],B_etaQ1[nBProc],B_etaQ2[nBProc],B_phiQ1[nBProc],B_phiQ2[nBProc];
   Float_t B_BtgQ1_DeepCSV[nBProc],B_BtgQ1_DeepJet[nBProc],B_BtgQ2_DeepCSV[nBProc],B_BtgQ2_DeepJet[nBProc],B_QtgQ1[nBProc],B_QtgQ2[nBProc];
   Float_t B_MinvQQ[nBProc],B_dETAQQ[nBProc],B_dPHIQQ[nBProc],B_NJ_InGap[nBProc],B_ptAll[nBProc],B_pzAll[nBProc],B_dETA_bb_qq[nBProc],B_alphabb[nBProc],B_alphaqq[nBProc];
   Float_t B_NJ_20[nBProc],B_E_rest_20[nBProc],B_HTT_rest_20[nBProc],B_NJ_30[nBProc],B_E_rest_30[nBProc],B_HTT_rest_30[nBProc];
   Float_t B_zep[nBProc],B_phiA_bb_qq[nBProc],B_minR[nBProc],B_dR_leadqH[nBProc],B_dR_subleadqH[nBProc],B_norm_pt[nBProc],B_prod_pt[nBProc],B_ptHtoM[nBProc];
   Float_t B_MinvBB[nBProc];
   Float_t B_MinvBB_Reg[nBProc];

   for(int k=0;k<nBProc;k++)
   {
    ntpb[k]->SetBranchAddress("T_weight",        &B_weight[k]);

    ntpb[k]->SetBranchAddress("T_pt0",           &B_pt0[k]);
    ntpb[k]->SetBranchAddress("T_pt1",           &B_pt1[k]);
    ntpb[k]->SetBranchAddress("T_pt2",           &B_pt2[k]);
    ntpb[k]->SetBranchAddress("T_pt3",           &B_pt3[k]);
    ntpb[k]->SetBranchAddress("T_pt4",           &B_pt4[k]);
    ntpb[k]->SetBranchAddress("T_pt5",           &B_pt5[k]);

    ntpb[k]->SetBranchAddress("T_ptb1",          &B_ptB1[k]);
    ntpb[k]->SetBranchAddress("T_ptb2",          &B_ptB2[k]);
    ntpb[k]->SetBranchAddress("T_etab1",         &B_etaB1[k]);
    ntpb[k]->SetBranchAddress("T_etab2",         &B_etaB2[k]);
    ntpb[k]->SetBranchAddress("T_phib1",         &B_phiB1[k]);
    ntpb[k]->SetBranchAddress("T_phib2",         &B_phiB2[k]);
    ntpb[k]->SetBranchAddress("T_btgb1_DeepCSV", &B_BtgB1_DeepCSV[k]);
    ntpb[k]->SetBranchAddress("T_btgb1_DeepJet", &B_BtgB1_DeepJet[k]);
    ntpb[k]->SetBranchAddress("T_btgb2_DeepCSV", &B_BtgB2_DeepCSV[k]);
    ntpb[k]->SetBranchAddress("T_btgb2_DeepJet", &B_BtgB2_DeepJet[k]);
    ntpb[k]->SetBranchAddress("T_Ordb1",         &B_OrdB1[k]);
    ntpb[k]->SetBranchAddress("T_Ordb2",         &B_OrdB2[k]);

    ntpb[k]->SetBranchAddress("T_ptq1",          &B_ptQ1[k]);
    ntpb[k]->SetBranchAddress("T_ptq2",          &B_ptQ2[k]);
    ntpb[k]->SetBranchAddress("T_etaq1",         &B_etaQ1[k]);
    ntpb[k]->SetBranchAddress("T_etaq2",         &B_etaQ2[k]);
    ntpb[k]->SetBranchAddress("T_phiq1",         &B_phiQ1[k]);
    ntpb[k]->SetBranchAddress("T_phiq2",         &B_phiQ2[k]);
    ntpb[k]->SetBranchAddress("T_btgq1_DeepCSV", &B_BtgQ1_DeepCSV[k]);
    ntpb[k]->SetBranchAddress("T_btgq1_DeepJet", &B_BtgQ1_DeepJet[k]);
    ntpb[k]->SetBranchAddress("T_btgq2_DeepCSV", &B_BtgQ2_DeepCSV[k]);
    ntpb[k]->SetBranchAddress("T_btgq2_DeepJet", &B_BtgQ2_DeepJet[k]);
    ntpb[k]->SetBranchAddress("T_qtgq1",         &B_QtgQ1[k]);
    ntpb[k]->SetBranchAddress("T_qtgq2",         &B_QtgQ2[k]);

    ntpb[k]->SetBranchAddress("T_mqq",           &B_MinvQQ[k]);
    ntpb[k]->SetBranchAddress("T_dETAqq",        &B_dETAQQ[k]);
    ntpb[k]->SetBranchAddress("T_dPHIqq",        &B_dPHIQQ[k]);
    ntpb[k]->SetBranchAddress("T_NJ_ingap",      &B_NJ_InGap[k]);
    ntpb[k]->SetBranchAddress("T_ptAll",         &B_ptAll[k]);
    ntpb[k]->SetBranchAddress("T_pzAll",         &B_pzAll[k]);
    ntpb[k]->SetBranchAddress("T_NJ_20",         &B_NJ_20[k]);
    ntpb[k]->SetBranchAddress("T_E_rest_20",     &B_E_rest_20[k]);
    ntpb[k]->SetBranchAddress("T_HTT_rest_20",   &B_HTT_rest_20[k]);
    ntpb[k]->SetBranchAddress("T_NJ_30",         &B_NJ_30[k]);
    ntpb[k]->SetBranchAddress("T_E_rest_30",     &B_E_rest_30[k]);
    ntpb[k]->SetBranchAddress("T_HTT_rest_30",   &B_HTT_rest_30[k]);
    ntpb[k]->SetBranchAddress("T_HTT_rest_30",   &B_HTT_rest_30[k]);
    ntpb[k]->SetBranchAddress("T_dETA_bb_qq",    &B_dETA_bb_qq[k]);
    ntpb[k]->SetBranchAddress("T_alphabb",       &B_alphabb[k]);
    ntpb[k]->SetBranchAddress("T_alphaqq",       &B_alphaqq[k]);

    ntpb[k]->SetBranchAddress("T_zep",           &B_zep[k]);
    ntpb[k]->SetBranchAddress("T_phiA_bb_qq",    &B_phiA_bb_qq[k]);
    ntpb[k]->SetBranchAddress("T_minR",          &B_minR[k]);
    ntpb[k]->SetBranchAddress("T_dR_leadqH",     &B_dR_leadqH[k]);
    ntpb[k]->SetBranchAddress("T_dR_subleadqH",  &B_dR_subleadqH[k]);
    ntpb[k]->SetBranchAddress("T_norm_pt",       &B_norm_pt[k]);
    ntpb[k]->SetBranchAddress("T_prod_pt",       &B_prod_pt[k]);
    ntpb[k]->SetBranchAddress("T_ptHtoM",        &B_ptHtoM[k]);

    ntpb[k]->SetBranchAddress("T_mbb",           &B_MinvBB[k]);
    ntpb[k]->SetBranchAddress("T_reg_mbb",       &B_MinvBB_Reg[k]);

   }
//***************************************************************************************
  TFile *f = new TFile("Mass_in_BDT_Output.root","RECREATE");
  gStyle->SetOptFile(0);
  gStyle->SetOptStat("mre");
  gStyle->SetPaintTextFormat("1.2e");

  double bdt_Cut_MIN=0.71;
  double bdt_Cut_MAX=1.0;
  double xmin=90;
  double xmax=170;
  int    NBINS=15;


  TH1D * _Mass_in_BDT_VBF_Hbb  = new TH1D("Mass_in_BDT_VBF_Hbb","Mass_in_BDT_VBF_Hbb",NBINS,xmin,xmax);
  TH1D * _Mass_in_BDT_ggF_Hbb  = new TH1D("Mass_in_BDT_ggF_Hbb","Mass_in_BDT_ggF_Hbb",NBINS,xmin,xmax);
  TH1D * _Mass_in_BDT_QCD  = new TH1D("Mass_in_BDT_QCD","Mass_in_BDT_QCD",NBINS,xmin,xmax);
  TH1D * _Mass_in_BDT_Single_Top  = new TH1D("Mass_in_BDT_Single_Top","Mass_in_BDT_Single_Top",NBINS,xmin,xmax);
  TH1D * _Mass_in_BDT_tt  = new TH1D("Mass_in_BDT_tt","Mass_in_BDT_tt",NBINS,xmin,xmax);
  TH1D * _Mass_in_BDT_WJets  = new TH1D("Mass_in_BDT_WJets","Mass_in_BDT_WJets",NBINS,xmin,xmax);
  TH1D * _Mass_in_BDT_ZJets  = new TH1D("Mass_in_BDT_ZJets","Mass_in_BDT_ZJets",NBINS,xmin,xmax);
  TH1D * _Mass_in_BDT_DATA   = new TH1D("Mass_in_BDT_DATA","Mass_in_BDT_DATA",NBINS,xmin,xmax);
//***************************************************************************************

  for(int k=0;k<nSProc;k++)
   for(int i=0;i<ntps[k]->GetEntries();i++)
    {
     ntps[k]->GetEntry(i);

     pt0=S_pt0[k]; pt1=S_pt1[k]; pt2=S_pt2[k]; pt3=S_pt3[k]; pt4=S_pt4[k]; pt5=S_pt5[k];
     ptB1=S_ptB1[k]; ptB2=S_ptB2[k]; etaB1=S_etaB1[k]; etaB2=S_etaB2[k]; phiB1=S_phiB1[k]; phiB2=S_phiB2[k];
     BtgB1_DeepCSV=S_BtgB1_DeepCSV[k]; BtgB1_DeepJet=S_BtgB1_DeepJet[k]; BtgB2_DeepCSV=S_BtgB2_DeepCSV[k]; BtgB2_DeepJet=S_BtgB2_DeepJet[k]; OrdB1=S_OrdB1[k]; OrdB2=S_OrdB2[k];
     ptQ1=S_ptQ1[k]; ptQ2=S_ptQ2[k]; etaQ1=S_etaQ1[k]; etaQ2=S_etaQ2[k]; phiQ1=S_phiQ1[k]; phiQ2=S_phiQ2[k];
     BtgQ1_DeepCSV=S_BtgQ1_DeepCSV[k]; BtgQ1_DeepJet=S_BtgQ1_DeepJet[k]; BtgQ2_DeepCSV=S_BtgQ2_DeepCSV[k]; BtgQ2_DeepJet=S_BtgQ2_DeepJet[k]; QtgQ1=S_QtgQ1[k]; QtgQ2=S_QtgQ2[k];
     MinvQQ=S_MinvQQ[k]; dETAQQ=S_dETAQQ[k]; dPHIQQ=S_dPHIQQ[k]; NJ_InGap=S_NJ_InGap[k]; ptAll=S_ptAll[k]; pzAll=S_pzAll[k];
     NJ_20=S_NJ_20[k]; E_rest_20=S_E_rest_20[k]; HTT_rest_20=S_HTT_rest_20[k]; NJ_30=S_NJ_30[k]; E_rest_30=S_E_rest_30[k]; HTT_rest_30=S_HTT_rest_30[k];
     dETA_bb_qq=S_dETA_bb_qq[k]; alphabb=S_alphabb[k]; alphaqq=S_alphaqq[k];
     zep=S_zep[k]; phiA_bb_qq=S_phiA_bb_qq[k]; minR=S_minR[k]; dR_leadqH=S_dR_leadqH[k]; dR_subleadqH=S_dR_subleadqH[k]; norm_pt=S_norm_pt[k]; prod_pt=S_prod_pt[k]; ptHtoM=S_ptHtoM[k];
     MinvBB=S_MinvBB[k];
     MinvBB_Reg=S_MinvBB_Reg[k];

     if(reader->EvaluateMVA("TMVAClassification_BDTG")<bdt_Cut_MIN) continue;
     if(reader->EvaluateMVA("TMVAClassification_BDTG")>bdt_Cut_MAX) continue;
     if (QtgQ1!=QtgQ1) continue;

     if(k==0)
      _Mass_in_BDT_VBF_Hbb->Fill(S_MinvBB_Reg[k],S_weight[k]);
     if(k==1)
      _Mass_in_BDT_ggF_Hbb->Fill(S_MinvBB_Reg[k],S_weight[k]);
    }

  for(int k=0;k<nBProc;k++)
   for(int i=0;i<ntpb[k]->GetEntries();i++)
    {
     ntpb[k]->GetEntry(i);

     pt0=B_pt0[k]; pt1=B_pt1[k]; pt2=B_pt2[k]; pt3=B_pt3[k]; pt4=B_pt4[k]; pt5=B_pt5[k];
     ptB1=B_ptB1[k]; ptB2=B_ptB2[k]; etaB1=B_etaB1[k]; etaB2=B_etaB2[k]; phiB1=B_phiB1[k]; phiB2=B_phiB2[k];
     BtgB1_DeepCSV=B_BtgB1_DeepCSV[k]; BtgB1_DeepJet=B_BtgB1_DeepJet[k]; BtgB2_DeepCSV=B_BtgB2_DeepCSV[k]; BtgB2_DeepJet=B_BtgB2_DeepJet[k]; OrdB1=B_OrdB1[k]; OrdB2=B_OrdB2[k];
     ptQ1=B_ptQ1[k]; ptQ2=B_ptQ2[k]; etaQ1=B_etaQ1[k]; etaQ2=B_etaQ2[k]; phiQ1=B_phiQ1[k]; phiQ2=B_phiQ2[k];
     BtgQ1_DeepCSV=B_BtgQ1_DeepCSV[k]; BtgQ1_DeepJet=B_BtgQ1_DeepJet[k]; BtgQ2_DeepCSV=B_BtgQ2_DeepCSV[k]; BtgQ2_DeepJet=B_BtgQ2_DeepJet[k]; QtgQ1=B_QtgQ1[k]; QtgQ2=B_QtgQ2[k];
     MinvQQ=B_MinvQQ[k]; dETAQQ=B_dETAQQ[k]; dPHIQQ=B_dPHIQQ[k]; NJ_InGap=B_NJ_InGap[k]; ptAll=B_ptAll[k]; pzAll=B_pzAll[k];
     NJ_20=B_NJ_20[k]; E_rest_20=B_E_rest_20[k]; HTT_rest_20=B_HTT_rest_20[k]; NJ_30=B_NJ_30[k]; E_rest_30=B_E_rest_30[k]; HTT_rest_30=B_HTT_rest_30[k];
     dETA_bb_qq=B_dETA_bb_qq[k]; alphabb=B_alphabb[k]; alphaqq=B_alphaqq[k];
     zep=B_zep[k]; phiA_bb_qq=B_phiA_bb_qq[k]; minR=B_minR[k]; dR_leadqH=B_dR_leadqH[k]; dR_subleadqH=B_dR_subleadqH[k]; norm_pt=B_norm_pt[k]; prod_pt=B_prod_pt[k]; ptHtoM=B_ptHtoM[k];
     MinvBB=B_MinvBB[k];
     MinvBB_Reg=B_MinvBB_Reg[k];

     if(reader->EvaluateMVA("TMVAClassification_BDTG")<bdt_Cut_MIN) continue;
     if(reader->EvaluateMVA("TMVAClassification_BDTG")>bdt_Cut_MAX) continue;
     if (QtgQ1!=QtgQ1) continue;

     if(k<8)
      _Mass_in_BDT_QCD->Fill(B_MinvBB_Reg[k],B_weight[k]);
     else if(k>=8 && k<12)
      _Mass_in_BDT_Single_Top->Fill(B_MinvBB_Reg[k],B_weight[k]);
     else if(k>=12 && k<15)
      _Mass_in_BDT_tt->Fill(B_MinvBB_Reg[k],B_weight[k]);
     else if(k>=15 && k<19)
      _Mass_in_BDT_WJets->Fill(B_MinvBB_Reg[k],B_weight[k]);
     else if(k>=19 && k<23)
      _Mass_in_BDT_ZJets->Fill(B_MinvBB_Reg[k],B_weight[k]);
     else if(k==23)
      _Mass_in_BDT_DATA->Fill(B_MinvBB_Reg[k],B_weight[k]);
    }

  f->Write();
  delete f00;
  delete f01;

  delete f0;
  delete f1;
  delete f2;
  delete f3;
  delete f4;
  delete f5;
  delete f6;
  delete f7;
  delete f8;
  delete f9;
  delete f10;
  delete f11;
  delete f12;
  delete f13;
  delete f14;
  delete f15;
  delete f16;
  delete f17;
  delete f18;
  delete f19;
  delete f20;
  delete f21;
  delete f22;
  delete f23;
}