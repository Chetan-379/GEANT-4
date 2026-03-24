#define AnalyzeLightBSM_cxx
#include "AnalyzeLightBSM.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <cstring>
#include <string>
#include <fstream>
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "Math/Vector3D.h"
#include <algorithm> // for std::shuffle
#include <random>    // for std::default_random_engine
#include <array>

#pragma link C++ class std::vector < std::vector> + ;
#pragma link C++ class std::vector < TLorentzVector> + ;
#ifdef __CINT__
#pragma link C++ class std::vector < TLorentzVector> + ;
#endif
using namespace TMVA;
using Array3D = std::array<std::array<std::array<float, 8>, 8>, 7>;
int main(int argc, char *argv[])
{

  if (argc < 3)
    {
      cerr << "Please give 3 arguments " << "inputFileList " << " " << "outputFileName" << " " << "DetectorType" << endl;
      return -1;
    }
  const char *inputFileList = argv[1];
  const char *outFileName = argv[2];
  const char *detType = argv[3];

  AnalyzeLightBSM ana(inputFileList, outFileName, detType);
  cout << "Detector type " << detType << endl;
  cout << inputFileList << "\t" << outFileName << "\t" << detType << endl;
  ana.EventLoop(detType, inputFileList, outFileName);
  Tools::Instance();
  return 0;
}

void AnalyzeLightBSM::EventLoop(const char *detType, const char *inputFileList, const char *outFileName)
{
  TString InFile = TString(inputFileList);
  
  cout << "inside event loop" << endl;
  if (fChain == 0)
    return;
  Long64_t nentries = fChain->GetEntriesFast();
  cout << "nentries " << nentries << endl;
  cout << "Analyzing data from " << detType << " " << endl;

  Long64_t nbytes = 0, nb = 0;
  int decade = 0;
  bool Debug = false;

  int nOpG_q12 =-99, nOpG_q13 =-99, nOpL_q12 =-99, nOpL_q13 =-99;
  int nTrig =0, ntwoHit=0, nFullEdep=0, nDiffClr=0, nDiffPix=0;
  int scat_idx =-999, abs_idx =-999, scat_idx_noCCal =-999, abs_idx_noCCal =-999;
  float scat_Edep = -99, abs_Edep = -99, scat_theta = -999, scat_Edep_noCCal = -99, abs_Edep_noCCal = -99, scat_theta_noCCal = -999;
  
  for (Long64_t jentry = 0; jentry < nentries; jentry++)
    {
      double progress = 10.0 * jentry / (1.0 * nentries);
      int k = int(progress);
      if (k > decade)
	cout << 10 * k << " %" << endl;
      decade = k;

      // ===============read this entry == == == == == == == == == == ==
    
      Long64_t ientry = LoadTree(jentry);
	
      if(Debug) cout<<"===load tree entry ==="<<jentry<<endl;
      if (ientry < 0) break;

      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if(Debug)	cout<<"===load tree entry ==="<<jentry<<endl;

      //Main analyzer code---------    
      h_nOp_GAGG->Fill(nOp_GAGG);
      h_nOp_LYSO->Fill(nOp_LYSO);

      h_nOp_GAGG_QE->Fill(nOp_GAGG_QE);
      h_nOp_LYSO_QE->Fill(nOp_LYSO_QE);

      h_nOp_GAGG_truth->Fill(nOp_GAGG_truth);
      h_nOp_LYSO_truth->Fill(nOp_LYSO_truth);

      h_nOp_GvsL_truth->Fill(nOp_GAGG_truth, nOp_LYSO_truth);
    
      h_Edep_truth->Fill(Edep_truth);
      h_nOp_GL_sum->Fill(nOp_GAGG + nOp_LYSO);

      h_EdepG_truth->Fill(Edep_GAGG_truth);
      h_EdepL_truth->Fill(Edep_LYSO_truth);

      h_nOp_GvsL->Fill(nOp_GAGG, nOp_LYSO);
      h_nOp_GvsL_truth->Fill(nOp_GAGG_truth, nOp_LYSO_truth);
    
      h_nOpG_vs_Edep->Fill(Edep_truth, nOp_GAGG);
      h_nOpL_vs_Edep->Fill(Edep_truth, nOp_LYSO);
      h_nOpSumGL_vs_Edep->Fill(Edep_truth, nOp_GAGG + nOp_LYSO);

      h_nOpG_vs_EdepG->Fill(Edep_GAGG_truth, nOp_GAGG);
      h_nOpG_truth_vs_EdepG->Fill(Edep_GAGG_truth, nOp_GAGG_truth);
    
      h_nOpL_vs_EdepL->Fill(Edep_LYSO_truth, nOp_LYSO);
      h_nOpL_truth_vs_EdepL->Fill(Edep_LYSO_truth, nOp_LYSO_truth);

      //after QE
      h_nOpG_QE_vs_EdepG->Fill(Edep_GAGG_truth, nOp_GAGG_QE);
      h_nOpL_QE_vs_EdepL->Fill(Edep_LYSO_truth, nOp_LYSO_QE);

      //filling histograms quadrants wise
      for (int iqd=0; iqd < 4; iqd++){
	h_nOpG_quad[iqd]->Fill(nOpGAGG_quad[iqd]);
	h_nOpL_quad[iqd]->Fill(nOpLYSO_quad[iqd]);
      }

      //comparing diagonal vs adjecent
      nOpG_q12 = nOpGAGG_quad[qd1]+nOpGAGG_quad[qd2];
      nOpG_q13 = nOpGAGG_quad[qd1]+nOpGAGG_quad[qd3];

      nOpL_q12 = nOpLYSO_quad[qd1]+nOpLYSO_quad[qd2];
      nOpL_q13 = nOpLYSO_quad[qd1]+nOpLYSO_quad[qd3];

      h_nOpG_q12->Fill(nOpG_q12);
      h_nOpG_q13->Fill(nOpG_q13);
      h_nOpL_q12->Fill(nOpL_q12);
      h_nOpL_q13->Fill(nOpL_q13);
    
      h_nOpG_q12_vs_q13->Fill(nOpG_q12, nOpG_q13);
      h_nOpL_q12_vs_q13->Fill(nOpL_q12, nOpL_q13);

      //######################################################################
      //                   Analysing the Info w/o opical photons
      //######################################################################
    
      float Edep_truth_cryst = 0., scat_theta_sim = -99999.;
      int clr_Idx = -99, row_Idx = -99, col_Idx = -99, nCompt = -99, nPhoto = -99;
    
      vector<double> Hit_Edep, Hit_scat_theta_sim;
      vector<int> Hit_rowIdx, Hit_colIdx, Hit_clrIdx, Hit_nCompt, Hit_nPhoto;

      for(int i =0; i<Edep_truth_vec->size(); i++){
	Edep_truth_cryst = (*Edep_truth_vec)[i];
	clr_Idx = (*Det_clr_Idx)[i];
	row_Idx = (*Det_row_Idx)[i];
	col_Idx = (*Det_col_Idx)[i];
	nCompt = (*Det_nCompt)[i];
	nPhoto = (*Det_nPhoto)[i];
	scat_theta_sim = (*Det_scat_theta)[i];

	h_Edep_truth_cryst[noCut]->Fill(Edep_truth_cryst);
	h_theta_Sim[noCut]->Fill(Deg(scat_theta_sim));

	if(clr_Idx ==1) h_EdepG_truth_cryst->Fill(Edep_truth_cryst);
	if(clr_Idx ==2) h_EdepL_truth_cryst->Fill(Edep_truth_cryst);

	h_rowVScol->Fill(col_Idx, row_Idx);
      

	//-------------
	if(Edep_truth_cryst >= 0.06 && Edep_truth_cryst <=0.405){        
	  Hit_Edep.push_back(Edep_truth_cryst);
	  Hit_rowIdx.push_back(row_Idx);
	  Hit_colIdx.push_back(col_Idx);
	  Hit_clrIdx.push_back(clr_Idx);
	  Hit_nCompt.push_back(nCompt);
	  Hit_nPhoto.push_back(nPhoto);
	  Hit_scat_theta_sim.push_back(scat_theta_sim);

	  h_Edep_truth_cryst[Trig]->Fill(Edep_truth_cryst);
	  h_theta_Sim[Trig]->Fill(Deg(scat_theta_sim));
	}
      }

      if (Hit_Edep.size()>0) nTrig++; 
    
      if (Hit_Edep.size() ==2){
	ntwoHit++;
	h_Edep_truth_cryst[nhit2]->Fill(Hit_Edep[0]);
	h_Edep_truth_cryst[nhit2]->Fill(Hit_Edep[1]);
	h_theta_Sim[nhit2]->Fill(Deg(Hit_scat_theta_sim[0]));
	h_theta_Sim[nhit2]->Fill(Deg(Hit_scat_theta_sim[1]));
      
	if(Hit_Edep[0]+Hit_Edep[1] >0.5109 && Hit_Edep[0]+Hit_Edep[1] < 0.5111) {
	  nFullEdep++;
	  h_Edep_truth_cryst[fullEdep]->Fill(Hit_Edep[0]);
	  h_Edep_truth_cryst[fullEdep]->Fill(Hit_Edep[1]);

	  h_theta_Sim[fullEdep]->Fill(Deg(Hit_scat_theta_sim[0]));
	  h_theta_Sim[fullEdep]->Fill(Deg(Hit_scat_theta_sim[1]));
	
	  if(Hit_clrIdx[0] != Hit_clrIdx[1]) {
	    nDiffClr++;
	    h_Edep_truth_cryst[diffClr]->Fill(Hit_Edep[0]);
	    h_Edep_truth_cryst[diffClr]->Fill(Hit_Edep[1]);

	    h_theta_Sim[diffClr]->Fill(Deg(Hit_scat_theta_sim[0]));
	    h_theta_Sim[diffClr]->Fill(Deg(Hit_scat_theta_sim[1]));
	  
	    if((Hit_rowIdx[0] != Hit_rowIdx[1]) || (Hit_colIdx[0] != Hit_colIdx[1])){
	      nDiffPix++;

	      //--------
	      if(Hit_clrIdx[0] == 1 && Hit_clrIdx[1] ==2){
		//hypothesis 1: 0 is scatterer
		if (Hit_Edep[0] <= 0.255){    //forward scattering
		  scat_idx = 0;
		  abs_idx = 1;
		}

		//hypothesis 2: 1 is scatterer
		else if (Hit_Edep[1] > 0.255 && Hit_Edep[0] <= 0.34067){   //back scattering
		  scat_idx = 1;
		  abs_idx = 0;		
		}
	      }
	    
	      else if(Hit_clrIdx[0] == 2 && Hit_clrIdx[1] ==1){
		//hypothesis 1: 0 is scatterer
		if (Hit_Edep[0] > 0.255 && Hit_Edep[0] <= 0.34067){    //back scattering
		  scat_idx = 0;
		  abs_idx = 1;
		}

		//hypothesis 2: 1 is scatterer
		else if (Hit_Edep[1] <= 0.255){  //forward scattering
		  scat_idx = 1;
		  abs_idx = 0;		
		}
	      }
	      //---------------

	      if (scat_idx >=0){
		scat_Edep = Hit_Edep[scat_idx];

		if(scat_Edep < 0.34067){	      
		  scat_theta = GetScatAngle(scat_Edep);	      
	      
		  h_theta_ana->Fill(Deg(scat_theta));
		  h_theta_sim->Fill(Deg(Hit_scat_theta_sim[scat_idx]));
		  h_theta_SimVsAna->Fill(Deg(Hit_scat_theta_sim[scat_idx]), Deg(scat_theta));

		  double dTheta = Hit_scat_theta_sim[scat_idx]-scat_theta;
		  if(scat_theta >=0 && Hit_scat_theta_sim[scat_idx] >= 0) h_dTheta->Fill(Deg(abs(dTheta)));

		  if(Hit_nCompt[scat_idx] ==1 && Hit_nPhoto[scat_idx] == 0){		
		    h_theta_ana_nInt->Fill(Deg(scat_theta));
		    h_theta_sim_nInt->Fill(Deg(Hit_scat_theta_sim[scat_idx]));
		    h_theta_SimVsAna_nInt->Fill(Deg(Hit_scat_theta_sim[scat_idx]), Deg(scat_theta));				  		
		    if(scat_theta >=0 && Hit_scat_theta_sim[scat_idx] >= 0) h_dTheta_nInt->Fill(Deg(abs(dTheta)));
		  }
		}	      
	      }	      
	    }	  
	  }
      
      
	  //-----------------Indentifying the scatterer for NoCCal case----------------	
	  if(InFile.Contains("noCCal")) {
	
	    if (Hit_Edep[0] < Hit_Edep[1]){	  	   
	      scat_idx_noCCal = 0;
	      abs_idx_noCCal = 1;
	    }	     
    	    
	    else if (Hit_Edep[0] > Hit_Edep[1]){	  	   
	      scat_idx_noCCal = 1;
	      abs_idx_noCCal = 0;
	    }

	    if (scat_idx_noCCal >=0){
	      scat_Edep_noCCal = Hit_Edep[scat_idx_noCCal];
	  
	      if(scat_Edep_noCCal < 0.34067){	      
		scat_theta_noCCal = GetScatAngle(scat_Edep_noCCal);	      
	    
		h_theta_ana_noCCal->Fill(Deg(scat_theta_noCCal));
		h_theta_sim_noCCal->Fill(Deg(Hit_scat_theta_sim[scat_idx_noCCal]));
		h_theta_SimVsAna_noCCal->Fill(Deg(Hit_scat_theta_sim[scat_idx_noCCal]), Deg(scat_theta_noCCal));
		double dTheta_noCCal = Hit_scat_theta_sim[scat_idx_noCCal]-scat_theta_noCCal;
		if(scat_theta_noCCal >=0 && Hit_scat_theta_sim[scat_idx_noCCal] >= 0) h_dTheta_noCCal->Fill(Deg(abs(dTheta_noCCal)));
	    
		if(Hit_nCompt[scat_idx_noCCal] ==1 && Hit_nPhoto[scat_idx_noCCal] == 0){		
		  h_theta_ana_nInt_noCCal->Fill(Deg(scat_theta_noCCal));
		  h_theta_sim_nInt_noCCal->Fill(Deg(Hit_scat_theta_sim[scat_idx_noCCal]));
		  h_theta_SimVsAna_nInt_noCCal->Fill(Deg(Hit_scat_theta_sim[scat_idx_noCCal]), Deg(scat_theta_noCCal));				  		
		  if(scat_theta_noCCal >=0 && Hit_scat_theta_sim[scat_idx_noCCal] >= 0) h_dTheta_nInt_noCCal->Fill(Deg(abs(dTheta_noCCal)));
		}
	      }	  
	    }
	  }
        
	}
      }
    
    } // jentry loop end

  //normalizing the difference hists------------
  if(InFile.Contains("_CCal_")) {
    h_dTheta->Scale(1.0 / h_dTheta->Integral());
    h_dTheta_nInt->Scale(1.0 / h_dTheta_nInt->Integral());
  }

  if(InFile.Contains("noCCal")) {
    h_dTheta_noCCal->Scale(1.0 / h_dTheta_noCCal->Integral());
    h_dTheta_nInt_noCCal->Scale(1.0 / h_dTheta_nInt_noCCal->Integral());
  }
  //--------------------------------------------
  
  cout << "Total Triggerd evts: " << nTrig/1000. << "%" << endl;
  cout << "two hit evts: " << ntwoHit*100./nTrig << "%" << endl;
  cout << "Full Edep evts in two hits: " << nFullEdep*100./nTrig << "%" << endl;
  cout << "Events with full energy deposits in different clrs: " << nDiffClr*100./nTrig << "%" << endl;
  cout << "Events with full energy deposits in different clrs and pixels: " << nDiffPix*100./nTrig << "%" << endl;

}

string AnalyzeLightBSM::MatName(int DetId)
{
  string mat = " ";
  if (DetId < 0)
    mat = "error";
  if (DetId / 4000 >= 1)
    mat = "BGO";
  else if (DetId / 1000 >= 1 && DetId / 4000 < 1)
    mat = "Scint";

  return mat;
}

float AnalyzeLightBSM::GetScatAngle(float E_dep)
{
  double scatAngle = acos(1 - (0.511 * ((1 / (0.511 - E_dep)) - (1 / 0.511))));
  return scatAngle;
}

float AnalyzeLightBSM::GetScatElecE(float theta)
{
  double elec_E = 0.511 - (0.511 / 2 - cos(theta * Pi / 180));
  return elec_E;
}

float AnalyzeLightBSM::Deg(float Rad)
{
  double degree = Rad * 180 / Pi;
  return degree;
}

float AnalyzeLightBSM::Rad(float Deg)
{
  double radian = Deg * Pi / 180;
  return radian;
}









//================================================================================================================================================================================================================

// some baks
