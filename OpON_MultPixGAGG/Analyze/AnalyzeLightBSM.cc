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
#include"TGraphErrors.h"
#include"TGraphAsymmErrors.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "Math/Vector3D.h"
#include <algorithm>   // for std::shuffle
#include <random>      // for std::default_random_engine
#include <array>

#pragma link C++ class std::vector< std::vector >+; 
#pragma link C++ class std::vector< TLorentzVector >+;
#ifdef __CINT__
#pragma link C++ class std::vector<TLorentzVector>+;
#endif
using namespace TMVA;
using Array3D = std::array<std::array<std::array<float, 8>, 8>, 7>;
int main(int argc, char* argv[])
{

  if (argc < 3) {
    cerr << "Please give 3 arguments " << "inputFileList " << " " << "outputFileName" << " " << "DetectorType" <<endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *detType       = argv[3];

  AnalyzeLightBSM ana(inputFileList, outFileName,detType);
  cout << "Detector type " << detType << endl;
  cout<<inputFileList <<"\t"<<outFileName<<"\t"<<detType<<endl;
  ana.EventLoop(detType,inputFileList,outFileName);
  Tools::Instance();
  return 0;
}

void AnalyzeLightBSM::EventLoop(const char *detType,const char *inputFileList, const char *outFileName) {
  cout<<"inside event loop"<<endl;
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  cout << "nentries " << nentries << endl;
  cout << "Analyzing data from " << detType << " " << endl;

  Long64_t nbytes = 0, nb = 0;
  int decade = 0;
  bool Debug=false;
  int Full_edep_evts = 0;
  int Comp2_evts = 0;
  int nDiffAngle = 0, nHits = -1, nRayl = 0;


  int nEvts2Hits = 0, check_top =0, check_bottom =0;
  

  vector<vector<float>> phi_vec, theta_vec;
  phi_vec.resize(2);
  theta_vec.resize(2);

  int Trig =0, full_Edep = 0, sing_pix =0, two_pix =0, mult_pix =0;
  
  vector<vector<float>> left_Hit_Edep_vec, left_Hit_DetX_vec, left_Hit_DetY_vec, right_Hit_Edep_vec, right_Hit_DetX_vec, right_Hit_DetY_vec;

  vector<vector<float>> left_Hit_Edep_shuff, left_Hit_DetX_shuff, left_Hit_DetY_shuff, right_Hit_Edep_shuff, right_Hit_DetX_shuff, right_Hit_DetY_shuff;

  int check_reco =0, check_true =0;

  right_Hit_Edep_vec.clear();
  right_Hit_DetX_vec.clear();
  right_Hit_DetY_vec.clear();
  right_Hit_Edep_shuff.clear();
  right_Hit_DetX_shuff.clear();
  right_Hit_DetY_shuff.clear();

  left_Hit_Edep_vec.clear();
  left_Hit_DetX_vec.clear();
  left_Hit_DetY_vec.clear();
  left_Hit_Edep_shuff.clear();
  left_Hit_DetX_shuff.clear();
  left_Hit_DetY_shuff.clear();
   
  for (Long64_t jentry=0; jentry<nentries;jentry++){
    //cout << "----------event No.: " << jentry+1 << "-----------" << endl;
    int nHitComp = 0;
    double progress = 10.0 * jentry / (1.0 * nentries);
    int k = int (progress);
    if (k > decade)
      cout << 10 * k << " %" << endl;
    decade = k;

    // ===============read this entry == == == == == == == == == == ==                                                                        
    Long64_t ientry = LoadTree(jentry);
	
    if(Debug) cout<<"===load tree entry ==="<<jentry<<endl;
    if (ientry < 0) break;

    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if(Debug)	cout<<"===load tree entry ==="<<jentry<<endl;

    h_Compt_Edep-> Fill(Edep_Compt);
    h_Photo_Edep-> Fill(Edep_Photo);
    h_ComptVsPhoto_Edep->Fill(Edep_Compt, Edep_Photo);
    h_Total_Edep-> Fill(Total_Edep);

    nHits = Hit_Time->size();

    ROOT::Math::XYZVector v_comp, vin, vout, v1, v2;


      
    //--------------------Analysing the Scintillation photons at the end---------------------------

    //------------ Making an array of modules-----------------
    
    float Modules[2][8][8][7];
    for (int i =0; i<8; i++){
      for (int j =0; j<8; j++){
	for (int k =0; k<7; k++){	    
	  Modules[r][i][j][k] = right_module[i][j][k];
	  Modules[l][i][j][k] = left_module[i][j][k];
	}
      }
    }
    //--------------------------------------------------------
      
    int nOp_mdl[2] = {0}, nOp_gen_mdl[2] = {0};
    int nOp_crys[2] = {0}, nOp_gen_crys[2] = {0};
    float EDep_truth_mdl[2] = {0.}, EDep_reco_mdl[2] = {0.}, EDep_crys[2] = {0.}, Edep_reco_crys[2] = {0.}, DetX[2] = {0.}, DetY[2] = {0.};
    float Edep_diff[2] = {0.};      
    float E_scatElec[2] = {0.};

    int nHits_mdl_reco[2] = {0}, nHits_mdl_truth[2] = {0};

    float sigma =0;

    TRandom3 rng(12345);
    float sigma_spread = 0.;
       


    vector<vector<float>> Hit_Edep_reco, Hit_Edep_truth, Hit_Scat_elecE, Hit_DetX_reco, Hit_DetY_reco, Hit_DetX_truth, Hit_DetY_truth, Hit_nOp_indv, Hit_nOp_gen_indv;
    Hit_Edep_reco.resize(2);
    Hit_Edep_truth.resize(2);
    Hit_Scat_elecE.resize(2);
    Hit_DetX_reco.resize(2);
    Hit_DetY_reco.resize(2);
    Hit_DetX_truth.resize(2);
    Hit_DetY_truth.resize(2);
    Hit_nOp_indv.resize(2);
    Hit_nOp_gen_indv.resize(2);
     
    for (int mdl =0; mdl <2; mdl++){
      for (int i =0; i<8; i++){
	for (int j =0; j<8; j++){
	  nOp_crys[mdl] = Modules[mdl][i][j][nOpPho];
	  nOp_gen_crys[mdl] = Modules[mdl][i][j][nGenOp];
	  EDep_crys[mdl] = Modules[mdl][i][j][Edep];
	  DetX[mdl] = Modules[mdl][i][j][DetPosX];
	  DetY[mdl] = Modules[mdl][i][j][DetPosY];


	  
	  Edep_reco_crys[mdl] = (1/21926.3) * nOp_crys[mdl];     //from optical photons at rear surface

	  // sigma_spread = sqrt(Edep_reco_crys[mdl]/0.511)*0.1;
	  // if(EDep_crys[mdl] > 0.) Edep_reco_crys[mdl] = rng.Gaus(EDep_crys[mdl], sigma_spread);     //spreading around the true edep
	  // else Edep_reco_crys[mdl] =0;
	  
	  E_scatElec[mdl] = Modules[mdl][i][j][E_elec];
	  
	  if(nOp_crys[mdl] != 0) {
	  //if(Edep_reco_crys[mdl] > 0.){
	    h_DetX_crys[mdl]->Fill(DetX[mdl]);
	    h_nOpPho_crys_inc[mdl]->Fill(nOp_crys[mdl]);
	    h_nOpPho_gen_crys_inc[mdl]->Fill(nOp_gen_crys[mdl]);
	      
	    h_nOp_vs_Edep_crys[mdl]->Fill(EDep_crys[mdl], nOp_crys[mdl]);
	    //Edep_reco_crys[mdl] = (1/21988.2) * nOp_crys[mdl];
	    //Edep_reco_crys[mdl] = (1/21926.3) * nOp_crys[mdl];
	    

	    h_Edep_truth_crys[mdl]->Fill(EDep_crys[mdl]);
	    h_Edep_reco_crys[mdl]->Fill(Edep_reco_crys[mdl]);
	    h_Edep_true_Vs_reco_crys[mdl]->Fill(EDep_crys[mdl], Edep_reco_crys[mdl]);

	    h_test2d->Fill(EDep_crys[mdl], EDep_crys[mdl]);
	    }
			  
	  nOp_mdl[mdl] += nOp_crys[mdl];
	  nOp_gen_mdl[mdl] += nOp_gen_crys[mdl];
	    
	  //EDep_truth_mdl[mdl] += EDep_crys[mdl];
	  EDep_reco_mdl[mdl] += Edep_reco_crys[mdl];	  

	  //-------------------getting the relavant hits using true Edep---------------------
	  if(EDep_crys[mdl] > 0.06 && EDep_crys[mdl] < 0.405){	    
	  //if(EDep_crys[mdl] > 0.06 && EDep_crys[mdl]< 0.340){
	  //if(EDep_crys[mdl] > 0.02 && EDep_crys[mdl]< 0.451){
	  //if(EDep_crys[mdl] > 0.06){
	    EDep_truth_mdl[mdl] += EDep_crys[mdl];
	    nHits_mdl_truth[mdl] ++;
	    Hit_Edep_truth[mdl].push_back(EDep_crys[mdl]);
	    Hit_Scat_elecE[mdl].push_back(E_scatElec[mdl]);
	      
	    Hit_DetX_truth[mdl].push_back(DetX[mdl]);
	    Hit_DetY_truth[mdl].push_back(DetY[mdl]);

	    h_elecE_vs_Edep_crys[mdl]->Fill(EDep_crys[mdl], E_scatElec[mdl]);
	    //h_Edep_truth_crys[mdl]->Fill(EDep_crys[mdl]);
	    h_DetXvsDetY_truth[mdl]->Fill(DetX[mdl], DetY[mdl]);

	    Hit_nOp_indv[mdl].push_back(nOp_crys[mdl]);
	    Hit_nOp_gen_indv[mdl].push_back(nOp_gen_crys[mdl]);
	      
	    h_nOpPho_crys_indv[mdl]->Fill(nOp_crys[mdl]);
	    h_nOpPho_gen_crys_indv[mdl]->Fill(nOp_gen_crys[mdl]);	   
	  }

	  //-------------------getting the relavant hits using nOpPho at rear surface---------------------
	  if(Edep_reco_crys[mdl] > 0.059 && Edep_reco_crys[mdl]< 0.407){
	    nHits_mdl_reco[mdl] ++;
	    Hit_Edep_reco[mdl].push_back(Edep_reco_crys[mdl]);	        

	    Hit_DetX_reco[mdl].push_back(DetX[mdl]);
	    Hit_DetY_reco[mdl].push_back(DetY[mdl]);
	    
	    h_DetXvsDetY_reco[mdl]->Fill(DetX[mdl], DetY[mdl]);	  	      
	  }
	  
	}  // j loop end
      }    // i loop end
      //---------------------------------------------------------------------------------------------------
      //truth hists
      h_Edep_truth_mdl[mdl] ->Fill(EDep_truth_mdl[mdl]);
      h_nHits_truth_mdl[mdl] ->Fill(nHits_mdl_truth[mdl]);

      //reco hists
      h_nOpPho_mdl_inc[mdl]->Fill(nOp_mdl[mdl]);
      h_nOpPho_gen_mdl_inc[mdl]->Fill(nOp_gen_mdl[mdl]);
      h_Edep_reco_mdl[mdl] ->Fill(EDep_reco_mdl[mdl]);
      h_nHits_reco_mdl[mdl]->Fill(nHits_mdl_reco[mdl]);

      //for calibration (using both truth and reco info at the same time)
      h_nOp_vs_Edep_mdl[mdl]->Fill(EDep_truth_mdl[mdl], nOp_mdl[mdl]);        
    }  // end mdl loop

    //cout << "working till here!!" << endl;

    if (nHits_mdl_truth[r]==2){
      if(Hit_Edep_truth[r][0]+Hit_Edep_truth[r][1] <0.512 && Hit_Edep_truth[r][0]+Hit_Edep_truth[r][1] >=0.510){
      right_Hit_Edep_vec.push_back(Hit_Edep_truth[r]);
      right_Hit_DetX_vec.push_back(Hit_DetX_truth[r]);
      right_Hit_DetY_vec.push_back(Hit_DetY_truth[r]);
      }
    }
    
    else if (nHits_mdl_truth[l]==2){
      if(Hit_Edep_truth[l][0]+Hit_Edep_truth[l][1] <0.512 && Hit_Edep_truth[l][0]+Hit_Edep_truth[l][1] >=0.510){
      left_Hit_Edep_vec.push_back(Hit_Edep_truth[l]);
      left_Hit_DetX_vec.push_back(Hit_DetX_truth[l]);
      left_Hit_DetY_vec.push_back(Hit_DetY_truth[l]);
      }
    }

    //--------------some Optical photon distributions
       //-------------for getting the sigma in the sum case

    // if(nHits_mdl_truth[r]==1){
    //   h_nOpPho_mdl_sum[r]->Fill(Hit_nOp_indv[r][0]);
    //   h_nOpPho_gen_mdl_sum[r]->Fill(Hit_nOp_gen_indv[r][0]);
    // }

    // if(nHits_mdl_truth[r]==2 && Hit_Edep_truth[r][0]+Hit_Edep_truth[r][1] <0.5111 && Hit_Edep_truth[r][0]+Hit_Edep_truth[r][1] >0.5109){
    //   h_nOpPho_mdl_sum[r]->Fill(Hit_nOp_indv[r][0]+Hit_nOp_indv[r][1]);
    //   h_nOpPho_gen_mdl_sum[r]->Fill(Hit_nOp_gen_indv[r][0]+Hit_nOp_gen_indv[r][1]);
    // }
    
    // if(nHits_mdl_truth[r]==2){
    //   h_nOpPho_mdl_sum[r]->Fill(Hit_nOp_indv[r][0]+Hit_nOp_indv[r][1]);
    //   h_nOpPho_gen_mdl_sum[r]->Fill(Hit_nOp_gen_indv[r][0]+Hit_nOp_gen_indv[r][1]);
    // }

    // if(nHits_mdl_truth[r].size()>=1){
    //    if(Hit_Edep_truth[r][0]+Hit_Edep_truth[r][1] <0.5111 && Hit_Edep_truth[r][0]+Hit_Edep_truth[r][1] >=0.5109){
    // 	 h_Edep_sum_reco->Fill(Hit_Edep_reco[r][0]+Hit_Edep_reco[r][1]);
    //    }
    //   //cout << Hit_Edep_reco[r][0]+Hit_Edep_reco[r][1] << endl;
    // }

    if (nHits_mdl_truth[r]>0) Trig ++;
    
    if(EDep_truth_mdl[r] > 0.5109 && EDep_truth_mdl[r] < 0.5111) {
      full_Edep ++;

      if (nHits_mdl_truth[r]==1) sing_pix++;
      else if (nHits_mdl_truth[r]==2) two_pix++;
      else if (nHits_mdl_truth[r]>2) mult_pix++;
    }
    

    



    //===============================Calculating theta using truth Edep====================================
    if(nHits_mdl_truth[r]==2 && nHits_mdl_truth[l]==2){
      if(Hit_Edep_truth[r][0]+Hit_Edep_truth[r][1] <= 0.512 && Hit_Edep_truth[r][0]+Hit_Edep_truth[r][1] >=0.510){	
	if(Hit_Edep_truth[l][0]+Hit_Edep_truth[l][1] <0.512 && Hit_Edep_truth[l][0]+Hit_Edep_truth[l][1] >=0.510){
	  check_true++;
	         	
	  check_top++;
	  float scat_Edep[2] ={0}, abs_Edep[2] ={0}, theta[2] = {9999};
	  int scat_idx[2] ={99}, abs_idx[2] ={99};
	  double phi[2] ={999}, dPix[2] = {99999};
	    
	  for (int Mdl =0; Mdl < 2; Mdl ++){
	    if (Hit_Edep_truth[Mdl][0] < Hit_Edep_truth[Mdl][1]){	  	   
	      scat_idx[Mdl] = 0;
	      abs_idx[Mdl] = 1;
	    }
	      
	    else if (Hit_Edep_truth[Mdl][1] < Hit_Edep_truth[Mdl][0]){
	      scat_idx[Mdl] = 1;
	      abs_idx[Mdl] = 0;	    
	    }
	      
	    if(scat_idx[Mdl] != 99 || abs_idx[Mdl] != 99){
	      scat_Edep[Mdl] = Hit_Edep_truth[Mdl][scat_idx[Mdl]];
	      abs_Edep[Mdl] = Hit_Edep_truth[Mdl][abs_idx[Mdl]];
		
	      h_scatElec_E[Mdl]->Fill(Hit_Scat_elecE[Mdl][scat_idx[Mdl]]);	  

	      
	      theta[Mdl] = GetScatAngle(scat_Edep[Mdl]);		

	      h_theta_truth[Mdl]->Fill(Deg(theta[Mdl]));	     
	      h_elecE_vs_scatEdep[Mdl]->Fill(scat_Edep[Mdl], Hit_Scat_elecE[Mdl][scat_idx[Mdl]]);
			       		
	      h_Edep_scat_truth[Mdl]->Fill(scat_Edep[Mdl]);
	      h_Edep_abs_truth[Mdl]->Fill(abs_Edep[Mdl]);
	      h_Edep_scatVsabs_truth[Mdl]->Fill(scat_Edep[Mdl], abs_Edep[Mdl]);  	  	  
		
	      //------------------Calculating Phi---------------------------------------
	      float abs_DetX = Hit_DetX_truth[Mdl][abs_idx[Mdl]];
	      float abs_DetY = Hit_DetY_truth[Mdl][abs_idx[Mdl]];
		
	      float scat_DetX = Hit_DetX_truth[Mdl][scat_idx[Mdl]];
	      float scat_DetY = Hit_DetY_truth[Mdl][scat_idx[Mdl]];	      
			      
	      phi[Mdl] = atan2(abs_DetY-scat_DetY, abs_DetX-scat_DetX);
	      
	      h_phi_truth[Mdl]->Fill(Deg(phi[Mdl]));
	      
		
	      //h_phi_truth[Mdl]->Fill(phi_truth[Mdl]*180/ TMath::Pi());
	      dPix[Mdl]= sqrt(((abs_DetX-scat_DetX)*(abs_DetX-scat_DetX))+((abs_DetY-scat_DetY)*(abs_DetY-scat_DetY)));
		
	      //cout << "pix distance: " << Pix_distance[Mdl] << endl;
	      h_dPix_truth[Mdl]->Fill(dPix[Mdl]);


	      phi_vec[Mdl].push_back(phi[Mdl]);
	      theta_vec[Mdl].push_back(theta[Mdl]);	      	      
	    }
	  }
	  	  	  
	  h_Phi1_vs_Phi2_truth->Fill(Deg(phi[r]), Deg(phi[l]));
	  
	  float dPhi = phi[r] - phi[l];

	  if(dPhi > Pi) dPhi -= 2*Pi;
	  else if (dPhi < -Pi) dPhi += 2*Pi;
	  
	  //Correlated dPhi distribution for theta inc-------------------
	  h_dPhi_truth_inc->Fill(Deg(dPhi));

	  if (dPix[r] > 3.2 && dPix[l] > 3.2) h_dPhi_truth_inc_dexc->Fill(Deg(dPhi));

	  //Correlated dPhi distribution for theta 70 to 90 degs-------------------
	  if (theta[r]< Rad(90) && theta[r]> Rad(72) && theta[l]< Rad(90) && theta[l]> Rad(72)){	    
	    //h_phi_truth[r]->Fill(Deg(phi[r]));
	    //h_phi_truth[l]->Fill(Deg(phi[l]));
	    //h_Phi1_vs_Phi2_truth->Fill(Deg(phi[r]), Deg(phi[l]));	    	  
	    h_dPhi_truth_70_90->Fill(Deg(dPhi));

	    if (dPix[r] > 3.2 && dPix[l] > 3.2) h_dPhi_truth_70_90_dexc->Fill(Deg(dPhi));
	    
	  }	  
	}
      }
    }

    //##############################################################################################################

        //===============================Calculating theta using Reco Edep====================================
    if(nHits_mdl_reco[r]==2 && nHits_mdl_reco[l]==2){            
      if(Hit_Edep_reco[r][0]+Hit_Edep_reco[r][1] < 0.55 && Hit_Edep_reco[r][0]+Hit_Edep_reco[r][1] >=0.493){ 
    	if(Hit_Edep_reco[l][0]+Hit_Edep_reco[l][1] < 0.55 && Hit_Edep_reco[l][0]+Hit_Edep_reco[l][1] >=0.493){
	  check_reco ++;
    	  float scat_Edep[2] ={0.}, abs_Edep[2] ={0.}, theta[2] = {9999.};
    	  int scat_idx[2] ={99}, abs_idx[2] ={99};
    	  double phi[2] ={999.}, dPix[2] = {99999.};
	    
    	  for (int Mdl =0; Mdl < 2; Mdl ++){
    	    if (Hit_Edep_reco[Mdl][0] < Hit_Edep_reco[Mdl][1]){	  	   
    	      scat_idx[Mdl] = 0;
    	      abs_idx[Mdl] = 1;
    	    }
	      
    	    else if (Hit_Edep_reco[Mdl][1] < Hit_Edep_reco[Mdl][0]){
    	      scat_idx[Mdl] = 1;
    	      abs_idx[Mdl] = 0;	    
    	    }
	      
    	    if(scat_idx[Mdl] != 99 || abs_idx[Mdl] != 99){
    	      scat_Edep[Mdl] = Hit_Edep_reco[Mdl][scat_idx[Mdl]];
    	      abs_Edep[Mdl] = Hit_Edep_reco[Mdl][abs_idx[Mdl]];
			     		
    	      theta[Mdl] = GetScatAngle(scat_Edep[Mdl]);	  
		
    	      h_theta_reco[Mdl]->Fill(theta[Mdl]*180/ TMath::Pi());	     
				
    	      h_Edep_scat_reco[Mdl]->Fill(scat_Edep[Mdl]);
    	      h_Edep_abs_reco[Mdl]->Fill(abs_Edep[Mdl]);
    	      h_Edep_scatVsabs_reco[Mdl]->Fill(scat_Edep[Mdl], abs_Edep[Mdl]);
	    
	    	  	  	 	  
    	      //------------------Calculating Phi---------------------------------------
    	      float abs_DetX = Hit_DetX_reco[Mdl][abs_idx[Mdl]];
    	      float abs_DetY = Hit_DetY_reco[Mdl][abs_idx[Mdl]];
	  
    	      float scat_DetX = Hit_DetX_reco[Mdl][scat_idx[Mdl]];
    	      float scat_DetY = Hit_DetY_reco[Mdl][scat_idx[Mdl]];

    	      //phi_truth[Mdl] = atan((abs_DetY_truth-scat_DetY_truth)/(abs_DetX_truth-scat_DetX_truth));	
    	      phi[Mdl] = atan2(abs_DetY-scat_DetY, abs_DetX-scat_DetX);
	  
    	      //h_phi[Mdl]->Fill(phi[Mdl]*180/ TMath::Pi());
    	      dPix[Mdl]= sqrt(((abs_DetX-scat_DetX)*(abs_DetX-scat_DetX))+((abs_DetY-scat_DetY)*(abs_DetY-scat_DetY)));

    	      //cout << "pix distance: " << Pix_distance[Mdl] << endl;
    	      h_dPix_reco[Mdl]->Fill(dPix[Mdl]);	      
    	    }
    	  }	  
	      	  
    	    h_phi_reco[r]->Fill(Deg(phi[r]));
    	    h_phi_reco[l]->Fill(Deg(phi[l]));
    	    h_Phi1_vs_Phi2_reco->Fill(Deg(phi[r]), Deg(phi[l]));
	    
    	    float dPhi = phi[r] - phi[l];

    	    if(dPhi > Pi) dPhi -= 2*Pi;
	    else if (dPhi < -Pi) dPhi += 2*Pi;
    	    
	    h_dPhi_reco_inc->Fill(Deg(dPhi));
    	    
	    if (theta[r]< Rad(90) && theta[r]> Rad(72) && theta[l]< Rad(90) && theta[l]> Rad(72)){
	      h_dPhi_reco_70_90->Fill(Deg(dPhi));	    
	    }
	
    	  //}
    	}
      }
    }


    
  } //jentry loop end

  // cout << "\ncheck true: " << check_true << endl;
  // cout << "check_reco: " << check_reco << endl;

  cout << "\n Total Triggered: " << Trig << endl;
  cout << "full Edep: " << 1.*full_Edep/Trig << endl;
  cout << "1 pix fired " << 1.*sing_pix/Trig << endl;
  cout << "2 pix fired " << 1.*two_pix/Trig << endl;
  cout << "more than fired " << 1.*mult_pix/Trig << endl;
  
  //===================================making a proxy to uncorrelated by shuffling=====================================================
  std::vector<size_t> idx_r(right_Hit_Edep_vec.size());
  std::iota(idx_r.begin(), idx_r.end(), 0);
  
  std::mt19937 rng_r(12345);
  std::shuffle(idx_r.begin(), idx_r.end(), rng_r);

  
  std::vector<size_t> idx_l(left_Hit_Edep_vec.size());
  std::iota(idx_l.begin(), idx_l.end(), 0);
  
  std::mt19937 rng_l(54321);
  std::shuffle(idx_l.begin(), idx_l.end(), rng_l);

  //ucomment from here  
  // std::vector<float> right_phi_shuff, right_theta_shuff, left_phi_shuff, left_theta_shuff;

  // for (auto i : idx_r){
  //   right_phi_shuff.push_back(phi_vec[r][i]);
  //   right_theta_shuff.push_back(theta_vec[r][i]);
  // }

  // for (auto i : idx_l){
  //   left_phi_shuff.push_back(phi_vec[l][i]);
  //   left_theta_shuff.push_back(theta_vec[l][i]);
  // }
  
  
  // float phi1 = -99999, phi2 = -9999, dPhi_mix = -9999;
  
  // // for (int i =0; i< right_phi_shuff.size()-1; i++){
  // //   if (Deg(right_theta_shuff[i])>70 && Deg(right_theta_shuff[i])<90 && Deg(right_theta_shuff[i+1])>70 && Deg(right_theta_shuff[i+1])<90){      
  // //   phi1 = right_phi_shuff[i];

  // //   h_phi_truth_EvtMix[r]->Fill(Deg(phi1));
    
  // //   phi2 = right_phi_shuff[i+1];
  // //   float dPhi_mix = phi1 - phi2;
    
  // //   if (abs(dPhi_mix) <= Pi) h_dPhi_truth_EvtMix->Fill(dPhi_mix*180/ Pi);
  // //   else if(abs(dPhi_mix) > Pi && dPhi_mix>0) h_dPhi_truth_EvtMix->Fill((dPhi_mix-(2*Pi))*180/ Pi);
  // //   else if(abs(dPhi_mix) > Pi && dPhi_mix<0) h_dPhi_truth_EvtMix->Fill((dPhi_mix+(2*Pi))*180/ Pi);
  // //   }
  // // }

  // for (int i =0; i< right_phi_shuff.size(); i++){
  //   phi1 = right_phi_shuff[i];
  //   phi2 = left_phi_shuff[i];
  //   //phi2 = right_phi_shuff[i+1];

  //   h_phi_truth_EvtMix[r]->Fill(Deg(phi1));
    
  //   dPhi_mix = phi1 - phi2;
    
  //   //dPhi distribuitions for theta inc----------------
  //   if (abs(dPhi_mix) <= Pi) h_dPhi_truth_EvtMix_inc->Fill(Deg(dPhi_mix));
  //   else if(abs(dPhi_mix) > Pi && dPhi_mix>0) h_dPhi_truth_EvtMix_inc->Fill(Deg(dPhi_mix-(2*Pi)));
  //   else if(abs(dPhi_mix) > Pi && dPhi_mix<0) h_dPhi_truth_EvtMix_inc->Fill(Deg(dPhi_mix+(2*Pi)));

    
  //   //dPhi for theta 70 to 90 degs--------------------
  //    if (Deg(right_theta_shuff[i]) >72 && Deg(right_theta_shuff[i])<90 && Deg(left_theta_shuff[i])>72 && Deg(left_theta_shuff[i])<90){

  //   // if (Deg(right_theta_shuff[i]) >72 && Deg(right_theta_shuff[i])<90 && Deg(right_theta_shuff[i+1])>72 && Deg(right_theta_shuff[i+1])<90){
      
  //     if (abs(dPhi_mix) <= Pi) h_dPhi_truth_EvtMix_70_90->Fill(Deg(dPhi_mix));
  //     else if(abs(dPhi_mix) > Pi && dPhi_mix>0) h_dPhi_truth_EvtMix_70_90->Fill(Deg(dPhi_mix-(2*Pi)));
  //     else if(abs(dPhi_mix) > Pi && dPhi_mix<0) h_dPhi_truth_EvtMix_70_90->Fill(Deg(dPhi_mix+(2*Pi)));          }
  // }
    
  // TH1F* h_dPhi_tmp_corr_70_90 = (TH1F*) h_dPhi_truth_70_90->Clone("tmp2");
  // h_dPhi_tmp_corr_70_90->SetDirectory(nullptr);

  // TH1F* h_dPhi_tmp_EvtMix_70_90 = (TH1F*) h_dPhi_truth_EvtMix_70_90->Clone("tmp3");
  // h_dPhi_tmp_EvtMix_70_90->SetDirectory(nullptr); 
  

  // h_dPhi_tmp_corr_70_90->Rebin(5);
  // h_dPhi_tmp_EvtMix_70_90->Rebin(5);
 
  // h_dPhi_tmp_corr_70_90->Divide(h_dPhi_tmp_EvtMix_70_90);
  
  // h_dPhi_truth_corrected_70_90->Add(h_dPhi_tmp_corr_70_90);

  // cout << "No. of events with exactly 2 Hits: " << nEvts2Hits << endl;

  // cout << "\nTheta: " << Pi << endl;
  // cout << "\nScat ElecE: " << GetScatElecE(90) << endl;

  //------------- A different approach-----------------------------------
  
  for (auto i : idx_r){
    right_Hit_Edep_shuff.push_back(right_Hit_Edep_vec[i]);
    right_Hit_DetX_shuff.push_back(right_Hit_DetX_vec[i]);
    right_Hit_DetY_shuff.push_back(right_Hit_DetY_vec[i]);
  }
  
  for (auto i : idx_l){
    // left_Hit_Edep_shuff.push_back(left_Hit_Edep_vec[i]);
    // left_Hit_DetX_shuff.push_back(left_Hit_DetX_vec[i]);
    // left_Hit_DetY_shuff.push_back(left_Hit_DetY_vec[i]);

    left_Hit_Edep_shuff.push_back(right_Hit_Edep_vec[i]);
    left_Hit_DetX_shuff.push_back(right_Hit_DetX_vec[i]);
    left_Hit_DetY_shuff.push_back(right_Hit_DetY_vec[i]);

  }

  int nEvtMix = 0;

  if (right_Hit_Edep_shuff.size() < left_Hit_Edep_shuff.size()) nEvtMix = right_Hit_Edep_shuff.size();
  else  nEvtMix = left_Hit_Edep_shuff.size();

  vector<float> right_Edep_vec, right_DetX_vec, right_DetY_vec, left_Edep_vec, left_DetX_vec, left_DetY_vec;
  for (int Evt = 0; Evt < check_top; Evt++){
    right_Edep_vec = right_Hit_Edep_shuff[Evt];
    right_DetX_vec = right_Hit_DetX_shuff[Evt];
    right_DetY_vec = right_Hit_DetY_shuff[Evt];

    left_Edep_vec = left_Hit_Edep_shuff[Evt];
    left_DetX_vec = left_Hit_DetX_shuff[Evt];
    left_DetY_vec = left_Hit_DetY_shuff[Evt];

    //if (right_Edep_vec.size() ==2 && left_Edep_vec.size() ==2){
    //if(right_Edep_vec[0]+right_Edep_vec[1] <0.512 && right_Edep_vec[0]+right_Edep_vec[1] >=0.510){
    //if(left_Edep_vec[0]+left_Edep_vec[1] <0.512 && left_Edep_vec[0]+left_Edep_vec[1] >=0.510){
  	  check_bottom ++;
  	  float scat_Edep[2] ={0.}, abs_Edep[2] ={0.}, theta[2] = {9999};
  	  int scat_idx[2] ={99}, abs_idx[2] ={99};
  	  double phi[2] ={999}, dPix[2] = {99999};
  	  float scat_DetX[2] = {9999.}, abs_DetX[2] = {9999.}, scat_DetY[2] = {9999.}, abs_DetY[2] = {9999.};
	  
  	  if (right_Edep_vec[0] < right_Edep_vec[1]){	  	   
  	    scat_idx[r] = 0;
  	    abs_idx[r] = 1;
  	  }
	  
  	  else if (right_Edep_vec[0] > right_Edep_vec[1]){	  	   
  	    scat_idx[r] = 1;
  	    abs_idx[r] = 0;
  	  }
	  
  	  //---------
  	  if (left_Edep_vec[0] < left_Edep_vec[1]){	  	   
  	    scat_idx[l] = 0;
  	    abs_idx[l] = 1;
  	  }
	  
  	  else if (left_Edep_vec[0] > left_Edep_vec[1]){	  	   
  	    scat_idx[l] = 1;
  	    abs_idx[l] = 0;
  	  }
	  
  	  scat_Edep[r] = right_Edep_vec[scat_idx[r]];
  	  abs_Edep[r] = right_Edep_vec[abs_idx[r]];
	  
  	  scat_Edep[l] = left_Edep_vec[scat_idx[l]];
  	  abs_Edep[l] = left_Edep_vec[abs_idx[l]];
	  
	  
  	  //calculating theta-------------------------------	  	      
  	  theta[r] = GetScatAngle(scat_Edep[r]);
  	  theta[l] = GetScatAngle(scat_Edep[l]);
	  
  	  //------------------Calculating Phi---------------------------------------
  	  abs_DetX[r] = right_DetX_vec[abs_idx[r]];
  	  abs_DetY[r] = right_DetY_vec[abs_idx[r]];
	  
  	  scat_DetX[r] = right_DetX_vec[scat_idx[r]];
  	  scat_DetY[r] = right_DetY_vec[scat_idx[r]];
	  
  	  //---
  	  abs_DetX[l] = left_DetX_vec[abs_idx[l]];
  	  abs_DetY[l] = left_DetY_vec[abs_idx[l]];
	  
  	  scat_DetX[l] = left_DetX_vec[scat_idx[l]];
  	  scat_DetY[l] = left_DetY_vec[scat_idx[l]];
	  
	  
  	  phi[r] = atan2(abs_DetY[r]-scat_DetY[r], abs_DetX[r]-scat_DetX[r]);
  	  phi[l] = atan2(abs_DetY[l]-scat_DetY[l], abs_DetX[l]-scat_DetX[l]);	  
	  
  	  float dPhi_mix = phi[r] - phi[l];
	  
  	  h_phi_truth_EvtMix[r]->Fill(Deg(phi[r]));
  	  h_phi_truth_EvtMix[l]->Fill(Deg(phi[l]));
	 	        
  	  //h_test2d->Fill(Deg(phi[r]),Deg(phi[r]));
	  
  	  //dPhi distribuitions for theta inc----------------
  	  if (abs(dPhi_mix) <= Pi) h_dPhi_truth_EvtMix_inc->Fill(Deg(dPhi_mix));
  	  else if(dPhi_mix > Pi) h_dPhi_truth_EvtMix_inc->Fill(Deg(dPhi_mix-(2*Pi)));
  	  else if(dPhi_mix < -Pi) h_dPhi_truth_EvtMix_inc->Fill(Deg(dPhi_mix+(2*Pi)));

	  //dPhi for theta 70 to 90 degs--------------------
	  if (Deg(theta[r]) >=72 && Deg(theta[r]) <=90 && Deg(theta[l]) >=72 && Deg(theta[l]) <=90){
	    if (abs(dPhi_mix) <= Pi) h_dPhi_truth_EvtMix_70_90->Fill(Deg(dPhi_mix));
	    else if(dPhi_mix > Pi) h_dPhi_truth_EvtMix_70_90->Fill(Deg(dPhi_mix-(2*Pi)));
	    else if(dPhi_mix < -Pi) h_dPhi_truth_EvtMix_70_90->Fill(Deg(dPhi_mix+(2*Pi)));
	  }
  
  }
  // }      
      //  }
    
//}

  TH1F* h_dPhi_tmp_corr_inc = (TH1F*) h_dPhi_truth_inc->Clone("tmp1_inc");
  h_dPhi_tmp_corr_inc->SetDirectory(nullptr);
  
  TH1F* h_dPhi_tmp_EvtMix_inc = (TH1F*) h_dPhi_truth_EvtMix_inc->Clone("tmp2_inc");
  h_dPhi_tmp_EvtMix_inc->SetDirectory(nullptr); 
  
  h_dPhi_tmp_corr_inc->Rebin(5);
  h_dPhi_tmp_EvtMix_inc->Rebin(5);
 
  h_dPhi_tmp_corr_inc->Divide(h_dPhi_tmp_EvtMix_inc);
  
  h_dPhi_truth_corrected_inc->Add(h_dPhi_tmp_corr_inc);

  //---------
   TH1F* h_dPhi_tmp_corr_70_90 = (TH1F*) h_dPhi_truth_70_90->Clone("tmp1_70_90");
  h_dPhi_tmp_corr_70_90->SetDirectory(nullptr);

  TH1F* h_dPhi_tmp_EvtMix_70_90 = (TH1F*) h_dPhi_truth_EvtMix_70_90->Clone("tmp2_70_90");
  h_dPhi_tmp_EvtMix_70_90->SetDirectory(nullptr); 
  

  h_dPhi_tmp_corr_70_90->Rebin(5);
  h_dPhi_tmp_EvtMix_70_90->Rebin(5);
 
  h_dPhi_tmp_corr_70_90->Divide(h_dPhi_tmp_EvtMix_70_90);
  
  h_dPhi_truth_corrected_70_90->Add(h_dPhi_tmp_corr_70_90);

  // cout << "check top: " << check_top << endl;
  // cout << "check bottom: " << check_bottom << endl;    
}

void AnalyzeLightBSM::FillHistogram(double time, double posX, double posY, double posZ, double theta_deg, double eta_deg, double Ein, double Eout){
  h_hitTime->Fill(time);
  h_hitPosX->Fill(posX);
  h_hitPosY->Fill(posY);
  h_hitPosZ->Fill(posZ);

  h_theta->Fill(theta_deg);

  h_eta->Fill(eta_deg);	    

  h_Ein->Fill(Ein);
  h_Eout->Fill(Eout);
}

string AnalyzeLightBSM::MatName(int DetId){
  string mat = " ";
  if (DetId < 0) mat = "error";
  if (DetId /4000 >=1) mat = "BGO"; 
  else if (DetId/ 1000 >=1 && DetId/ 4000 <1) mat = "Scint";

  return mat;
}

float AnalyzeLightBSM::GetScatAngle(float E_dep){
  double scatAngle = acos(1-(0.511*((1/(0.511-E_dep)) - (1/0.511))));
  return scatAngle;
}

float AnalyzeLightBSM::GetScatElecE(float theta){
  double elec_E = 0.511 -(0.511/2-cos(theta*Pi/180));
  return elec_E;
}
 
float AnalyzeLightBSM::Deg(float Rad){
  double degree = Rad * 180 /Pi;
  return degree;
}

float AnalyzeLightBSM::Rad(float Deg){
  double radian = Deg * Pi/180;
  return radian;
}
 



















//================================================================================================================================================================================================================












//some baks


// if (jentry == 195) {

//   cout << "\ntime(" << HitTime_sort[iHit] << "): " << HitPosX_sort[iHit] << ", " << HitPosY_sort[iHit] << ", " << HitPosZ_sort[iHit] << endl;
//   // cout << "detID" << (*Hit_DetId)[iHit] << endl;
//   cout << "Ana_theta: " << scat_theta*180.0 / TMath::Pi() << endl;
//   cout << "step_theta: " << scat_theta_SD*180.0 / TMath::Pi() << endl;
//   //cout << "v1: " << v1.X() << << endl;
//   //cout << "v2:" << v2 << endl;
//   cout << "v1: (" << v1.X() << ", " << v1.Y() << ", " << v1.Z() << ")" << endl;
//   cout << "v2: (" << v2.X() << ", " << v2.Y() << ", " << v2.Z() << ")" << endl;
//   cout << "vin: (" << vin.X() << ", " << vin.Y() << ", " << vin.Z() << ")" << endl;
//   cout << "vout: (" << vout.X() << ", " << vout.Y() << ", " << vout.Z() << ")" << "\n\n";
//   //cout << "Hit position Ana: (" << (*Hit_PositionX)[iHit] << ", " << (*Hit_PositionY)[iHit] << ", " << (*Hit_PositionZ)[iHit] << ")" << "\n\n";	
// }
	 	    
// cout << "scattering angle SD: " << scat_theta*180.0 / TMath::Pi() << endl;      
// cout << "scattering angle Ana: " << scat_theta_SD*180.0 / TMath::Pi() << "\n\n";



// ROOT::Math::XYZVector Photon1_Pol(Relv4Hits[A1][Pol_X], Relv4Hits[A1][Pol_Y], Relv4Hits[A1][Pol_Z]);
// ROOT::Math::XYZVector Photon2_Pol(Relv4Hits[A2][Pol_X], Relv4Hits[A2][Pol_Y], Relv4Hits[A2][Pol_Z]);

// double Angle_Pol1Pol2 = acos(Photon1_Pol.Unit().Dot(Photon2_Pol.Unit())) *180.0 / TMath::Pi();
// h_rough[0]->Fill(Angle_Pol1Pol2);
	  

















































//BACKUP

//***************************************
//convention followed for nHits ==2:
//S: Scintillator, B: BGO
//left for location of 1st hit;
//right for 2nd;
//1 and 2 in S or B for identifying same or different crystal of a material;
//no identifier when hits are in different materials;
//e.g. nHit2_S1S1 means 1st hit in scintillator 2nd hit in the same scinitillator
//***************************************
	  	  
// enum HitsLoc{all_inc =0, nHit0 =1, nHit1_S =2, nHit1_B =3, nHit2_S_AnyDiffDet =4, nHit2_S1S1 =5, nHit2_S1S2 =6, nHit2_SB =7, nHit2_B1B1 =8, nHit2_B1B2 =9, nHit2_BS =10, remain =11, test = 12};

// //Definitng the booleans for each category
// bool inc_all =false, noHit =false;
// bool hit1_S =false, hit1_B =false;
// bool hit2_S_AnyDiffDet = false, hit2_S1S1 =false, hit2_S1S2 =false, hit2_SB =false, hit2_B1B1 =false, hit2_B1B2 =false, hit2_BS =false, rest =false;
// //bool teSt = false;
// bool teSt = true;
	  	  
// int DetId1 = HitDetId_sort[0];
// int DetId2 = -1;
// if(nHitComp >=2) DetId2 = HitDetId_sort[1];

// if (nHitComp >0) inc_all = true;
	  
// //if (nHitComp > 0 && MatName(DetId1) == "BGO") {hit_BGO_inc = true; DetId1 = HitDetId_sort[0];}
// // if (nHitComp >=2 && DetId1 != DetId2 && MatName(DetId1) == "Scint") {
// //   teSt = true;
// // }
	  
// if (nHitComp ==0) noHit = true;
	  
// else if (nHitComp ==1){
//   n1++;
//   if (MatName(DetId1) =="Scint") {n1s++; hit1_S = true;} 
//   else if (MatName(DetId1) =="BGO") {n1b++; hit1_B = true;}
//   else rest = true;
// }

// // else if (nHitComp ==2){
// else if (nHitComp >=2){	  	    
//    if (DetId1 == DetId2){	   		    
//     if (MatName(DetId1) =="Scint") hit2_S1S1 = true; 
//     else if (MatName(DetId1) =="BGO") hit2_B1B1 = true;
//     else rest = true;	      
//   }	    
	    
//    else if (DetId1 != DetId2)
//      {if(MatName(DetId1) =="Scint") hit2_S_AnyDiffDet = true;
// 	 if(MatName(DetId1) =="Scint" && MatName(DetId2) =="Scint") {hit2_S1S2 = true;} //cout << "event with s1s2: " << jentry +1 << endl;}
// 	 else if(MatName(DetId1) =="BGO" && MatName(DetId2) =="BGO") hit2_B1B2 = true;
// 	 else if (MatName(DetId1) =="Scint" && MatName(DetId2) =="BGO") {hit2_SB = true;}// cout << "event with SB: " << jentry +1 << endl;}	       
// 	 else if (MatName(DetId1) =="BGO" && MatName(DetId2) == "Scint") hit2_BS = true;
// 	 else rest = true;
//      }
// }
	  
// else rest = true;

// //if (hit2_S1S1) cout << "two hits in same scint" << endl;

// //else if (DetId1 != DetId2 && MatName(DetId1) =="Scint" && MatName(DetId2) =="Scint") hit2_S1S2 == true;
// //else if (MatName(DetId1) =="Scint" && MatName(DetId2) =="BGO") hit2_SB == true;
// //else if (DetId1 == DetId2 && MatName(DetId1) =="BGO") hit2_B1B1 == true;
// //else if (DetId1 != DetId2 && MatName(DetId1) =="BGO" && MatName(DetId2) =="BGO") hit2_B1B2 == true;
// //else if (DetId1 != DetId2 && MatName(DetId1) =="BGO" && MatName(DetId2) == "Scint") hit2_BS == true;
	  

	  
	  
// if (procId !=0 && procId != 1){       //excluding rayleigh && phot (not required exactly but kept of safety purpose)
//   h_ScatAng_SDvsAna->Fill(scat_theta, scat_theta_SD);
	    
//   if (inc_all) FillHistogram(all_inc, Time, PosX, PosY, PosZ, scat_theta_SD_deg, EtaPolScat_deg, Pol0, Pol1, Pol2, Ein, Eout);
//   if (teSt && scat_theta_SD > 0) {FillHistogram(test, Time, PosX, PosY, PosZ, scat_theta_SD_deg, EtaPolScat_deg, Pol0, Pol1, Pol2, Ein, Eout);	     
//   //if (teSt) FillHistogram(test, Time, PosX, PosY, PosZ, scat_theta_SD_deg, EtaPolScat_deg, Pol0, Pol1, Pol2, Ein, Eout);
//   }
//   if (noHit) FillHistogram(nHit0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	    
//   if (hit1_S) FillHistogram(nHit1_S, Time, PosX, PosY, PosZ, scat_theta_SD_deg, EtaPolScat_deg, Pol0, Pol1, Pol2, Ein, Eout);
//   if (hit1_B) FillHistogram(nHit1_B, Time, PosX, PosY, PosZ, scat_theta_SD_deg, EtaPolScat_deg, Pol0, Pol1, Pol2, Ein, Eout);
	    
//   if (hit2_S_AnyDiffDet) FillHistogram(nHit2_S_AnyDiffDet, Time, PosX, PosY, PosZ, scat_theta_SD_deg, EtaPolScat_deg, Pol0, Pol1, Pol2, Ein, Eout);	      
//   if (hit2_S1S1) FillHistogram(nHit2_S1S1, Time, PosX, PosY, PosZ, scat_theta_SD_deg, EtaPolScat_deg, Pol0, Pol1, Pol2, Ein, Eout);	      
//   if (hit2_S1S2) FillHistogram(nHit2_S1S2, Time, PosX, PosY, PosZ, scat_theta_SD_deg, EtaPolScat_deg, Pol0, Pol1, Pol2, Ein, Eout);
//   if (hit2_SB) FillHistogram(nHit2_SB, Time, PosX, PosY, PosZ, scat_theta_SD_deg, EtaPolScat_deg, Pol0, Pol1, Pol2, Ein, Eout);
	    
//   if (hit2_B1B1) FillHistogram(nHit2_B1B1, Time, PosX, PosY, PosZ, scat_theta_SD_deg, EtaPolScat_deg, Pol0, Pol1, Pol2, Ein, Eout);
//   if (hit2_B1B2) FillHistogram(nHit2_B1B2, Time, PosX, PosY, PosZ, scat_theta_SD_deg, EtaPolScat_deg, Pol0, Pol1, Pol2, Ein, Eout);
//   if (hit2_BS) FillHistogram(nHit2_BS, Time, PosX, PosY, PosZ, scat_theta_SD_deg, EtaPolScat_deg, Pol0, Pol1, Pol2, Ein, Eout);
	    
//   if (rest) FillHistogram(remain, Time, PosX, PosY, PosZ, scat_theta_SD_deg, EtaPolScat_deg, Pol0, Pol1, Pol2, Ein, Eout);	    	    	    	    	  
// } //procId condition end	


// h_Ana_dt_A1S1->Fill(Ana_dt_A1S1);
// h_Ana_dt_A1S2->Fill(Ana_dt_A1S2);
// h_Ana_dt_A2S1->Fill(Ana_dt_A2S1);
// h_Ana_dt_A2S2->Fill(Ana_dt_A2S2);
// h_Ana_dt_A1A2->Fill(Ana_dt_A1A2);
// h_Ana_dt_S1S2->Fill(abs(Ana_dt_S1S2));




//############   IMPORTANT ################

    // //===============================Calculating theta using Reco Edep====================================
    // if(nHits_mdl_reco[r]==2 && nHits_mdl_reco[l]==2){      
    //   if(Hit_Edep_reco[r][0]+Hit_Edep_reco[r][1] < (0.511+(3*sigma)) && Hit_Edep_reco[r][0]+Hit_Edep_reco[r][1] >=(0.511-3*sigma)){	  
    // 	if(Hit_Edep_reco[l][0]+Hit_Edep_reco[l][1] < (0.511+(3*sigma)) && Hit_Edep_reco[l][0]+Hit_Edep_reco[l][1] >=(0.511-3*sigma)){
    // 	  float scat_Edep[2] ={0}, abs_Edep[2] ={0}, theta[2] = {9999};
    // 	  int scat_idx[2] ={99}, abs_idx[2] ={99};
    // 	  double phi[2] ={999}, dPix[2] = {99999};
	    
    // 	  for (int Mdl =0; Mdl < 2; Mdl ++){
    // 	    if (Hit_Edep_reco[Mdl][0] < Hit_Edep_reco[Mdl][1]){	  	   
    // 	      scat_idx[Mdl] = 0;
    // 	      abs_idx[Mdl] = 1;
    // 	    }
	      
    // 	    else if (Hit_Edep_reco[Mdl][1] < Hit_Edep_reco[Mdl][0]){
    // 	      scat_idx[Mdl] = 1;
    // 	      abs_idx[Mdl] = 0;	    
    // 	    }
	      
    // 	    if(scat_idx[Mdl] != 99 || abs_idx[Mdl] != 99){
    // 	      scat_Edep[Mdl] = Hit_Edep_reco[Mdl][scat_idx[Mdl]];
    // 	      abs_Edep[Mdl] = Hit_Edep_reco[Mdl][abs_idx[Mdl]];
			     		
    // 	      theta[Mdl] = GetScatAngle(scat_Edep[Mdl]);	  
		
    // 	      h_theta_reco[Mdl]->Fill(theta[Mdl]*180/ TMath::Pi());	     
				
    // 	      h_Edep_scat_reco[Mdl]->Fill(scat_Edep[Mdl]);
    // 	      h_Edep_abs_reco[Mdl]->Fill(abs_Edep[Mdl]);
    // 	      h_Edep_scatVsabs_reco[Mdl]->Fill(scat_Edep[Mdl], abs_Edep[Mdl]);
	    
	    	  	  	 	  
    // 	      //------------------Calculating Phi---------------------------------------
    // 	      float abs_DetX = Hit_DetX_reco[Mdl][abs_idx[Mdl]];
    // 	      float abs_DetY = Hit_DetY_reco[Mdl][abs_idx[Mdl]];
	  
    // 	      float scat_DetX = Hit_DetX_reco[Mdl][scat_idx[Mdl]];
    // 	      float scat_DetY = Hit_DetY_reco[Mdl][scat_idx[Mdl]];

    // 	      //phi_truth[Mdl] = atan((abs_DetY_truth-scat_DetY_truth)/(abs_DetX_truth-scat_DetX_truth));	
    // 	      phi[Mdl] = atan2(abs_DetY-scat_DetY, abs_DetX-scat_DetX);
	  
    // 	      //h_phi[Mdl]->Fill(phi[Mdl]*180/ TMath::Pi());
    // 	      dPix[Mdl]= sqrt(((abs_DetX-scat_DetX)*(abs_DetX-scat_DetX))+((abs_DetY-scat_DetY)*(abs_DetY-scat_DetY)));

    // 	      //cout << "pix distance: " << Pix_distance[Mdl] << endl;
    // 	      h_dPix_reco[Mdl]->Fill(dPix[Mdl]);	      
    // 	    }
    // 	  }	  
	  
    // 	  //if (theta[r]< Rad(90) && theta[r]> Rad(72) && theta[l]< Rad(90) && theta[l]> Rad(72)){
    // 	    h_phi_reco[r]->Fill(Deg(phi[r]));
    // 	    h_phi_reco[l]->Fill(Deg(phi[l]));
    // 	    h_Phi1_vs_Phi2_reco->Fill(Deg(phi[r]), Deg(phi[l]));
	    
    // 	    float dPhi = phi[r] - phi[l];
    // 	    if (abs(dPhi) <= Pi) h_dPhi_reco->Fill(dPhi*180/ Pi);
    // 	    else if(abs(dPhi) > Pi && dPhi>0) h_dPhi_reco->Fill((dPhi-(2*Pi))*180/ Pi);
    // 	    else if(abs(dPhi) > Pi && dPhi<0) h_dPhi_reco->Fill((dPhi+(2*Pi))*180/ Pi);	    
    // 	    //}

	
    // 	  //}
    // 	}
    //   }
    // }
