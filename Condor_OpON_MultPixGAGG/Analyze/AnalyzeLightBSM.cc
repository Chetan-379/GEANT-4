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

  int n1s = 0, n1 =0;
  int n1b =0;

  int nEvts2Hits = 0;
   
  for (Long64_t jentry=0; jentry<nentries;jentry++){
      // cout << "----------event No.: " << jentry+1 << "-----------" << endl;
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

      //sorting the hits in the ascending order of time
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      vector<double> HitPosX_sort, HitPosY_sort, HitPosZ_sort, HitTime_sort, HitScatAng_sort, HitEta_sort, // HitPol0_sort, HitPol1_sort, HitPol2_sort,
	HitEin_sort, HitEout_sort;// , HitScatMomX_sort, HitScatMomY_sort, HitScatMomZ_sort;
      vector<double> HitPosX_cpy, HitPosY_cpy, HitPosZ_cpy, HitTime_cpy, // HitPol0_cpy, HitPol1_cpy, HitPol2_cpy,
	HitScatAng_cpy, HitEta_cpy, HitEin_cpy, HitEout_cpy; // HitScatMomX_cpy, HitScatMomY_cpy, HitScatMomZ_cpy
	

      vector <int> HitProcId_sort, HitDetId_sort, HitGunId_sort;
      vector<int> HitProcId_cpy, HitDetId_cpy, HitGunId_cpy;
	   
      HitPosX_sort.clear();
      HitPosY_sort.clear();
      HitPosZ_sort.clear();
      HitTime_sort.clear();
      HitScatAng_sort.clear();
      HitEta_sort.clear();
      HitProcId_sort.clear();
      HitEout_sort.clear();
      HitEin_sort.clear();
      HitDetId_sort.clear();
      HitGunId_sort.clear();

	   
      HitPosX_cpy = *Hit_PositionX;
      HitPosY_cpy = *Hit_PositionY;
      HitPosZ_cpy = *Hit_PositionZ;
      HitTime_cpy = *Hit_Time;
      HitScatAng_cpy = *Hit_ScatAngle;
      HitEta_cpy = *Hit_Eta;
      HitProcId_cpy = *Hit_ProcId;
      HitEout_cpy = *Hit_Eout;
      HitEin_cpy = *Hit_Ein;
      HitDetId_cpy = *Hit_DetId;
      HitGunId_cpy = *Hit_GunId;

	
      for (int iHit=0 ; iHit< nHits; iHit++){                               	 
	auto min_id_ptr = std::min_element(HitTime_cpy.begin(), HitTime_cpy.end());
	size_t min_idx = std::distance(HitTime_cpy.begin(), min_id_ptr);

	HitPosX_sort.push_back(HitPosX_cpy[min_idx]);
	HitPosY_sort.push_back(HitPosY_cpy[min_idx]);
	HitPosZ_sort.push_back(HitPosZ_cpy[min_idx]);
	HitTime_sort.push_back(HitTime_cpy[min_idx]);
	HitScatAng_sort.push_back(HitScatAng_cpy[min_idx]);
	HitEta_sort.push_back(HitEta_cpy[min_idx]);
	HitProcId_sort.push_back(HitProcId_cpy[min_idx]);
	HitEout_sort.push_back(HitEout_cpy[min_idx]);
	HitEin_sort.push_back(HitEin_cpy[min_idx]);
	HitDetId_sort.push_back(HitDetId_cpy[min_idx]);
	HitGunId_sort.push_back(HitGunId_cpy[min_idx]);
	     
	HitTime_cpy.erase(HitTime_cpy.begin() + min_idx);
	HitPosX_cpy.erase(HitPosX_cpy.begin() + min_idx);
	HitPosY_cpy.erase(HitPosY_cpy.begin() + min_idx);
	HitPosZ_cpy.erase(HitPosZ_cpy.begin() + min_idx);
	HitScatAng_cpy.erase(HitScatAng_cpy.begin() + min_idx);
	HitEta_cpy.erase(HitEta_cpy.begin() + min_idx);
	HitProcId_cpy.erase(HitProcId_cpy.begin() + min_idx);
	HitEout_cpy.erase(HitEout_cpy.begin() + min_idx);
	HitEin_cpy.erase(HitEin_cpy.begin() + min_idx);
	HitDetId_cpy.erase(HitDetId_cpy.begin() + min_idx);
	HitGunId_cpy.erase(HitGunId_cpy.begin() + min_idx);
      }
      //cout << "\n\nWorking till here!" << endl;
      //Sorting done
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


      ///////////Hit Proc Id Convention used while generating the tree////////////
      //
      //0: Rayleigh
      //1: Photoelectric
      //2: Compton
      //
      ///////////////////////////////////////////////////////////////////////////

      vector<vector<double>> Hits;
      for (int i = 0; i< HitTime_sort.size(); i++){
        vector<double> info;
	
	info.push_back(HitTime_sort[i]);
	info.push_back(HitPosX_sort[i]);
	info.push_back(HitPosY_sort[i]);
	info.push_back(HitPosZ_sort[i]);
	info.push_back(HitScatAng_sort[i]*180.0 / TMath::Pi());	  
	info.push_back(HitEin_sort[i]);
	info.push_back(HitDetId_sort[i]);
	info.push_back(HitGunId_sort[i]);
	info.push_back(HitEta_sort[i]*180.0 / TMath::Pi());

	Hits.push_back(info);
      }
     
      enum hitInfo{time, PosX, PosY, PosZ, scat_theta, Ein, DetId, GunId, Eta};//, Pol_X, Pol_Y, Pol_Z, ScatMomX, ScatMomY, ScatMomZ};      
      h_nHits->Fill(nHits);  //before looking at the histogram check which kind of events you want to look for

      if(nHits>0){
	for (int iHit=0 ; iHit< 1; iHit++){	 
	  double Time = HitTime_sort[iHit];
	  double posX = HitPosX_sort[iHit];
	  double posY = HitPosY_sort[iHit];
	  double posZ = HitPosZ_sort[iHit];

	  double eout = HitEout_sort[iHit];	  
	  double ein = HitEin_sort[iHit];	 

	  int detId = HitDetId_sort[iHit];
	  int gunId = HitGunId_sort[iHit];
	    
	  double scat_theta_SD = HitScatAng_sort[iHit];
	  double scat_theta_SD_deg = scat_theta_SD*180.0 / TMath::Pi();

	  double scat_eta = HitEta_sort[iHit];
	  double scat_eta_deg = scat_eta*180.0 / TMath::Pi();

	  int procId = HitProcId_sort[iHit];
	  
	  string procName;
	  if(procId == 0) procName = "Rayl";
	  if(procId == 1) procName = "Photo";
	  if(procId == 2) procName = "Compt";

	  if(procName == "Compt" && gunId ==1 && posZ>0) FillHistogram(time, posX, posY, posZ, scat_theta_SD_deg, scat_eta_deg, ein, eout);
	} //end of iHit loop	  	 	    	 
      }   

      
      //--------------------Analysing the Scintillation photons at the end---------------------------

      //------------ Making an array of modules-----------------
      enum ModIdx{r=0, l=1};
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
	  
	    E_scatElec[mdl] = Modules[mdl][i][j][E_elec];
	  

	    if(nOp_crys[mdl] != 0) {	    
	      h_DetX_crys[mdl]->Fill(DetX[mdl]);
	      h_nOpPho_crys_inc[mdl]->Fill(nOp_crys[mdl]);
	      h_nOpPho_gen_crys_inc[mdl]->Fill(nOp_gen_crys[mdl]);
	      
	      h_nOp_vs_Edep_crys[mdl]->Fill(EDep_crys[mdl], nOp_crys[mdl]);
	      //Edep_reco_crys[mdl] = (1/21988.2) * nOp_crys[mdl];
	      Edep_reco_crys[mdl] = (1/21926.3) * nOp_crys[mdl];

	      h_Edep_truth_crys[mdl]->Fill(EDep_crys[mdl]);
	      h_Edep_reco_crys[mdl]->Fill(Edep_reco_crys[mdl]);
	      h_Edep_true_Vs_reco_crys[mdl]->Fill(EDep_crys[mdl], Edep_reco_crys[mdl]);
	    }
			  
	    nOp_mdl[mdl] += nOp_crys[mdl];
	    nOp_gen_mdl[mdl] += nOp_gen_crys[mdl];
	    
	    EDep_truth_mdl[mdl] += EDep_crys[mdl];
	    EDep_reco_mdl[mdl] += Edep_reco_crys[mdl];

	    //-------------------getting the relavant hits using true Edep---------------------
	    if(EDep_crys[mdl] > 0.06 && EDep_crys[mdl]< 0.405){	    
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
	    if(Edep_reco_crys[mdl] > 0.06-(3*sigma) && Edep_reco_crys[mdl]< 0.405+(3*sigma)){	    
	      nHits_mdl_reco[mdl] ++;
	      Hit_Edep_reco[mdl].push_back(Edep_reco_crys[mdl]);	        

	      Hit_DetX_reco[mdl].push_back(DetX[mdl]);
	      Hit_DetY_reco[mdl].push_back(DetY[mdl]);

	      //h_Edep_reco_crys[mdl]->Fill(Edep_reco_crys[mdl]);
	      h_DetXvsDetY_reco[mdl]->Fill(DetX[mdl], DetY[mdl]);	  	      
	    }
	  }
	}
	//---------------------------------------------------------------------------------------------------
	//truth hists
	h_Edep_truth_mdl[mdl] ->Fill(EDep_truth_mdl[mdl]);
	h_nHits_truth_mdl[mdl]->Fill(nHits_mdl_truth[mdl]);

	//reco hists
	h_nOpPho_mdl_inc[mdl]->Fill(nOp_mdl[mdl]);
	h_nOpPho_gen_mdl_inc[mdl]->Fill(nOp_gen_mdl[mdl]);
	h_Edep_reco_mdl[mdl] ->Fill(EDep_reco_mdl[mdl]);
	h_nHits_reco_mdl[mdl]->Fill(nHits_mdl_reco[mdl]);

	//for calibration (using both truth and reco info at the same time)
	h_nOp_vs_Edep_mdl[mdl]->Fill(EDep_truth_mdl[mdl], nOp_mdl[mdl]);	
      }       
      
      //-------------for getting the sigma in the sum case
      if(nHits_mdl_truth[r]==2 && Hit_Edep_truth[r][0]+Hit_Edep_truth[r][1] <0.5111 && Hit_Edep_truth[r][0]+Hit_Edep_truth[r][1] >0.5109){
      	  h_nOpPho_mdl_sum[r]->Fill(Hit_nOp_indv[r][0]+Hit_nOp_indv[r][1]);
	  h_nOpPho_gen_mdl_sum[r]->Fill(Hit_nOp_gen_indv[r][0]+Hit_nOp_gen_indv[r][1]);
      }

      if(nHits_mdl_truth[l]==2 && Hit_Edep_truth[l][0]+Hit_Edep_truth[l][1] <0.5111 && Hit_Edep_truth[l][0]+Hit_Edep_truth[l][1] >0.5109){
	h_nOpPho_mdl_sum[l]->Fill(Hit_nOp_indv[l][0]+Hit_nOp_indv[l][1]);
	h_nOpPho_gen_mdl_sum[l]->Fill(Hit_nOp_gen_indv[l][0]+Hit_nOp_gen_indv[l][1]);
      }

//===============================Calculating theta using truth Edep====================================
      if(nHits_mdl_truth[r]==2 && nHits_mdl_truth[l]==2){      
  	if(Hit_Edep_truth[r][0]+Hit_Edep_truth[r][1] <0.511 && Hit_Edep_truth[r][0]+Hit_Edep_truth[r][1] >=0.510){
	  if(Hit_Edep_truth[l][0]+Hit_Edep_truth[l][1] <0.512 && Hit_Edep_truth[l][0]+Hit_Edep_truth[l][1] >=0.510){
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
	    
	    //if(scat_Edep < 0.0476 && scat_Edep > 0.0474) cout << "look for Evt: " << jentry +1 << endl;
	  
	    theta[Mdl] = GetScatAngle(scat_Edep[Mdl]);	  
	    
	    //if (scat_Edep[Mdl] < 0.25){	    
	    h_theta_truth[Mdl]->Fill(theta[Mdl]*180/ TMath::Pi());	     
	    h_elecE_vs_scatEdep[Mdl]->Fill(scat_Edep[Mdl], Hit_Scat_elecE[Mdl][scat_idx[Mdl]]);
	    
	    //}

  	  //if (abs(scat_Edep-Hit_Scat_elecE[idx]) < 0.001) h_theta_truth[Mdl]->Fill(theta_truth*180/ TMath::Pi());

  	  // for (int i=0; i<Hit_Edep_truth.size(); i++){
  	  //   if (abs(Hit_Edep_truth[i] - Hit_Scat_elecE[i]) < 0.0001){
  	  //     theta_truth = GetScatAngle(Hit_Edep_truth[i]);
  	  //     h_theta_truth[Mdl]->Fill(theta_truth*180/ TMath::Pi());
  	  //   }
  	  // }

  	  
  	    h_Edep_scat_truth[Mdl]->Fill(scat_Edep[Mdl]);
  	    h_Edep_abs_truth[Mdl]->Fill(abs_Edep[Mdl]);
  	    h_Edep_scatVsabs_truth[Mdl]->Fill(scat_Edep[Mdl], abs_Edep[Mdl]);  	  	  
	  
  	  //------------------Calculating Phi---------------------------------------
  	  float abs_DetX = Hit_DetX_truth[Mdl][abs_idx[Mdl]];
  	  float abs_DetY = Hit_DetY_truth[Mdl][abs_idx[Mdl]];
	  
  	  float scat_DetX = Hit_DetX_truth[Mdl][scat_idx[Mdl]];
  	  float scat_DetY = Hit_DetY_truth[Mdl][scat_idx[Mdl]];

  	  //phi_truth[Mdl] = atan((abs_DetY_truth-scat_DetY_truth)/(abs_DetX_truth-scat_DetX_truth));	
	  phi[Mdl] = atan2(abs_DetY-scat_DetY, abs_DetX-scat_DetX);
	  
  	  //h_phi_truth[Mdl]->Fill(phi_truth[Mdl]*180/ TMath::Pi());
	  dPix[Mdl]= sqrt(((abs_DetX-scat_DetX)*(abs_DetX-scat_DetX))+((abs_DetY-scat_DetY)*(abs_DetY-scat_DetY)));

	  //cout << "pix distance: " << Pix_distance[Mdl] << endl;
	  h_dPix_truth[Mdl]->Fill(dPix[Mdl]);
	  }
	  }	  
	  
	  if (theta[r]< Rad(90) && theta[r]> Rad(72) && theta[l]< Rad(90) && theta[l]> Rad(72)){
	    h_phi_truth[r]->Fill(Deg(phi[r]));
	    h_phi_truth[l]->Fill(Deg(phi[l]));
	    h_Phi1_vs_Phi2_truth->Fill(Deg(phi[r]), Deg(phi[l]));
	    
	    float dPhi = phi[r] - phi[l];
	    if (abs(dPhi) <= Pi) h_dPhi_truth->Fill(dPhi*180/ Pi);
	    else if(abs(dPhi) > Pi && dPhi>0) h_dPhi_truth->Fill((dPhi-(2*Pi))*180/ Pi);
	    else if(abs(dPhi) > Pi && dPhi<0) h_dPhi_truth->Fill((dPhi+(2*Pi))*180/ Pi);	    
	    }
	  //}
	  }
	}
      }

      //===============================Calculating theta using Reco Edep====================================
      if(nHits_mdl_reco[r]==2 && nHits_mdl_reco[l]==2){      
  	if(Hit_Edep_reco[r][0]+Hit_Edep_reco[r][1] < (0.511+(3*sigma)) && Hit_Edep_reco[r][0]+Hit_Edep_reco[r][1] >=(0.511-3*sigma)){	  
	  if(Hit_Edep_reco[l][0]+Hit_Edep_reco[l][1] < (0.511+(3*sigma)) && Hit_Edep_reco[l][0]+Hit_Edep_reco[l][1] >=(0.511-3*sigma)){
	    float scat_Edep[2] ={0}, abs_Edep[2] ={0}, theta[2] = {9999};
	    int scat_idx[2] ={99}, abs_idx[2] ={99};
	    double phi[2] ={999}, dPix[2] = {99999};
	    
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
	  
	  if (theta[r]< Rad(90) && theta[r]> Rad(72) && theta[l]< Rad(90) && theta[l]> Rad(72)){
	    h_phi_reco[r]->Fill(Deg(phi[r]));
	    h_phi_reco[l]->Fill(Deg(phi[l]));
	    h_Phi1_vs_Phi2_reco->Fill(Deg(phi[r]), Deg(phi[l]));
	    
	    float dPhi = phi[r] - phi[l];
	    if (abs(dPhi) <= Pi) h_dPhi_reco->Fill(dPhi*180/ Pi);
	    else if(abs(dPhi) > Pi && dPhi>0) h_dPhi_reco->Fill((dPhi-(2*Pi))*180/ Pi);
	    else if(abs(dPhi) > Pi && dPhi<0) h_dPhi_reco->Fill((dPhi+(2*Pi))*180/ Pi);	    
	    }
	  //}
	  }
	}
      }

  


      // //Calculating theta using reco Edep
      // if(nHits_mdl[r]==2) {
      // 	if(Hit_Edep_reco[r][0]+Hit_Edep_reco[r][1] <0.512 && Hit_Edep_reco[r][0]+Hit_Edep_reco[r][1] >=0.510){	
      // 	  if (Hit_Edep_reco[r][0] < Hit_Edep_reco[r][1]) theta_reco[r] = GetScatAngle(Hit_Edep_reco[r][0]);
      // 	  else if (Hit_Edep_reco[r][1] < Hit_Edep_reco[r][0]) theta_reco[r] = GetScatAngle(Hit_Edep_reco[r][1]);  
      // 	  h_theta_reco[r]->Fill(theta_reco[r]*180/ TMath::Pi());	  
      // 	}
      // }
      
           
      } //jentry loop end

  cout << "No. of events with exactly 2 Hits: " << nEvts2Hits << endl;

  cout << "\nTheta: " << Pi << endl;
  cout << "\nScat ElecE: " << GetScatElecE(90) << endl;  
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
