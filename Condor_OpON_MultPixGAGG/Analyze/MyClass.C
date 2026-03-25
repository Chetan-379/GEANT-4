#define MyClass_cxx
#include "MyClass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void MyClass::Loop()
{
//   In a ROOT session, you can do:
//      root> .L MyClass.C
//      root> MyClass t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   Long64_t nbytes = 0, nb = 0;
   int decade = 0;
   bool Debug=false;
   int Full_edep_evts = 0;
   int Comp2_evts = 0;
   int nDiffAngle = 0, nHits = -1, nRayl = 0;
   
   int n1s = 0, n1 =0;
   int n1b =0;
   
   int nEvts2Hits = 0;
   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

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
      float EDep_mdl[2] = {0.}, EDep_crys[2] = {0.}, DetX[2] = {0.}, DetY[2] = {0.};
      float Edep_diff[2] = {0.};
      float Edep_reco_crys[2] = {0.};
      float E_scatElec[2] = {0.};

      int nHits_mdl[2] = {0}, nHits_mdl_truth[2] = {0};

      vector<vector<float>> Hit_Edep_reco, Hit_Edep_truth, Hit_Scat_elecE, Hit_DetX, Hit_DetY, Hit_DetX_truth, Hit_DetY_truth;
      Hit_Edep_reco.resize(2);
      Hit_Edep_truth.resize(2);
      Hit_Scat_elecE.resize(2);
      Hit_DetX.resize(2);
      Hit_DetY.resize(2);
      Hit_DetX_truth.resize(2);
      Hit_DetY_truth.resize(2);          

      float scat_Edep[2] = {0}, abs_Edep[2] = {0}, theta_reco[2] = {0}, theta_truth[2] = {99999};

      

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
	      h_nOpPho_crys[mdl]->Fill(nOp_crys[mdl]);	    
	      h_nOp_vs_Edep_crys[mdl]->Fill(EDep_crys[mdl], nOp_crys[mdl]);
	      //Edep_reco_crys[mdl] = (1/21988.2) * nOp_crys[mdl];
	      Edep_reco_crys[mdl] = (1/21929.7) * nOp_crys[mdl];
	      h_Edep_true_Vs_reco_crys[mdl]->Fill(EDep_crys[mdl], Edep_reco_crys[mdl]);
	      h_DetXvsDetY[mdl]->Fill(DetX[mdl], DetY[mdl]);
	    
	      nHits_mdl[mdl]++;

	      Hit_Edep_reco[mdl].push_back(Edep_reco_crys[mdl]);
	      Hit_DetX[mdl].push_back(DetX[mdl]);
	      Hit_DetY[mdl].push_back(DetY[mdl]);	   
	    }
			  
	    nOp_mdl[mdl] += nOp_crys[mdl];
	    nOp_gen_mdl[mdl] += nOp_gen_crys[mdl];
	    EDep_mdl[mdl] += EDep_crys[mdl];


	    //if(EDep_crys[mdl] > 0.0001){
	    if(EDep_crys[mdl] > 0.06 && EDep_crys[mdl]< 0.405){	    
	      nHits_mdl_truth[mdl] ++;
	      Hit_Edep_truth[mdl].push_back(EDep_crys[mdl]);
	      Hit_Scat_elecE[mdl].push_back(E_scatElec[mdl]);
	      h_elecE_vs_Edep_crys[mdl]->Fill(EDep_crys[mdl], E_scatElec[mdl]);
	      h_Edep_crys[mdl]->Fill(EDep_crys[mdl]);

	      //if(abs(EDep_crys - E_scatElec) < 0.001){
  	      //auto theta_truth = GetScatAngle(EDep_crys); 
  	      //h_theta_truth[mdl]->Fill(theta_truth*180/ TMath::Pi());
	      //}
	      Hit_DetX_truth[mdl].push_back(DetX[mdl]);
	      Hit_DetY_truth[mdl].push_back(DetY[mdl]);	      
	    }
	  }
	}
      }


      //cout << nHits_mdl_truth << endl;

      //Calculating theta using truth Edep
      if(nHits_mdl_truth[r]==2 && nHits_mdl_truth[l]==2){
	//if(nHits_mdl_truth[l]==2){
      
  	if(Hit_Edep_truth[r][0]+Hit_Edep_truth[r][1] <0.512 && Hit_Edep_truth[r][0]+Hit_Edep_truth[r][1] >=0.510){	  
	  if(Hit_Edep_truth[l][0]+Hit_Edep_truth[l][1] <0.512 && Hit_Edep_truth[l][0]+Hit_Edep_truth[l][1] >=0.510){

  	  // h_Edep_crys->Fill(Hit_Edep_truth[0]);
  	  // h_Edep_crys->Fill(Hit_Edep_truth[1]);

  	  float scat_Edep[2] ={0}, abs_Edep[2] ={0};
  	  int scat_idx[2] ={99}, abs_idx[2] ={99};
	  double phi_truth[2] ={999}, dPix[2] = {99999};

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
	  
	    theta_truth[Mdl] = GetScatAngle(scat_Edep[Mdl]);	  
	    
	    //if (scat_Edep[Mdl] < 0.25){	    
	    h_theta_truth[Mdl]->Fill(theta_truth[Mdl]*180/ TMath::Pi());	     
	    h_elecE_vs_scatEdep[Mdl]->Fill(scat_Edep[Mdl], Hit_Scat_elecE[Mdl][scat_idx[Mdl]]);
	    
	    //}

  	  //if (abs(scat_Edep-Hit_Scat_elecE[idx]) < 0.001) h_theta_truth[Mdl]->Fill(theta_truth*180/ TMath::Pi());

  	  // for (int i=0; i<Hit_Edep_truth.size(); i++){
  	  //   if (abs(Hit_Edep_truth[i] - Hit_Scat_elecE[i]) < 0.0001){
  	  //     theta_truth = GetScatAngle(Hit_Edep_truth[i]);
  	  //     h_theta_truth[Mdl]->Fill(theta_truth*180/ TMath::Pi());
  	  //   }
  	  // }

  	  if (Hit_Edep_truth[Mdl][0] < Hit_Edep_truth[Mdl][1]){
  	    h_Edep_scat[Mdl]->Fill(Hit_Edep_truth[Mdl][0]);
  	    h_Edep_abs[Mdl]->Fill(Hit_Edep_truth[Mdl][1]);
  	    h_Edep_scatVsabs[Mdl]->Fill(Hit_Edep_truth[Mdl][0], Hit_Edep_truth[Mdl][1]);
  	  }
	  
  	  else if (Hit_Edep_truth[Mdl][1] < Hit_Edep_truth[Mdl][0]){
  	    h_Edep_scat[Mdl]->Fill(Hit_Edep_truth[Mdl][1]);
  	    h_Edep_abs[Mdl]->Fill(Hit_Edep_truth[Mdl][0]);
  	    h_Edep_scatVsabs[Mdl]->Fill(Hit_Edep_truth[Mdl][1], Hit_Edep_truth[Mdl][0]);
  	  }
	  	  
	  
  	  //------------------Calculating Phi---------------------------------------
  	  float abs_DetX_truth = Hit_DetX_truth[Mdl][abs_idx[Mdl]];
  	  float abs_DetY_truth = Hit_DetY_truth[Mdl][abs_idx[Mdl]];
	  
  	  float scat_DetX_truth = Hit_DetX_truth[Mdl][scat_idx[Mdl]];
  	  float scat_DetY_truth = Hit_DetY_truth[Mdl][scat_idx[Mdl]];

  	  //phi_truth[Mdl] = atan((abs_DetY_truth-scat_DetY_truth)/(abs_DetX_truth-scat_DetX_truth));	
	  phi_truth[Mdl] = atan2(abs_DetY_truth-scat_DetY_truth, abs_DetX_truth-scat_DetX_truth);
	  
  	  //h_phi_truth[Mdl]->Fill(phi_truth[Mdl]*180/ TMath::Pi());
	  dPix[Mdl]= sqrt(((abs_DetX_truth-scat_DetX_truth)*(abs_DetX_truth-scat_DetX_truth))+((abs_DetY_truth-scat_DetY_truth)*(abs_DetY_truth-scat_DetY_truth)));

	  //cout << "pix distance: " << Pix_distance[Mdl] << endl;
	  h_dPix[Mdl]->Fill(dPix[Mdl]);
	  }
	  }	  
	  
	  if (theta_truth[r]< Rad(90) && theta_truth[r]> Rad(72) && theta_truth[l]< Rad(90) && theta_truth[l]> Rad(72)){
	    h_phi_truth[r]->Fill(Deg(phi_truth[r]));
	    h_phi_truth[l]->Fill(Deg(phi_truth[l]));
	    h_Phi1_vs_Phi2_truth->Fill(Deg(phi_truth[r]), Deg(phi_truth[l]));
	    
	    float dPhi_truth = phi_truth[r] - phi_truth[l];
	    if (abs(dPhi_truth) <= Pi) h_dPhi_truth->Fill(dPhi_truth*180/ Pi);
	    else if(abs(dPhi_truth) > Pi && dPhi_truth>0) h_dPhi_truth->Fill((dPhi_truth-(2*Pi))*180/ Pi);
	    else if(abs(dPhi_truth) > Pi && dPhi_truth<0) h_dPhi_truth->Fill((dPhi_truth+(2*Pi))*180/ Pi);	    
	    }
	  //}
	  }
	}
      }
      
  


      //Calculating theta using reco Edep
      if(nHits_mdl[r]==2) {
  	if(Hit_Edep_reco[r][0]+Hit_Edep_reco[r][1] <0.512 && Hit_Edep_reco[r][0]+Hit_Edep_reco[r][1] >=0.510){	
  	  if (Hit_Edep_reco[r][0] < Hit_Edep_reco[r][1]) theta_reco[r] = GetScatAngle(Hit_Edep_reco[r][0]);
  	  else if (Hit_Edep_reco[r][1] < Hit_Edep_reco[r][0]) theta_reco[r] = GetScatAngle(Hit_Edep_reco[r][1]);  
  	  h_theta_reco[r]->Fill(theta_reco[r]*180/ TMath::Pi());	  
  	}
      }
      
      h_nOpPho_mdl[r]->Fill(nOp_mdl[r]);
      h_nOpPho_gen_mdl[r]->Fill(nOp_gen_mdl[r]);
      h_Edep_mdl[r] ->Fill(EDep_mdl[r]);
      h_nOp_vs_Edep_mdl[r]->Fill(EDep_mdl[r], nOp_mdl[r]);             
      //h_nHits_mdl->Fill(nHits_mdl);
      h_nHits_mdl[r]->Fill(nHits_mdl_truth[r]);     
      

      // if (nOp > 0) h_nOpPho->Fill(nOp);
      // if (nOp_gen >0) h_nOpPho_gen->Fill(nOp_gen);

      //if(EDep_mdl >= 0.444 && EDep_mdl <= 0.477) cout << jentry << endl;

      
   } //jentry loop end 
}
