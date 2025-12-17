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

#pragma link C++ class std::vector< std::vector >+; 
#pragma link C++ class std::vector< TLorentzVector >+;
#ifdef __CINT__
#pragma link C++ class std::vector<TLorentzVector>+;
#endif
using namespace TMVA;
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
   
  for (Long64_t jentry=0; jentry<nentries;jentry++)
    {
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
	}	  	 	    	 
      }   //end of iHit loop

      
      //--------------------Analysing the Scintillation photons at the end---------------------------
      int nOp = 0, nOp_gen =0;
      float EDep = 0;
      for (int i =0; i<8; i++){
      	for (int j =0; j<8; j++){
      	  nOp += right_module[i][j][nOpPho];
      	  nOp_gen += right_module[i][j][nGenOp];
	  EDep += right_module[i][j][Edep];
	  	  
      	  if(right_module[i][j][nOpPho] != 0) {
	    h_DetX_crys->Fill(right_module[i][j][DetPosX]);
	    h_Edep_crys->Fill(right_module[i][j][Edep]);
	    h_nOp_vs_Edep_crys->Fill(right_module[i][j][nOpPho], right_module[i][j][Edep]);
	  }
      	}
      }
      
      // if (nOp > 0) h_nOpPho->Fill(nOp);
      // if (nOp_gen >0) h_nOpPho_gen->Fill(nOp_gen);
      h_nOpPho_mdl->Fill(nOp);
      h_nOpPho_gen_mdl->Fill(nOp_gen);

      h_Edep_mdl ->Fill(EDep);

      h_nOp_vs_Edep_mdl->Fill(nOp, EDep);
      
      
    } //jentry loop end

  // TF1 *f1 = new TF1("","[0]+[1]*x",0,1);
  // f1->SetParameters(0.,1.);
  // f1->SetLineColor(kRed);
  // h2->Fit(f1);
  // h2->Draw();
  // f1->Draw("same");
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
