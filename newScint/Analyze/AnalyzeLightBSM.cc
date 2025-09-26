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

      h_CreateZ->Fill(OpPho_Z);
      h_OpPho_lmbda->Fill(OpPho_Lambda);
      h_lmbdVsZ->Fill(OpPho_Z, OpPho_Lambda);
       

    } //jentry loop end
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
