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
  //int nRelvEvts =0;
   
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

      //if (nHits>2){	
      ROOT::Math::XYZVector v_comp, vin, vout, v1, v2;

      //sorting the hits in the ascending order of time
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      vector<double> HitPosX_sort, HitPosY_sort, HitPosZ_sort, HitTime_sort, HitScatAng_sort, HitEta_sort, HitPol0_sort, HitPol1_sort, HitPol2_sort, HitEin_sort, HitEout_sort;
      vector<double> HitPosX_cpy, HitPosY_cpy, HitPosZ_cpy, HitTime_cpy, HitScatAng_cpy, HitEta_cpy, HitPol0_cpy, HitPol1_cpy, HitPol2_cpy, HitEin_cpy, HitEout_cpy;

      vector <int> HitProcId_sort, HitDetId_sort, HitGunId_sort;
      vector<int> HitProcId_cpy, HitDetId_cpy, HitGunId_cpy;
	   
      HitPosX_sort.clear();
      HitPosY_sort.clear();
      HitPosZ_sort.clear();
      HitTime_sort.clear();
      HitScatAng_sort.clear();
      HitEta_sort.clear();
      HitProcId_sort.clear();
      HitPol0_sort.clear();
      HitPol1_sort.clear();
      HitPol2_sort.clear();	
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
      HitPol0_cpy = *Hit_Pol0;
      HitPol1_cpy = *Hit_Pol1;
      HitPol2_cpy = *Hit_Pol2;
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
	HitPol0_sort.push_back(HitPol0_cpy[min_idx]);
	HitPol1_sort.push_back(HitPol1_cpy[min_idx]);
	HitPol2_sort.push_back(HitPol2_cpy[min_idx]);
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
	HitPol0_cpy.erase(HitPol0_cpy.begin() + min_idx);
	HitPol1_cpy.erase(HitPol1_cpy.begin() + min_idx);
	HitPol2_cpy.erase(HitPol2_cpy.begin() + min_idx);
	HitEout_cpy.erase(HitEout_cpy.begin() + min_idx);
	HitEin_cpy.erase(HitEin_cpy.begin() + min_idx);
	HitDetId_cpy.erase(HitDetId_cpy.begin() + min_idx);
	HitGunId_cpy.erase(HitGunId_cpy.begin() + min_idx);
      }

      //Sorting done
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


      ///////////Hit Proc Id Convention used while generating the tree////////////
      //
      //0: Rayleigh
      //1: Photoelectric
      //2: Compton
      //
      ///////////////////////////////////////////////////////////////////////////
      
      bool incrs_nHit = false;
      //bool RejEvt = false;
	
      if (nHits>0 && HitProcId_sort[0] == 2) incrs_nHit = true;

      for (int i = 0; i< HitTime_sort.size(); i++){
	if(incrs_nHit && HitProcId_sort[i] != 0) nHitComp++;   //excluding rayleigh hit
	//if(incrs_nHit) nHitComp++;
	//if (HitProcId_sort[i] == 0) RejEvt = true;
      }

      //if (RejEvt) continue;   //rejecting the event the event with any rayl just for checking 

      h_nHits->Fill(nHitComp);

      
      //if(nHitComp >0){
      if(nHits >0){
      //if(nHitComp >1){
	vector<vector<double>> AnniPho, ScatPho, Relv4Hits, Ana_Relv4Hits, AS_ConnHits;
	int nHits_RaylEx = 0;
	for (int iHit=0 ; iHit< nHits; iHit++){     //cannot calculate theta at the last hit
	//for (int iHit=0 ; iHit < 1; iHit++){     //considering only the first hit
	  
	  v_comp.SetXYZ(0,0,0);
	  vin.SetXYZ(0,0,0);
	  vout.SetXYZ(0,0,0);
	  v1.SetXYZ(0,0,0);
	  v2.SetXYZ(0,0,0);

	  v_comp.SetXYZ(HitPosX_sort[iHit], HitPosY_sort[iHit], HitPosZ_sort[iHit]);

	  if(iHit == 0) v1.SetXYZ(0,0,0);            //for first hit initial point will be gun postion
	  else v1.SetXYZ(HitPosX_sort[iHit-1], HitPosY_sort[iHit-1], HitPosZ_sort[iHit-1]);

	  v2.SetXYZ(HitPosX_sort[iHit+1], HitPosY_sort[iHit+1], HitPosZ_sort[iHit+1]);

	  vin = v_comp - v1;
	  vout = v2 - v_comp;

	  double scat_theta = 10000;

	  scat_theta = acos(vin.Dot(vout)/(vin.R()*vout.R()));
	 
	  double Time = HitTime_sort[iHit];
	  double PosX = HitPosX_sort[iHit];
	  double PosY = HitPosY_sort[iHit];
	  double PosZ = HitPosZ_sort[iHit];

	  double Eout = HitEout_sort[iHit];	  
	  double Ein = HitEin_sort[iHit];

	  int DetId = HitDetId_sort[iHit];
	  int GunId = HitGunId_sort[iHit];

	    
	  double scat_theta_SD = HitScatAng_sort[iHit];
	  double Ang_diff = scat_theta - scat_theta_SD;
	  double scat_theta_deg = scat_theta*180.0 / TMath::Pi();
	  double scat_theta_SD_deg = scat_theta_SD*180.0 / TMath::Pi();	
	  //h_ScatAng_SDvsAna->Fill(scat_theta, scat_theta_SD);

	  string procName;
	  if(HitProcId_sort[iHit] == 0) procName = "Rayl";
	  if(HitProcId_sort[iHit] == 1) procName = "Photo";
	  if(HitProcId_sort[iHit] == 2) procName = "Compt";

	  float AnniPho_scat_theta_SD = -1000;

	  if (procName != "Rayl") nHits_RaylEx++ ;
	  if (procName == "Compt" && Ein == 0.511) AnniPho_scat_theta_SD = scat_theta_SD;  //be careful here, there can be rayl before compt

	  float AnniPho_scat_theta_SD_deg = AnniPho_scat_theta_SD*180.0 / TMath::Pi();

	  //vector<vector<double>> AnniPho;
	  vector<double> info;

	  info.push_back(Time);
	  info.push_back(PosX);
	  info.push_back(PosY);
	  info.push_back(PosZ);
	  info.push_back(scat_theta_SD_deg);
	  info.push_back(Ein);
	  info.push_back(DetId);
	  info.push_back(GunId);

	  if (AnniPho_scat_theta_SD > 0) {	    
	    AnniPho.push_back(info);
	  }

	  if (procName != "Rayl" && Ein != 0.511) ScatPho.push_back(info); 
	  	  	  	  
	  double EtaPolScat = HitEta_sort[iHit];
	  //double EtaPolScat = -100000;
	  //if (nHitComp >=2) EtaPolScat = HitEta_sort[1];
	  double EtaPolScat_deg = EtaPolScat*180.0 / TMath::Pi();
	  int procId = HitProcId_sort[iHit];

	  float Pol0 = HitPol0_sort[iHit];
	  float Pol1 = HitPol1_sort[iHit];
	  float Pol2 = HitPol2_sort[iHit];

	  
	 	    
	  h_diff_Ana_SD->Fill(Ang_diff);
	  h_posX_posY ->Fill(PosX, PosY);
	  h_polX_polY ->Fill(Pol0, Pol1);

	  //vector<double> Pol_Ana;
	  ROOT::Math::XYZVector Pol_Ana;
	  Pol_Ana = vin.Cross(vout)/(vin.R()*vout.R());

	  h_PolX_Ana_SD->Fill(Pol_Ana.X(), Pol0);
	  h_PolY_Ana_SD->Fill(Pol_Ana.Y(), Pol1);
	  h_PolZ_Ana_SD->Fill(Pol_Ana.Z(), Pol2);

	  //ROOT::Math::XYZVector Pol_SD(Pol0, Pol1, Pol2);

	  //if (EtaPolScat_deg == 90) 

	  h_Pol_Ana->Fill(Pol_Ana.Mag2());
	  

	  
	  //***************************************
	  //convention followed for nHits ==2:
	  //S: Scintillator, B: BGO
	  //left for location of 1st hit;
	  //right for 2nd;
	  //1 and 2 in S or B for identifying same or different crystal of a material;
	  //no identifier when hits are in different materials;
	  //e.g. nHit2_S1S1 means 1st hit in scintillator 2nd hit in the same scinitillator
	  //***************************************
	  	  
	  enum HitsLoc{all_inc =0, nHit0 =1, nHit1_S =2, nHit1_B =3, nHit2_S_AnyDiffDet =4, nHit2_S1S1 =5, nHit2_S1S2 =6, nHit2_SB =7, nHit2_B1B1 =8, nHit2_B1B2 =9, nHit2_BS =10, remain =11, test = 12};

	  //Definitng the booleans for each category
	  bool inc_all =false, noHit =false;
	  bool hit1_S =false, hit1_B =false;
	  bool hit2_S_AnyDiffDet = false, hit2_S1S1 =false, hit2_S1S2 =false, hit2_SB =false, hit2_B1B1 =false, hit2_B1B2 =false, hit2_BS =false, rest =false;
	  //bool teSt = false;
	  bool teSt = true;
	  	  
	  int DetId1 = HitDetId_sort[0];
	  int DetId2 = -1;
	  if(nHitComp >=2) DetId2 = HitDetId_sort[1];

	  if (nHitComp >0) inc_all = true;
	  
	  //if (nHitComp > 0 && MatName(DetId1) == "BGO") {hit_BGO_inc = true; DetId1 = HitDetId_sort[0];}
	  // if (nHitComp >=2 && DetId1 != DetId2 && MatName(DetId1) == "Scint") {
	  //   teSt = true;
	  // }
	  
	  if (nHitComp ==0) noHit = true;
	  
	  else if (nHitComp ==1){
	    n1++;
	    if (MatName(DetId1) =="Scint") {n1s++; hit1_S = true;} 
	    else if (MatName(DetId1) =="BGO") {n1b++; hit1_B = true;}
	    else rest = true;
	  }

	  // else if (nHitComp ==2){
	  else if (nHitComp >=2){	  	    
	     if (DetId1 == DetId2){	   		    
	      if (MatName(DetId1) =="Scint") hit2_S1S1 = true; 
	      else if (MatName(DetId1) =="BGO") hit2_B1B1 = true;
	      else rest = true;	      
	    }	    
	    
	     else if (DetId1 != DetId2)
	       {if(MatName(DetId1) =="Scint") hit2_S_AnyDiffDet = true;
		 if(MatName(DetId1) =="Scint" && MatName(DetId2) =="Scint") {hit2_S1S2 = true;} //cout << "event with s1s2: " << jentry +1 << endl;}
		 else if(MatName(DetId1) =="BGO" && MatName(DetId2) =="BGO") hit2_B1B2 = true;
		 else if (MatName(DetId1) =="Scint" && MatName(DetId2) =="BGO") {hit2_SB = true;}// cout << "event with SB: " << jentry +1 << endl;}	       
		 else if (MatName(DetId1) =="BGO" && MatName(DetId2) == "Scint") hit2_BS = true;
		 else rest = true;
	       }
	  }
	  
	  else rest = true;

	  //if (hit2_S1S1) cout << "two hits in same scint" << endl;

	  //else if (DetId1 != DetId2 && MatName(DetId1) =="Scint" && MatName(DetId2) =="Scint") hit2_S1S2 == true;
	  //else if (MatName(DetId1) =="Scint" && MatName(DetId2) =="BGO") hit2_SB == true;
	  //else if (DetId1 == DetId2 && MatName(DetId1) =="BGO") hit2_B1B1 == true;
	  //else if (DetId1 != DetId2 && MatName(DetId1) =="BGO" && MatName(DetId2) =="BGO") hit2_B1B2 == true;
	  //else if (DetId1 != DetId2 && MatName(DetId1) =="BGO" && MatName(DetId2) == "Scint") hit2_BS == true;
	  

	  
	  
	  if (procId !=0 && procId != 1){       //excluding rayleigh && phot (not required exactly but kept of safety purpose)
	    h_ScatAng_SDvsAna->Fill(scat_theta, scat_theta_SD);
	    
	    if (inc_all) FillHistogram(all_inc, Time, PosX, PosY, PosZ, scat_theta_SD_deg, EtaPolScat_deg, Pol0, Pol1, Pol2, Ein, Eout);
	    if (teSt && scat_theta_SD > 0) {FillHistogram(test, Time, PosX, PosY, PosZ, scat_theta_SD_deg, EtaPolScat_deg, Pol0, Pol1, Pol2, Ein, Eout);	     
	    //if (teSt) FillHistogram(test, Time, PosX, PosY, PosZ, scat_theta_SD_deg, EtaPolScat_deg, Pol0, Pol1, Pol2, Ein, Eout);
	    }
	    if (noHit) FillHistogram(nHit0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	    
	    if (hit1_S) FillHistogram(nHit1_S, Time, PosX, PosY, PosZ, scat_theta_SD_deg, EtaPolScat_deg, Pol0, Pol1, Pol2, Ein, Eout);
	    if (hit1_B) FillHistogram(nHit1_B, Time, PosX, PosY, PosZ, scat_theta_SD_deg, EtaPolScat_deg, Pol0, Pol1, Pol2, Ein, Eout);
	    
	    if (hit2_S_AnyDiffDet) FillHistogram(nHit2_S_AnyDiffDet, Time, PosX, PosY, PosZ, scat_theta_SD_deg, EtaPolScat_deg, Pol0, Pol1, Pol2, Ein, Eout);	      
	    if (hit2_S1S1) FillHistogram(nHit2_S1S1, Time, PosX, PosY, PosZ, scat_theta_SD_deg, EtaPolScat_deg, Pol0, Pol1, Pol2, Ein, Eout);	      
	    if (hit2_S1S2) FillHistogram(nHit2_S1S2, Time, PosX, PosY, PosZ, scat_theta_SD_deg, EtaPolScat_deg, Pol0, Pol1, Pol2, Ein, Eout);
	    if (hit2_SB) FillHistogram(nHit2_SB, Time, PosX, PosY, PosZ, scat_theta_SD_deg, EtaPolScat_deg, Pol0, Pol1, Pol2, Ein, Eout);
	    
	    if (hit2_B1B1) FillHistogram(nHit2_B1B1, Time, PosX, PosY, PosZ, scat_theta_SD_deg, EtaPolScat_deg, Pol0, Pol1, Pol2, Ein, Eout);
	    if (hit2_B1B2) FillHistogram(nHit2_B1B2, Time, PosX, PosY, PosZ, scat_theta_SD_deg, EtaPolScat_deg, Pol0, Pol1, Pol2, Ein, Eout);
	    if (hit2_BS) FillHistogram(nHit2_BS, Time, PosX, PosY, PosZ, scat_theta_SD_deg, EtaPolScat_deg, Pol0, Pol1, Pol2, Ein, Eout);
	    
	    if (rest) FillHistogram(remain, Time, PosX, PosY, PosZ, scat_theta_SD_deg, EtaPolScat_deg, Pol0, Pol1, Pol2, Ein, Eout);	    	    	    	    	  
	  } //procId condition end	

	}   //end of iHit loop

	enum hitInfo{time, PosX, PosY, PosZ, scat_theta, Ein, DetId, GunId};
	enum AnniScat{A1, A2, S1, S2};

	//convention followed:
	//A1: first compton hit from Annihilation
	//A2: second compton hit from Annihilation
	//S1: first hit from Scattered photon
	//S3: 2nd hit from Scattered photon
	  
	if(AnniPho.size()==2 && nHits_RaylEx ==4) {
	  h_nAnniPho->Fill(AnniPho.size());

	  Relv4Hits.push_back(AnniPho[0]);
	  Relv4Hits.push_back(AnniPho[1]);

	  Ana_Relv4Hits.push_back(AnniPho[0]);
	  Ana_Relv4Hits.push_back(AnniPho[1]);

	  AS_ConnHits.push_back(AnniPho[0]);
	  AS_ConnHits.push_back(AnniPho[1]);

	  //Relv4Hits.push_back(ScatPho[0]);
	  //Relv4Hits.push_back(ScatPho[1]);

	  // for (int iscat = 0; iscat<ScatPho.size(); iscat++){
	  //   if(AnniPho[0][GunId] == ScatPho[iscat][GunId]) Relv4Hits.insert(Relv4Hits.begin() + 2, ScatPho[iscat]);
	  //   if(AnniPho[1][GunId] == ScatPho[iscat][GunId]) Relv4Hits.insert(Relv4Hits.begin() + 3, ScatPho[iscat]);
	    
	  // }

	  double GunId_A1 = AnniPho[0][GunId];
	  double GunId_A2 = AnniPho[1][GunId];
	  
	  double GunId_S1 = ScatPho[0][GunId];
	  double GunId_S2 = ScatPho[1][GunId];
	  
	   
	  if(GunId_A1 == GunId_S1 && GunId_A2 == GunId_S2) {
	    Relv4Hits.push_back(ScatPho[0]);
	    Relv4Hits.push_back(ScatPho[1]);	      
	  }

	  else if(GunId_A1 == GunId_S2 && GunId_A2 == GunId_S1) {
	    Relv4Hits.push_back(ScatPho[1]);
	    Relv4Hits.push_back(ScatPho[0]);	      
	  }

	  if(GunId_S1 != GunId_S2) {
	    Ana_Relv4Hits.push_back(ScatPho[0]);
	    Ana_Relv4Hits.push_back(ScatPho[1]);
	  }
	    
	  
	} //nAnniPho condition end

	if (Relv4Hits.size() == 4){
	 
	  float dt_A1S1 = Relv4Hits[S1][time] - Relv4Hits[A1][time];
	  float dt_A1S2 = Relv4Hits[S2][time] - Relv4Hits[A1][time];
	  float dt_A2S1 = Relv4Hits[S1][time] - Relv4Hits[A2][time];
	  float dt_A2S2 = Relv4Hits[S2][time] - Relv4Hits[A2][time];

	  float dt_A1A2 = Relv4Hits[A2][time] - Relv4Hits[A1][time];
	  float dt_S1S2 = Relv4Hits[S2][time] - Relv4Hits[S1][time];

	  h_dt_A1S1->Fill(dt_A1S1);
	  h_dt_A1S2->Fill(dt_A1S2);
	  h_dt_A2S1->Fill(abs(dt_A2S1));
	  h_dt_A2S2->Fill(dt_A2S2);
	  h_dt_A1A2->Fill(dt_A1A2);
	  h_dt_S1S2->Fill(dt_S1S2);	  
	}

	//if (Ana_Relv4Hits.size() == 4 && Ana_Relv4Hits.size() == 4){
	if (Ana_Relv4Hits.size() == 4){	  	  
	  // float Ana_dt_A1S1 = Ana_Relv4Hits[S1][time] - Ana_Relv4Hits[A1][time];	  
	  // float Ana_dt_A1S2 = Ana_Relv4Hits[S2][time] - Ana_Relv4Hits[A1][time];
	  // float Ana_dt_A2S1 = Ana_Relv4Hits[S1][time] - Ana_Relv4Hits[A2][time];
	  // float Ana_dt_A2S2 = Ana_Relv4Hits[S2][time] - Ana_Relv4Hits[A2][time];

	  // float Ana_dt_A1A2 = Ana_Relv4Hits[A2][time] - Ana_Relv4Hits[A1][time];
	  // float Ana_dt_S1S2 = Ana_Relv4Hits[S2][time] - Ana_Relv4Hits[S1][time];

	  float Ana_dt_A1S1 = Relv4Hits[S1][time] - Relv4Hits[A1][time];	  
	  float Ana_dt_A1S2 = Relv4Hits[S2][time] - Relv4Hits[A1][time];
	  float Ana_dt_A2S1 = Relv4Hits[S1][time] - Relv4Hits[A2][time];
	  float Ana_dt_A2S2 = Relv4Hits[S2][time] - Relv4Hits[A2][time];

	  float Ana_dt_A1A2 = Relv4Hits[A2][time] - Relv4Hits[A1][time];
	  float Ana_dt_S1S2 = Relv4Hits[S2][time] - Relv4Hits[S1][time];

	  //creating STij
	  ROOT::Math::XYZVector PosVec_A1, PosVec_A2, PosVec_S1, PosVec_S2, test;
	  // PosVec_A1.SetXYZ(Ana_Relv4Hits[A1][PosX], Ana_Relv4Hits[A1][PosY], Ana_Relv4Hits[A1][PosZ]);
	  // PosVec_A2.SetXYZ(Ana_Relv4Hits[A2][PosX], Ana_Relv4Hits[A2][PosY], Ana_Relv4Hits[A2][PosZ]);
	  // PosVec_S1.SetXYZ(Ana_Relv4Hits[S1][PosX], Ana_Relv4Hits[S1][PosY], Ana_Relv4Hits[S1][PosZ]);
	  // PosVec_S2.SetXYZ(Relv4Hits[S2][PosX], Relv4Hits[S2][PosY], Relv4Hits[S2][PosZ]);

	  PosVec_A1.SetXYZ(Relv4Hits[A1][PosX], Relv4Hits[A1][PosY], Relv4Hits[A1][PosZ]);
	  PosVec_A2.SetXYZ(Relv4Hits[A2][PosX], Relv4Hits[A2][PosY], Relv4Hits[A2][PosZ]);
	  PosVec_S1.SetXYZ(Relv4Hits[S1][PosX], Relv4Hits[S1][PosY], Relv4Hits[S1][PosZ]);
	  PosVec_S2.SetXYZ(Relv4Hits[S2][PosX], Relv4Hits[S2][PosY], Relv4Hits[S2][PosZ]);

	  double d_A1S1, d_A1S2, d_A2S1, d_A2S2;
	  d_A1S1 = (PosVec_A1 - PosVec_S1).R();
	  d_A1S2 = (PosVec_A1 - PosVec_S2).R();
	  d_A2S1 = (PosVec_A2 - PosVec_S1).R();
	  d_A2S2 = (PosVec_A2 - PosVec_S2).R();
	  

	  double ST11, ST12, ST21, ST22;
	  double c = TMath::C() * (1e-6);    //converting from m/s to mm/ns
	  

	  ST11 = Ana_dt_A1S1 - (d_A1S1/c);
	  ST12 = Ana_dt_A1S2 - (d_A1S2/c);
	  ST21 = Ana_dt_A2S1 - (d_A2S1/c);
	  ST22 = Ana_dt_A2S2 - (d_A2S2/c);

	  //h_rough->Fill(ST21);
	 
	  bool A1toS1 = false, A1toS2 = false;

	  // if (abs(ST11 < ST21)) {A1toS1 = true;}
	  // else if (abs(ST11 > ST21)) {A1toS2 = true;}
	  //if (abs(ST11) > abs(ST12)) cout << "fsaldkjkf" << endl; 

	  if (abs(ST11) < abs(ST12)) {A1toS1 = true;}
	  else if (abs(ST11) > abs(ST12))  {A1toS2 = true;}

	  
	  // if (A1toS1) {
	  //   AS_ConnHits.push_back(Ana_Relv4Hits[S1]);
	  //   AS_ConnHits.push_back(Ana_Relv4Hits[S2]);
	  // }

	  // else if (A1toS2) {
	  //   AS_ConnHits.push_back(Ana_Relv4Hits[S2]);
	  //   AS_ConnHits.push_back(Ana_Relv4Hits[S1]);
	  // }

	  if (AS_ConnHits.size() != 2) cout << "\n\ncaution\n\n" << endl;
	

	  if (A1toS1) {
	    AS_ConnHits.push_back(Relv4Hits[S1]);
	    AS_ConnHits.push_back(Relv4Hits[S2]);

	    //if(abs(AS_ConnHits[2][time] - Relv4Hits[2][time]) > 1) cout << "hfsakll" << endl; 
	  }

	  else if (A1toS2) {
	    AS_ConnHits.push_back(Relv4Hits[S2]);
	    AS_ConnHits.push_back(Relv4Hits[S1]);
	  }

	  //if (abs(ST11 ST21)) cout << "" << endl;

	  
	  h_rough->Fill(AS_ConnHits[S1][time]-Relv4Hits[S1][time]);
	  //h_rough->Fill(ST11);
	  //h_rough->Fill(Relv4Hits[A1][time]-Ana_Relv4Hits[A1][time]);

	  h_Ana_dt_A1S1->Fill(Ana_dt_A1S1);
	  h_Ana_dt_A1S2->Fill(Ana_dt_A1S2);
	  h_Ana_dt_A2S1->Fill(Ana_dt_A2S1);
	  h_Ana_dt_A2S2->Fill(Ana_dt_A2S2);
	  h_Ana_dt_A1A2->Fill(Ana_dt_A1A2);
	  h_Ana_dt_S1S2->Fill(abs(Ana_dt_S1S2));		  
	}
	  
      }       // nHit >0 condition brkt end
    } //jentry loop end

  cout << "events with 2 or more comps: " << Comp2_evts << endl;
  cout << "n1: " << n1 << endl;
  cout << "n1s: " << n1s << endl;
  cout << "n1b: " << n1b << endl;
}

void AnalyzeLightBSM::FillHistogram(int cat, double time, double posX, double posY, double posZ, double theta_deg, double eta_deg, double pol0, double pol1, double pol2, double Ein, double Eout){
  h_hitTime[cat]->Fill(time);
  h_hitPosX[cat]->Fill(posX);
  h_hitPosY[cat]->Fill(posY);
  h_hitPosZ[cat]->Fill(posZ);

  h_theta[cat]->Fill(theta_deg);

  h_eta[cat]->Fill(eta_deg);	    
  h_theta_eta[cat]->Fill(eta_deg, theta_deg);

  h_pol0[cat]->Fill(pol0);
  h_pol1[cat]->Fill(pol1);
  h_pol2[cat]->Fill(pol2);

  h_Ein[cat]->Fill(Ein);
  h_Eout[cat]->Fill(Eout);
  
  h_theta_Eout[cat]->Fill(theta_deg, Eout);
  
  if(theta_deg > 50 && theta_deg < 100) h_compSigEta[cat]->Fill(eta_deg);    //taking large values of theta to see maximum variation
}

string AnalyzeLightBSM::MatName(int DetId){
  string mat = " ";
  if (DetId < 0) mat = "error";
  if (DetId /4000 >=1) mat = "BGO"; 
  else if (DetId/ 1000 >=1 && DetId/ 4000 <1) mat = "Scint";

  return mat;
}

// double distance(vector<double> v1, vector<double> v2){
//   double ds;
//   //ds = sqrt()
// }

// vector<float> AnalyzeLightBSM::SortInTime(vector<float> v){
//   vector<float> v_cpy, v_sort;
//   v_cpy = v;
  
// }




















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
	   
