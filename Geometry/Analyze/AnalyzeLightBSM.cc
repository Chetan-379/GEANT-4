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
  const char *detType          = argv[3];

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
   
    for (Long64_t jentry=0; jentry<nentries;jentry++)
      {
	double progress = 10.0 * jentry / (1.0 * nentries);
	int k = int (progress);
	if (k > decade)
	  cout << 10 * k << " %" << endl;
	decade = k;



// 	// ===============read this entry == == == == == == == == == == ==                                                                        
 	Long64_t ientry = LoadTree(jentry);
	
	if(Debug)
	  cout<<"===load tree entry ==="<<jentry<<endl;
	if (ientry < 0) break;


 	nb = fChain->GetEntry(jentry);   nbytes += nb;
		
	//cout << "\n\nh_theta: " << h_theta << endl;
	if(Debug)
	  cout<<"===load tree entry ==="<<jentry<<endl;

	// do something with tree branches here	
	h_Compt_Edep-> Fill(Edep_Compt);
	h_Photo_Edep-> Fill(Edep_Photo);
	h_ComptVsPhoto_Edep->Fill(Edep_Compt, Edep_Photo);
	h_Total_Edep-> Fill(Total_Edep);

	nHits = Hit_Time->size();
	h_nHits->Fill(nHits);

	//ROOT::Math::XYZVector v_inc(0,0,0);
	if(Hit_Time->size() >1) {
	  //cout << "\n============Event: " << jentry << endl;
	  Comp2_evts++;

	   ROOT::Math::XYZVector v_comp, vin, vout, v1, v2;

	   //sorting the hits in the ascending order of time
	   vector<double> HitPosX_sort, HitPosY_sort, HitPosZ_sort, HitTime_sort, HitScatAng_sort, HitEta_sort;
	   vector<double> HitPosX_cpy, HitPosY_cpy, HitPosZ_cpy, HitTime_cpy, HitScatAng_cpy, HitEta_cpy;
	   vector <int> HitProcId_sort;
	   vector<int> HitProcId_cpy;
	   
	   HitPosX_sort.clear(); HitPosY_sort.clear(); HitPosZ_sort.clear(); HitTime_sort.clear(), HitScatAng_sort.clear(), HitEta_sort.clear(), HitProcId_sort.clear();
	   HitPosX_cpy.clear(); HitPosY_cpy.clear(); HitPosZ_cpy.clear(); HitTime_cpy.clear(), HitScatAng_cpy.clear(), HitEta_cpy.clear(), HitProcId_cpy.clear(); 
	   
	   //if (jentry == 159){
	     HitPosX_cpy = *Hit_PositionX;
	     HitPosY_cpy = *Hit_PositionY;
	     HitPosZ_cpy = *Hit_PositionZ;
	     HitTime_cpy = *Hit_Time;
	     HitScatAng_cpy = *Hit_ScatAngle;
	     HitEta_cpy = *Hit_Eta;
	     HitProcId_cpy = *Hit_ProcId;
	     
	     for (int iHit=0 ; iHit< nHits; iHit++){                               	 
	     auto min_id_ptr = std::min_element(HitTime_cpy.begin(), HitTime_cpy.end());
	     size_t min_idx = std::distance(HitTime_cpy.begin(), min_id_ptr);

	     //cout << "type max time: " << typeid(max_time).name() << endl;
	     
	     // cout << "min time: " << *min_id_ptr << endl;
	     // cout << "min idx: " << min_idx << endl;

	     HitPosX_sort.push_back(HitPosX_cpy[min_idx]);
	     HitPosY_sort.push_back(HitPosY_cpy[min_idx]);
	     HitPosZ_sort.push_back(HitPosZ_cpy[min_idx]);
	     HitTime_sort.push_back(HitTime_cpy[min_idx]);
	     HitScatAng_sort.push_back(HitScatAng_cpy[min_idx]);
	     HitEta_sort.push_back(HitEta_cpy[min_idx]);
	     HitProcId_sort.push_back(HitProcId_cpy[min_idx]);

	     //cout << "original size of Hit vector: " << HitTime_cpy.size() << endl;
	      HitTime_cpy.erase(HitTime_cpy.begin() + min_idx);
	      HitPosX_cpy.erase(HitPosX_cpy.begin() + min_idx);
	      HitPosY_cpy.erase(HitPosY_cpy.begin() + min_idx);
	      HitPosZ_cpy.erase(HitPosZ_cpy.begin() + min_idx);
	      HitScatAng_cpy.erase(HitScatAng_cpy.begin() + min_idx);
	      HitEta_cpy.erase(HitEta_cpy.begin() + min_idx);
	      HitProcId_cpy.erase(HitProcId_cpy.begin() + min_idx);
	      
	      //cout << "size after cutting: " << HitTime_cpy.size() << "\n\n";	     
	     }

	     //counting the events with atleast one rayleigh
	     for (int i = 0; i< HitTime_sort.size(); i++){
	       if ((*Hit_ProcId)[i] == 0) { nRayl++; break;}    //procId 0 for rayleigh
	     }      

	   for (int iHit=0 ; iHit< nHits-1; iHit++){     //cannot calculate theta at the last hit
	     
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
	    double scat_theta_deg = scat_theta*180.0 / TMath::Pi();
	    
	    double scat_theta_SD = HitScatAng_sort[iHit];
	    //cout << "\n\nWorkgin till here!! " << "\n\n" << endl;
	    double Ang_diff = scat_theta - scat_theta_SD;
	    h_ScatAng_SDvsAna->Fill(scat_theta, scat_theta_SD);

	    // if (Ang_diff > 0.0001) {
	    //   nDiffAngle++;
	    //   //cout << "spurious event: " << jentry+1 << endl;}

	    double EtaPolScat = HitEta_sort[iHit];
	    double EtaPolScat_deg = EtaPolScat*180.0 / TMath::Pi();
	    int procId = HitProcId_sort[iHit];
	    
	    
	    //if (abs(Ang_diff) > 0.5) cout << "EventID of diff angle: " << jentry << endl;
	    
	    h_diff_Ana_SD->Fill(Ang_diff);
	    
	    if (procId !=0 && procId !=1)       //excluding rayleigh && phot
	      {
	        h_theta->Fill(scat_theta_deg);
		h_eta->Fill(EtaPolScat_deg);
		h_theta_eta->Fill(EtaPolScat_deg, scat_theta_deg);
		
		if(scat_theta_deg > 50 && scat_theta_deg < 100) h_compSigEta->Fill(EtaPolScat*180.0 / TMath::Pi());      //taking large values of theta to see more variation
	      }	
	     
	   }   //end of iHit loop
	}      // nHit >1 condition brkt end
	
	
      } //jentry loop end

    cout << "events with 2 or more comps: " << Comp2_evts << endl;
    cout << "no. of hits with non zero diff in theta are: " << nDiffAngle << endl;
    cout << "no. of events with atleast one rayleigh: " << nRayl << endl;

 }

//==================================================================================================================================================================================================================












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
	   
