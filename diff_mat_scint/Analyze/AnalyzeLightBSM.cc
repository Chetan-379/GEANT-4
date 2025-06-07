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
  
    for (Long64_t jentry=0; jentry<nentries;jentry++)
      {
	double progress = 10.0 * jentry / (1.0 * nentries);
	int k = int (progress);
	if (k > decade)
	  cout << 10 * k << " %" << endl;
	decade = k;
	
	// ===============read this entry == == == == == == == == == == ==                                                                        
	Long64_t ientry = LoadTree(jentry);
	if(Debug)
	  cout<<"===load tree entry ==="<<jentry<<endl;
	if (ientry < 0) break;
	nb = fChain->GetEntry(jentry);   nbytes += nb;
	if(Debug)
	  cout<<"===load tree entry ==="<<jentry<<endl;

	  // do something with tree branches here

	h_Compt_Edep-> Fill(Edep_Compt);
	h_Photo_Edep-> Fill(Edep_Photo);
	h_ComptVsPhoto_Edep->Fill(Edep_Compt, Edep_Photo);
	h_Total_Edep-> Fill(Total_Edep);
	//h_Total_Edep_fine_binned->Fill(Total_Edep);

	for (int i=0; i< OptPho_PosX->size(); i++){
	  h_OptPho_PosX->Fill((*OptPho_PosX)[i]);
	  h_OptPho_PosY->Fill((*OptPho_PosY)[i]);
	  h_OptPho_PosZ->Fill((*OptPho_PosZ)[i]);

	// if ((*OptPho_PosZ)[i] == 150) h_OptPho_XvsY->Fill((*OptPho_PosX)[i], (*OptPho_PosY)[i]);
	 }

	h_nOptPho-> Fill(nOpticalPhotons);
	//=====
	// h_nOptPho-> Fill(nOpticalPhotons);
	// h_nOptPho_Edep-> Fill(Total_Edep, nOpticalPhotons);
	
	
	// int nOptPhotOnDet = 0;
	// for (int i =0; i< positionX->size(); i++){
	//   if (abs((*OptPho_PosX)[i]) < 40 && abs((*OptPho_PosY)[i]) < 40 && (*OptPho_PosZ)[i] == 150.) nOptPhotOnDet++;
	//   //cout << 1240e-6/(*OptPho_Energy)[i] << endl;	  
	// }

	for (int i =0; i< OptPho_Energy->size(); i++){
	  h_OptPho_lmbda-> Fill((1240e-6)/(*OptPho_Energy)[i]);
	  //h_OptPho_lmbda-> Fill((*OptPho_Energy)[i]);
	  h_OptPho_time-> Fill((*OptPho_Time)[i]);
	}

	// h_OptPhoOnDet -> Fill(nOptPhotOnDet);
	// h_nOptPhoOnDet_genOptPho-> Fill(nOptPhotOnDet, nOpticalPhotons);
	//=============

	double Einc = 1.;
	if (Total_Edep > (Einc-0.011)) Full_edep_evts++;
	//if (jentry < 100)cout << "Edep: " << Total_Edep << endl; 
      }

    cout << "Events with full Edep: " << Full_edep_evts << endl;
    //cout << "total Integral is: " << h_Total_Edep->Integral() << endl;
    cout << "Efficiency is: " << Full_edep_evts/(h_Total_Edep->Integral()) << endl; 
}
