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

    //Analysing truth info
    float Edep_truth_cryst = 0;
    int clr_Idx = -99;
    for(int i =0; i<Edep_truth_vec->size(); i++){
      Edep_truth_cryst = (*Edep_truth_vec)[i];
      clr_Idx = (*Det_clr_Idx)[i];

      h_Edep_truth_cryst->Fill(Edep_truth_cryst);

      if(clr_Idx ==1) h_EdepG_truth_cryst->Fill(Edep_truth_cryst);
      if(clr_Idx ==2) h_EdepL_truth_cryst->Fill(Edep_truth_cryst);	
    }

  } // jentry loop end
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
