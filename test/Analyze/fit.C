
double fitf(double *x,double *par) {
  double fitval = (par[0]*TMath::Exp(par[2]*x[0])) + (par[1]*TMath::Exp(par[3]*x[0]));
  return fitval;
}


void fit(){
  TCanvas *canvas_n1 = new TCanvas("tag_name", "tag_name",900,750);//600,600,1200,1200);
  canvas_n1->Range(-60.25,-0.625,562.25,0.625);
  canvas_n1->SetFillColor(0);
  canvas_n1->SetBorderMode(0);
  canvas_n1->SetBorderSize(2);
  gStyle->SetOptStat(0);

  
  TF1 *func = new TF1("fit",fitf,1,40,4);
  
  TFile *f = new TFile("out_root_files/Opt_Pho_SY200_SY1_0.8_SY2_0.2_Emission_time_PbWO4_out.root");
  
  string hist_name = "OptPho_time"; 
  TH1F *hist = (TH1F*)f->Get(hist_name.c_str());

  cout << hist->GetName() << endl;
  
  //func->SetParameters(500,hpx->GetMean(),hpx->GetRMS());
  func->SetParameters(0.8, 0.2, -1, -10);//, 0.1, -20);
  //func->SetParameters(1, -50); //, -10, -1);
  //cout << "working" << endl;
  // Give the parameters names.
  //func->SetParNames("Constant","Mean_value","Sigma");
  func->SetParNames("A", "B", "decay Const1", "decay Const2");//, "C", "decay Const3");
  //func->SetParNames("A","decay Const1");
  
  // Call TH1::Fit with the name of the TF1 object.
  hist->Fit("fit");
  canvas_n1->Draw();
  hist->Draw("hist");
  func->Draw("SAME");
  //f->Close();
}
