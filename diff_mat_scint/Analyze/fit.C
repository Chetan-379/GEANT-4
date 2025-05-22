double fitf(double *x,double *par) {
  double fitval = (par[0]*TMath::Exp(par[2]*x[0])) + (par[1]*TMath::Exp(par[3]*x[0]));
  return fitval;
}


void fit(string rootfile){
  TCanvas *canvas_n1 = new TCanvas("tag_name", "tag_name",900,750);
  canvas_n1->Range(-60.25,-0.625,562.25,0.625);
  canvas_n1->SetFillColor(0);
  canvas_n1->SetBorderMode(0);
  canvas_n1->SetBorderSize(2);
  gStyle->SetOptStat(1);

  TFile *f = new TFile(rootfile.c_str());
  
  string hist_name = "OptPho_time"; 
  TH1F *hist = (TH1F*)f->Get(hist_name.c_str());
  hist->GetXaxis()->SetTitle("time (ns)");
  cout << hist->GetName() << endl;

  //TPaveStats *ptstats1 = new TPaveStats(0.50,0.50,0.965,0.9,"brNDC");
  //gPad->Update();

  

  TF1 *func = new TF1("fit",fitf,1,100,4); 
  func->SetParameters(0.2, 0.8, -10, -1); 
  func->SetParNames("A", "B", "decay Const1", "decay Const2");
  //func->SetStat(0);
  
  hist->Fit("fit","R");

  canvas_n1->Draw();
  
  hist->Draw("hist");
  func->Draw("same");

  TPaveStats *ptstats1 = (TPaveStats*) hist->FindObject("stats");
  if (ptstats1){
    ptstats1->SetX1NDC(0.68);
    ptstats1->SetY1NDC(0.75);
    ptstats1->SetX2NDC(0.9);
    ptstats1->SetY2NDC(0.9);
    ptstats1->SetBorderSize(1);
    ptstats1->SetFillColor(0);
    ptstats1->SetLineColor(hist->GetLineColor());
    ptstats1->SetTextAlign(12);
    ptstats1->SetTextColor(kBlue);
    ptstats1->SetTextFont(42);
    ptstats1->SetTextSize(0.03);
    ptstats1->SetOptStat(110011);
    //ptstats1->SetOptFit(1);
    hist->GetListOfFunctions()->Add(ptstats1);
    ptstats1->SetParent(hist);
  }
  canvas_n1->Update();  
  
  TLatex *latex = new TLatex();
  latex->SetTextSize(0.03);
  // char* fitted_func = new char[500];
  // sprintf(fitted_func,"%1f exp(%fx) + %.1fexp(%fx)",func->GetParameter(0),func->GetParameter(2),func->GetParameter(1),func->GetParameter(3));
  char* param1 = new char[50];
  char* param2 = new char[50];
  char* param3 = new char[50];
  char* param4 = new char[50];
  char* chi2ndf = new char[50];

  float delTau1=(1/pow(func->GetParameter(3),2))*func->GetParError(3);
  float delTau2=(1/pow(func->GetParameter(2),2))*func->GetParError(2);
  
  sprintf(param1, "A: %0.2e #pm %0.2f", func->GetParameter(1), func->GetParError(1));
  sprintf(param2, "B: %0.2e #pm %0.2f", func->GetParameter(0), func->GetParError(0));
  sprintf(param3, "tau1: %0.2f #pm %0.2e", -1/func->GetParameter(3), delTau1);
  sprintf(param4, "tau2: %0.2f #pm %0.2e", -1/func->GetParameter(2), delTau2);
  sprintf(chi2ndf, "\\chi^{2}/NDF: %0.2e", func->GetChisquare()/func->GetNDF());

  latex->DrawLatexNDC(0.66,0.70, param1);
  latex->DrawLatexNDC(0.66,0.65, param2);
  latex->DrawLatexNDC(0.66,0.60, param3);
  latex->DrawLatexNDC(0.66,0.55, param4);
  latex->DrawLatexNDC(0.66,0.50, chi2ndf);  

 
  char* canvas_name = new char[100];
  sprintf(canvas_name,"%s.png","plots/Scintillation/Opt_Pho_Ana/PbWO4");
  canvas_n1->SaveAs(canvas_name);

  gPad->SetLogy();
  //gPad->Update();
  sprintf(canvas_name,"%s_logY.png","plots/fitting");
  canvas_n1->SaveAs(canvas_name);
  
}


