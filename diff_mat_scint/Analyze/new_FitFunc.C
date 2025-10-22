int line_width[12] = {2,2,2,2,2,2,2,2,2,2,2,2};
int line_style[12] = {1,1,1,1,1,1,1,1,1,1,1,1};
int marker_style[12] = {8,22,34,47,33,22,34,39,41,43,45,23};
int line_color[9] = {kBlue,kRed,kGreen+2,kViolet+2,kCyan+2,kYellow+1,kGray+2,kMagenta,kBlue+2};
int line_color1[9]= {kBlue,kGreen+2,kGray+1,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta,kBlue+2};
int line_color2[9] = {kGreen+2,kBlue,kViolet,kGray,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta};
vector<int> col={kGreen+2,kBlue,kViolet,kGray,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta,kBlue+2,kMagenta,kCyan};
vector<int> Style={3008,1001,3008,1001};
void decorate(TH1F*,int);
void decorate(TH1F* hist,int i){
  hist->SetLineWidth(3);
}
void setLastBinAsOverFlow(TH1F*);
TH1F* setMyRange(TH1F*,double,double);

// TH1F* setMyRange(TH1F *h1,double xLow,double xHigh){
//   //call it after setting last bin as overflow
//   double err=0;
//   if(xHigh > 13000) return h1;
//   if(xLow < -13000) return h1;
//   int nMax=h1->FindBin(xHigh);
//   h1->SetBinContent(nMax,h1->IntegralAndError(nMax,h1->GetNbinsX(),err));
//   h1->SetBinError(nMax,err);
//   for(int i=nMax+1;i<=h1->GetNbinsX()+1;i++){
//     h1->SetBinContent(i,0);
//     h1->SetBinError(i,0);

//   }
//   h1->GetXaxis()->SetRangeUser(xLow,xHigh);
//   cout<<xLow<<"\t"<<xHigh<<"\t"<<"set range"<<endl;
//   return h1;
// }

TH1F* DrawOverflow(TH1F*);
TH1F* DrawOverflow(TH1F* h,int xmin, int xrange){

  UInt_t nx    = h->GetNbinsX()+1;
  Double_t *xbins= new Double_t[nx+1];
  for (UInt_t i=0;i<nx;i++)
    xbins[i]=h->GetBinLowEdge(i+1);
  xbins[nx]=xbins[nx-1]+h->GetBinWidth(nx);
  char *tempName= new char[strlen(h->GetName())+10];
  sprintf(tempName,"%swtOverFlow",h->GetName());
  h->GetXaxis()->SetLimits(xmin,xrange);
  // Book a temporary histogram having ab extra bin for overflows
  TH1F *htmp = new TH1F(tempName, h->GetTitle(), nx, xbins);
  htmp->GetXaxis()->SetRange(xmin,xrange);
  // Reset the axis labels
  htmp->SetXTitle(h->GetXaxis()->GetTitle());
  htmp->SetYTitle(h->GetYaxis()->GetTitle());
  // Fill the new hitogram including the extra bin for overflows
  for (UInt_t i=1; i<=nx; i++)
    htmp->Fill(htmp->GetBinCenter(i), h->GetBinContent(i));
  // Fill the underflows
  htmp->Fill(h->GetBinLowEdge(1)-1, h->GetBinContent(0));

  htmp->SetEntries(h->GetEntries());

  return htmp;
}


void setLastBinAsOverFlow(TH1F* h_hist){
  double lastBinCt =h_hist->GetBinContent(h_hist->GetNbinsX()),overflCt =h_hist->GetBinContent(h_hist->GetNbinsX()+1);
  double lastBinErr=h_hist->GetBinError(h_hist->GetNbinsX()),  overflErr=h_hist->GetBinError(h_hist->GetNbinsX()+1);

  if(lastBinCt!=0 && overflCt!=0)
    lastBinErr = (lastBinCt+overflCt)* (sqrt( ((lastBinErr/lastBinCt)*(lastBinErr/lastBinCt)) + ((overflErr/overflCt)*(overflErr/overflCt)) ) );

  else if(lastBinCt==0 && overflCt!=0)
    lastBinErr = overflErr;
  else if(lastBinCt!=0 && overflCt==0)
    lastBinErr = lastBinErr;
  else lastBinErr=0;

  lastBinCt = lastBinCt+overflCt;
  h_hist->SetBinContent(h_hist->GetNbinsX(),lastBinCt);
  h_hist->SetBinError(h_hist->GetNbinsX(),lastBinErr);
  cout<<lastBinCt<<"\t"<<"Last bin values"<<endl;

}
#include "TF1.h"
#include "TMath.h"
#include <cmath>

// Define the Gaussian function
double Gauss(double *x, double *par) {
  double A = par[0];   // Amplitude
  double mu = par[1];  // Mean
  double sigma = par[2]; // Std. deviation
  double arg = ((x[0] - mu) * (x[0] - mu)) / (sigma*sigma);
  
  return A * TMath::Exp(-0.5 * arg);
}

const int nfiles=100;
TFile *f[nfiles];

vector<double> generate_1Dplot(vector<TH1*> hist, char const *tag_name="", int xmax=-1,int xmin=-1,char const *leg_head="",
			       bool normalize=false, bool log_flag=true, bool DoRebin=false, bool save_canvas=true, char const *title="", const char* xtitile="", const char* ytitile="", int rebin=-1){

  TCanvas *canvas_n1 = new TCanvas(tag_name, tag_name, 950, 850);
  canvas_n1->SetLeftMargin(0.135);
  canvas_n1->SetRightMargin(0.035);
  canvas_n1->SetTopMargin(0.067);
  canvas_n1->SetBottomMargin(0.1);
  gStyle->SetOptStat(0);

  TLegend *legend = new TLegend(0.66, 0.84, 0.80, 0.90);
  legend->SetTextSize(0.035);
  legend->SetLineColor(kWhite);
  legend->SetNColumns(2);
  legend->SetHeader(title);
  double MaxY=0;

  TH1F* hist_uncut;
  
  for(int i = 0; i < (int)hist.size(); i++) {
    hist_uncut = (TH1F*)hist.at(i)->Clone("hist_uncut");
    if(normalize)   hist.at(i)->Scale(1.0 / hist.at(i)->Integral());
    hist.at(i)->Rebin(rebin);
    hist.at(i)->GetXaxis()->SetRangeUser(xmin,xmax);

    if(hist.at(i)->GetMaximum() > MaxY){MaxY=hist.at(i)->GetMaximum();}
  }
  
  vector<vector<TLatex>> latexVec(hist.size(),vector<TLatex>(3));
  vector<TF1*> Fit_func(hist.size());
  vector<double> results;
  
  for(int i = 0; i < (int)hist.size(); i++) {         
    if(normalize) {
      hist.at(i)->Scale(1.0 / hist.at(i)->Integral());
      hist.at(i)->GetYaxis()->SetTitle("Normalized");
    }

    else {
      hist.at(i)->GetYaxis()->SetTitle("Entries");
    }

    hist.at(i)->GetXaxis()->SetTitle("nOpticalPhotons");

    // Skiping the peak at 0
    int maxBin = -1;
    double maxContent = -1;
    for (int j = 1; j <= hist.at(i)->GetNbinsX(); ++j) {
      double center = hist.at(i)->GetBinCenter(j);
      double content = hist.at(i)->GetBinContent(j);

      if (center >= -0.5 && center <= +0.5) continue;

      if (content > maxContent) {
	maxContent = content;
	maxBin = j;
      }
    }

    int peakBin = maxBin;
    //double peakValue = hist.at(i)->GetXaxis()->GetBinCenter(peakBin);
    double peakValue = hist.at(i)->GetBinContent(hist.at(i)->GetMaximumBin());
    
    // double fitMin = hist.at(i)->GetXaxis()->GetBinCenter(binLeft);
    // double fitMax = hist.at(i)->GetXaxis()->GetBinCenter(binRight);
    
    double fitMin = 0, fitMax = 0;
    double old_mean=hist.at(i)->GetMean();
    double old_sigma = hist.at(i)->GetRMS();

    fitMin = old_mean - 3*old_sigma;
    fitMax = old_mean + 3*old_sigma;


    //Doing Iterative Fitting
    for (int itr = 0; itr < 400; ++itr) {
      Fit_func[i] = new TF1("GaussFit", Gauss, fitMin, fitMax, 3);
      Fit_func[i]->SetParameters(peakValue, old_mean, old_sigma);
      hist.at(i)->Fit(Fit_func[i], "ERQ");

      double mu = Fit_func[i]->GetParameter(1);
      double sigma = Fit_func[i]->GetParameter(2);
      double chi2byNDF = Fit_func[i]->GetChisquare()/Fit_func[i]->GetNDF();

      double rangeMin = mu - 2.0 * abs(sigma);
      double rangeMax = mu + 2.0 * abs(sigma);

      if (fabs(fitMin-rangeMin)/rangeMin <= 0.0002) break;
      fitMin = rangeMin;
      fitMax = rangeMax;

      Fit_func[i]->SetRange(fitMin, fitMax);

      old_mean = mu;
      old_sigma = abs(sigma);
    }

    hist.at(i)->GetXaxis()->SetRangeUser(xmin, xmax);
    hist.at(i)->Draw("hist");
    Fit_func[i]->Draw("SAME");

    double A, Mu, Sigma;

    A = Fit_func[i]->GetParameter(0);
    Mu = Fit_func[i]->GetParameter(1);
    Sigma = Fit_func[i]->GetParameter(2);

    TLatex *latex = new TLatex();
    latex->SetTextSize(0.03);

    char* param1 = new char[50];
    char* param2 = new char[50];
    char* param3 = new char[50];
    char* chisqrNDF = new char[50];

    sprintf(param1, "A: %0.2e #pm %0.4f", A, Fit_func[i]->GetParError(0));
    sprintf(param2, "\\mu: %0.3f \\pm %0.4f", Mu, Fit_func[i]->GetParError(1));
    sprintf(param3, "\\sigma: %0.4f \\pm %0.4f", abs(Sigma), Fit_func[i]->GetParError(2));
    sprintf(chisqrNDF, "\\chi^{2}/NDF: %.4f", Fit_func[i]->GetChisquare()/Fit_func[i]->GetNDF());

    latex->DrawLatexNDC(0.70,0.80, param1);
    latex->DrawLatexNDC(0.70,0.75, param2);
    latex->DrawLatexNDC(0.70,0.70, param3);
    latex->DrawLatexNDC(0.70,0.65, chisqrNDF);

    
    //calculating the efficiency
    int BinMax=0, BinMin=0;
    double MAX=0, MIN=0;
    double GaussArea=0, HistArea = 0, Efficiency =0;

    MAX = Mu + 3 * abs(Sigma);
    MIN = Mu - 3 * abs(Sigma);

    BinMax = hist_uncut->FindBin(MAX);
    BinMin = hist_uncut->FindBin(MIN);

    cout << "\n\n" << tag_name << endl;
    cout << "xMax: " << MAX << endl;
    cout << "xMin: " << MIN << endl;

    cout << "\nbinMax: " << BinMax << endl;
    cout << "binMin: " << BinMin << endl;

    HistArea = hist_uncut->Integral();

    //calculating the integral under Gaussian from mu-3sigma to mu+3sigma
    for (int ibin = BinMin; ibin <= BinMax; ibin++){
      GaussArea += hist_uncut->GetBinContent(ibin);
    }

    Efficiency = GaussArea/HistArea;
    cout << "\nGaussian Area is: " << GaussArea << endl;
    cout << "Hist Area is: " << HistArea << endl;
    cout << "Efficiency is: " << Efficiency << endl;

    results.push_back(Mu);
    results.push_back(Fit_func[i]->GetParError(1));
    results.push_back(Sigma);
    results.push_back(Fit_func[i]->GetParError(2));
    results.push_back(Efficiency);
    
  }  //overlay hist loop end

  legend->Draw();
  gPad->Modified();
  gPad->Update();

  if (save_canvas) canvas_n1->SaveAs(Form("%s.png", tag_name));
  //if (save_canvas) canvas_n1->SaveAs(Form("%s.pdf", tag_name));

  canvas_n1->Close();
  
  return results;
}

// Defining the fit function
double sqrtE(double *x, double *par) {

  double a = par[0];   // coefficient of 1/sqrtE term
  double b = par[1];   // coefficient of 1/E term
  double c = par[2];   // constant term

  return (sqrt((pow(a/sqrt(x[0]), 2)) + (pow(b/x[0],2)) + (pow(c,2))));
}

double Linear(double *x, double *par) {
  double a = par[0];   // slope  
  return a*x[0];
}

void FitOverlay(vector<TF1*> funcs, vector<TGraphErrors*> graphs, string MatName, string VarName) {
  TCanvas* c1 = new TCanvas("c1", "Resolution fit Overlay", 800, 600);
    
  funcs[0]->SetMinimum(0.07); // set lower Y bound
  //funcs[0]->SetMaximum(0.045);  // set upper Y bound
 
  for (int i = 0; i < funcs.size(); ++i) {
    TString name = Form("f%d", i);
    funcs[i]->SetLineColor(i + 1);
    funcs[i]->SetLineWidth(2);
    funcs[i]->GetXaxis()->SetTitleSize(0.05);
    funcs[i]->GetXaxis()->SetTitleOffset(0.86);

    funcs[i]->GetYaxis()->SetTitleOffset(1.2);
    funcs[i]->GetYaxis()->SetTitleSize(0.05);   
  }

  funcs[0]->Draw();
  
  for (int i = 1; i < funcs.size(); ++i) {
    funcs[i]->Draw("SAME");
  }
  
  TH1* func_hist = funcs[0]->GetHistogram();
  func_hist->GetXaxis()->SetTitleSize(0.07);
  func_hist->GetYaxis()->SetTitleSize(0.07);

  c1->Update();

  auto legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  for (int i = 0; i < funcs.size(); ++i) {
    legend->AddEntry(funcs[i], funcs[i]->GetName(), "l");
  }

  legend->Draw();

  string OvrlImgName1;
  if(VarName == "Resol") {
    OvrlImgName1 = MatName + "_resolution_function_overlay.png";
    c1->SaveAs(OvrlImgName1.c_str());

    // OvrlImgName1 = MatName + "_resolution_function_overlay.pdf";
    // c1->SaveAs(OvrlImgName1.c_str());
  }
  
  c1->Close();
  
 
   
  TCanvas* c2 = new TCanvas("c2", "Resolution graphs Overlay", 800, 600);
  c2->SetLeftMargin(0.15); 
  c2->SetRightMargin(0.05);

  int colors[] = {kRed, kBlue, kGreen + 2, kBlack};
  for (int i = 0; i < graphs.size(); ++i)
    {
      graphs[i]->SetLineWidth(2);
      graphs[i]->SetMarkerStyle(22+i);
      graphs[i]->SetMarkerSize(1.5);
      graphs[i]->SetMarkerColor(colors[i % 4]);
      
      graphs[i]->GetXaxis()->SetTitleSize(0.05);
      graphs[i]->GetXaxis()->SetTitleOffset(0.86);
      
      if (VarName == "Resol") graphs[i]->GetYaxis()->SetTitleOffset(1.2);
      if (VarName == "Eff") graphs[i]->GetYaxis()->SetTitleOffset(0.8);
      
      graphs[i]->GetYaxis()->SetTitleSize(0.05);   
    }
  
  graphs[0]->Draw("AP");
  for (int i = 1; i < graphs.size(); ++i)
    {
      graphs[i]->Draw("P SAME");
    }
    
  TLegend* leg2;
  if (VarName == "Resol")  leg2 = new TLegend(0.7, 0.7, 0.9, 0.9);
  else if (VarName == "Eff") leg2 = new TLegend(0.2, 0.12, 0.4, 0.32);
  
  for (size_t i = 0; i < graphs.size(); ++i)
    leg2->AddEntry(graphs[i], graphs[i]->GetName(), "AP");
  leg2->Draw();

  string OvrlImgName2;
  if (VarName == "Resol") OvrlImgName2 = MatName + "_resolution_points_overlay.png";
  else if (VarName == "Eff") OvrlImgName2 = MatName + "_Efficiency_points_overlay.png";
  c2->SaveAs(OvrlImgName2.c_str());

  // OvrlImgName2 = MatName + "_resolution_points_overlay.pdf";
  // c2->SaveAs(OvrlImgName2.c_str());

  c2->Close();  
}

vector<TF1*> Fit_array;
vector<TGraphErrors*> Graph_array_Resol, Graph_array_Eff;

void data_plot(string MatName, string txtFilePath, string PlotsFolder) {
  vector<string> data_files;
  
  //filling all the txt files in an array
  TSystemDirectory dir(txtFilePath.c_str(), txtFilePath.c_str());
  TList* files = dir.GetListOfFiles();
  if (!files) return;

  TIter next(files);
  TObject* obj;

  while ((obj = next())) {
    TString fname = obj->GetName();
    if (!obj->IsA()->InheritsFrom("TSystemFile")) continue;
    if (fname.EndsWith(".txt")) {
      data_files.push_back(txtFilePath + string(fname.Data()));
    }
  }

  cout << "\n\nSize of the files vector is: " << data_files.size() << endl;

  for (int i=0; i<data_files.size(); i++){
    string xtitle, ytitle;    
    TString FileName = data_files[i];
    string ImgId;    
    ifstream infile;
    vector<double> x_vals, y_vals, x_err_vals, y_err_vals;
    double x, y, y_err;
    string WidthId;

    xtitle = "E_{incident}(keV)";
    
    if (FileName.Contains("1cm")) WidthId = "1cm";
    if (FileName.Contains("2.5cm")) WidthId = "2.5cm";
    if (FileName.Contains("_5cm")) WidthId = "5cm";
    if (FileName.Contains("_10cm")) WidthId = "10cm";
    

    // Create canvas
    data_files[i].resize(data_files[i].length() - 4);        //removing .txt part

    TCanvas *canvas = new TCanvas(data_files[i].c_str(), data_files[i].c_str(), 800, 600);
    canvas->SetLeftMargin(0.135);
    canvas->SetRightMargin(0.035);
    canvas->SetTopMargin(0.057);
    canvas->SetBottomMargin(0.11);
    
    TGraph *graph;
    if (FileName.Contains("Emeasured_vs_Einc") || FileName.Contains("Ngamma_vs_Einc")){
      if (FileName.Contains("Ngamma_vs_Einc")) {ImgId = "Ngamma_vs_Einc"; ytitle = "N_{\\gamma}";}
      if (FileName.Contains("Emeasured_vs_Einc")) {ImgId = "Emeasured_vs_Einc"; ytitle = "E_{Measured}(keV)";} 

      infile.open(FileName);
      while (infile >> x >> y){
        x_vals.push_back(x);
        y_vals.push_back(y);}

      if (x_vals.empty()) {
        cerr << "No data found!" << endl;
        return;}

      graph = new TGraph(x_vals.size(), &x_vals[0], &y_vals[0]);

      graph->SetTitle(" ");
      graph->SetMarkerStyle(20);
      graph->SetMarkerSize(1.0);
      graph->SetMarkerColor(kBlue);
      graph->Draw("AP");

      graph->GetXaxis()->SetTitle(xtitle.c_str());
      graph->GetXaxis()->SetTitleSize(0.055);
      graph->GetXaxis()->SetTitleOffset(0.8);
      
      graph->GetYaxis()->SetTitle(ytitle.c_str());
      graph->GetYaxis()->SetTitleSize(0.055);
      graph->GetYaxis()->SetTitleOffset(0.9);
      

      TF1 *fit = new TF1("line_fit", Linear, 100, 1000,1);
      fit ->SetParameters(1.,1);
      fit ->SetRange(100.,1000.);

      graph->Fit(fit, "EQ");   
      fit->SetLineColor(kRed);
      fit->Draw("same");
      
      double a = fit->GetParameter(0);  //slope
      
      TString eq;
      if (ImgId == "Ngamma_vs_Einc") eq.Form("N_{\\gamma} = (%.2f \\pm %.1e)E_{inc}", a, fit->GetParError(0));
      else eq.Form("E_{Measured} = (%.2f \\pm %.1e)E_{inc}", a, fit->GetParError(0));
      

      TLatex latex;
      latex.SetNDC();         // Use normalized coordinates
      latex.SetTextSize(0.05);
      latex.DrawLatex(0.25, 0.85, eq);

      canvas->Update();          
    }

    if (FileName.Contains("Response_vs_Einc")){
      ImgId = "Response_vs_Einc";
      ytitle = "\\frac{E_{Measured}}{E_{incident}}";

      infile.open(FileName);
      while (infile >> x >> y){
        x_vals.push_back(x);
        y_vals.push_back(y);}
      
      if (x_vals.empty()) {
        cerr << "No data found!" << endl;
        return;}
      
      graph = new TGraphErrors(x_vals.size(), &x_vals[0], &y_vals[0], &x_err_vals[0], &y_err_vals[0]);

      graph->SetTitle(" ");
      graph->SetMarkerStyle(20);
      graph->SetMarkerSize(1.0);
      graph->SetMarkerColor(kBlue);
      graph->Draw("AP");
          
      graph->GetXaxis()->SetTitle(xtitle.c_str());
      graph->GetXaxis()->SetTitleSize(0.055);
      graph->GetXaxis()->SetTitleOffset(0.8);
      
      graph->GetYaxis()->SetTitle(ytitle.c_str());
      graph->GetYaxis()->SetTitleSize(0.055);
      graph->GetYaxis()->SetTitleOffset(0.9);

      graph->GetYaxis()->SetRangeUser(0.5,1.2);

      TLine *line = new TLine(10, 1, 1090, 1);
      line->SetLineColor(kRed);
      line->SetLineWidth(2);
      line->SetLineStyle(2);
      line->Draw("same");
      
      canvas->Update();
    }

    if (FileName.Contains("Resolution_vs_Einc")){
      ImgId = "Resolution_vs_Einc";
      ytitle = "\\frac{\\sigma_{N_{\\gamma}}}{\\mu_{N_{\\gamma}}}";

      infile.open(FileName);
      while (infile >> x >> y >> y_err){
        x_vals.push_back(x);
        y_vals.push_back(y);
    	x_err_vals.push_back(0);
    	y_err_vals.push_back(y_err);
      }

      if (x_vals.empty()) {
        cerr << "No data found!" << endl;
        return;
      }

      auto EGraph = new TGraphErrors(x_vals.size(), &x_vals[0], &y_vals[0], &x_err_vals[0], &y_err_vals[0]);
      EGraph->SetTitle(" ");
      EGraph->SetName(WidthId.c_str());
      EGraph->SetMarkerStyle(20);
      EGraph->SetMarkerSize(1.0);
      EGraph->SetMarkerColor(kBlue);
      EGraph->Draw("AP");
      
      EGraph->GetXaxis()->SetTitle(xtitle.c_str());
      EGraph->GetXaxis()->SetTitleSize(0.055);
      EGraph->GetXaxis()->SetTitleOffset(0.8);

      EGraph->GetYaxis()->SetTitle(ytitle.c_str());
      EGraph->GetYaxis()->SetTitleSize(0.055);
      EGraph->GetYaxis()->SetTitleOffset(1.0);

    
      TF1 *fit_func = new TF1(WidthId.c_str(), sqrtE, 100, 1000,3);
      fit_func->SetParameters(1,1/2,0);    
      fit_func->SetRange(92, 1020);
      
      EGraph->Fit(fit_func, "EQR0");     //0 is to store the fit function but not draw  
      
      fit_func->SetLineColor(kRed);
      fit_func->SetTitle("Resolution; E_{inc}; Resolution (#sigma)");      
      fit_func->Draw("same");
      
      Fit_array.push_back(fit_func);
      Graph_array_Resol.push_back(EGraph);

      double a = fit_func->GetParameter(0);
      double b = fit_func->GetParameter(1);
      double c = fit_func->GetParameter(2);
        
      TString eq;
      eq.Form("R = %.2fE^{-\\frac{1}{2}} \\oplus %.2fE^{-1} \\oplus %.2f", a, b, c);

      char* chisqrNDF = new char[50];
      char *a_err = new char[50];
      char *b_err = new char[50];
      char *c_err = new char[50];

      sprintf(chisqrNDF, "\\chi^{2}/NDF: %.2e", fit_func->GetChisquare()/fit_func->GetNDF());
      sprintf(a_err, "\\Delta a: %.2e", fit_func->GetParError(0));
      sprintf(b_err, "\\Delta b: %.2e", fit_func->GetParError(1));
      sprintf(c_err, "\\Delta c: %.2e", fit_func->GetParError(2));


      TLatex latex;
      latex.SetNDC();
      latex.SetTextSize(0.045);
      latex.DrawLatex(0.45, 0.85, eq);
      latex.DrawLatexNDC(0.45,0.78, chisqrNDF);
      latex.DrawLatexNDC(0.45,0.71, a_err);
      latex.DrawLatexNDC(0.45,0.64, b_err);
      latex.DrawLatexNDC(0.45,0.57, c_err);

      canvas->Update();
    }

    if (FileName.Contains("Efficiency_vs_Einc")){
      ImgId = "Efficiency_vs_Einc";
      ytitle = "Efficiency";

      infile.open(FileName);
      while (infile >> x >> y){
        x_vals.push_back(x);
        y_vals.push_back(y);}
      
      if (x_vals.empty()) {
        cerr << "No data found!" << endl;
        return;}
      
      auto graph = new TGraphErrors(x_vals.size(), &x_vals[0], &y_vals[0], &x_err_vals[0], &y_err_vals[0]);

      graph->SetTitle(" ");
      graph->SetMarkerStyle(20);
      graph->SetMarkerSize(1.0);
      graph->SetMarkerColor(kBlue);
      graph->SetName(WidthId.c_str());
      graph->Draw("AP");
          
      graph->GetXaxis()->SetTitle(xtitle.c_str());
      graph->GetXaxis()->SetTitleSize(0.055);
      graph->GetXaxis()->SetTitleOffset(0.5);
      
      graph->GetYaxis()->SetTitle(ytitle.c_str());
      graph->GetYaxis()->SetTitleSize(0.055);
      graph->GetYaxis()->SetTitleOffset(0.5);

      graph->GetYaxis()->SetRangeUser(0.0,1.1);

      Graph_array_Eff.push_back(graph);
      
      canvas->Update();
    }


    infile.close();

    //Defining the folder names:
    string MainFolder, folder, ImgName;
    MainFolder = PlotsFolder + "Scintillation/Fits/Resol_Resp";
    string Width, ImageId;

    if (FileName.Contains("1cm")) Width = "1cm";
    else if (FileName.Contains("2.5cm")) Width = "2.5cm";
    else if (FileName.Contains("_5cm")) Width = "5cm";
    else if (FileName.Contains("10cm")) Width = "10cm";

    folder = MainFolder + "/" + MatName + "/" + Width + "/";
    
    ImgName = MatName + "_" + Width + "_" + ImgId;      
    
    //saving the canvas
    string image_name = folder + ImgName + ".png";
    cout << image_name << endl;
    canvas->SaveAs(image_name.c_str());

    // image_name = folder + ImgName + ".pdf";
    // canvas->SaveAs(image_name.c_str());
  }
}


void FitFunc(string pathname)
{
  char *hist_name = new char[100];
  char *full_path = new char[1000];
  char *leg_head = new char[100];
  string filetag, folder, WidthId;

  double Mu_calib = 0, E_calib = 0;

  // vector<TString> Folders = {"PbWO4_out_root_files", "BGO_out_root_files", "LYSO_out_root_files", "GAGG_out_root_files", "LaBr3_out_root_files", "Plastic_out_root_files"};

  vector<TString> Folders = {"LYSO_out_root_files"};

  for (int iFolder = 0; iFolder < Folders.size(); iFolder++){
    Graph_array_Resol.clear();
    Graph_array_Eff.clear();
    Fit_array.clear();
    vector<string> f;
   
    //string FileFolder = "BGO_out_root_files";
    //string FileFolder = "PbWO4_out_root_files";
    TString FileFolder = Folders[iFolder];
    //string FileFolder = "LYSO_out_root_files";
    //string FileFolder = "LaBr3_out_root_files";
    //string FileFolder = "out_root_files";
    TSystemDirectory dir(FileFolder.Data(), FileFolder.Data());
    TList* files = dir.GetListOfFiles();
    if (!files) return;

    TIter next(files);
    TObject* obj;

    while ((obj = next())) {
      TString fname = obj->GetName();
      if (!obj->IsA()->InheritsFrom("TSystemFile")) continue;
      if (fname.EndsWith(".root")) {
	f.push_back(string(FileFolder.Data()) + "/" + string(fname.Data()));
      }
    }

    cout << "\n\nSize of the files vector is: " << f.size() << endl;

    string name = "nOptical_Photons";
    sprintf(hist_name,"%s",name.c_str());  

    string PlotFolder = "test_plots/";
    string txtPath = PlotFolder + "Scintillation/Fits/txtFiles/";
    string MainFolder = PlotFolder + "Scintillation/Fits/ItrFit";;

    string Energy, Width, txtFileId, Material;
    int MatIdx, rebin;
  
    vector<string> Energies = {"511keV", "100keV", "150keV", "300keV", "450keV", "600keV", "800keV", "1000keV"};
    vector<string> Materials = {"PbWO4", "BGO", "LYSO", "GAGG", "LaBr3"};
    vector<string> Widths = {"1cm", "2.5cm", "5cm", "10cm"};
    vector<int> Rebin_Vals = {1, 2, 5, 5, 10};

    for(int i=0; i< Materials.size(); i++)
      {
	if(FileFolder.Contains(Materials[i])) {Material = Materials[i]; MatIdx = i;}	      	     
      }

    rebin = Rebin_Vals[MatIdx];
    cout << "\n==========rebin of " << Material << ": " << rebin << "\n\n";
  
    for(int iFile=0; iFile < f.size(); iFile++){   
      TFile *root_file = new TFile(f[iFile].c_str());
      vector<TH1F> HistList;
      TH1F* hist = (TH1F*)root_file->Get(hist_name);
      int xmax, xmin;
      vector<TH1*> hist_list;

      hist_list.push_back(hist);

      int xrange=0.0;
      string xtitle, ytitle;

      xtitle = hist->GetXaxis()->GetTitle();
      ytitle = hist->GetYaxis()->GetTitle();

    
      TString FileId = f[iFile];

      for(int i=0; i< Energies.size(); i++)
	{
	  if(FileId.Contains(Energies[i])) Energy = Energies[i];		  
	}

      for(int i=0; i< Widths.size(); i++)
	{
	  if(FileId.Contains("_" + Widths[i])) Width = Widths[i];	      
	}

      
      filetag = Material + "_" + Energy + "_" + Width;
      //folder = "Trans_20cmX20cm/" + MainFolder + "/" + Material + "/depth_" + Width;
      folder = MainFolder + "/" + Material + "/" + Width;
      txtFileId = Material + "_" + Width;


    

      cout << "\n==================Details for energy: " << Energy << endl;
      //defininig the initial range for fitting
      if(Material == "BGO"){
	if(Width == "5cm") {
	  if (Energy == "100keV" || Energy == "150keV") {xmin = 600; xmax = 1500;}
	  if (Energy == "300keV") {xmin = 1500; xmax = 3000;}
	  if (Energy == "450keV" ) {xmin = 3000; xmax = 4500;}
	  if (Energy == "511keV" ) {xmin = 3600; xmax = 4500;}
	  if (Energy == "600keV") {xmin = 4000; xmax = 6500;}
	  if (Energy == "800keV") {xmin = 6000 ; xmax = 8000;}
	  if (Energy == "1000keV") {xmin = 7200; xmax=9000;}
	}

	if(Width == "10cm") {
	  if (Energy == "100keV" || Energy == "150keV") {xmin = 600; xmax = 1500;}
	  if (Energy == "300keV") {xmin = 1500; xmax = 3000;}
	  if (Energy == "450keV" || Energy == "511keV") {xmin = 3000; xmax = 4500;}
	  if (Energy == "600keV") {xmin = 4000; xmax = 6500;}
	  if (Energy == "800keV") {xmin = 6000 ; xmax = 8000;}
	  if (Energy == "1000keV") {xmin = 7200; xmax=9000;}
	}

	if(Width == "2.5cm") {
	  if (Energy == "100keV") {xmin = 600; xmax = 900;}
	  if (Energy == "150keV") {xmin = 1000; xmax = 1300;}		
	  if (Energy == "300keV") {xmin = 2000; xmax = 3000;}
	  if (Energy == "450keV") {xmin = 3300; xmax = 4500;}
	  if (Energy == "511keV") {xmin = 3700; xmax = 5000;}
	  if (Energy == "600keV") {xmin = 4500; xmax = 5100;}
	  if (Energy == "800keV") {xmin = 6000 ; xmax = 7000;}
	  if (Energy == "1000keV") {xmin = 7500; xmax=8500;}
	}           	
      }

      if(Material == "PbWO4"){      
	if (Energy != "800keV"){xmax = 250; xmin = 1;}
	else {xmax = 200; xmin = 100;}    
      }

    
      if(Material == "GAGG"){
	if(Width == "5cm") {
	  if (Energy == "100keV") {xmin = 2400; xmax = 3200;}
	  if (Energy == "150keV") {xmin = 3500; xmax = 5000;}
	  if (Energy == "300keV") {xmin = 7500; xmax = 9500;}
	  if (Energy == "450keV") {xmin = 11000; xmax = 14000;}
	  if (Energy == "511keV") {xmin = 13500; xmax = 15500;}
	  if (Energy == "600keV") {xmin = 15500; xmax = 18000;}
	  if (Energy == "800keV") {xmin = 21000 ; xmax = 24000;}
	  if (Energy == "1000keV") {xmin = 26000; xmax = 30000;}
	}

	if(Width == "10cm") {
	  if (Energy == "100keV") {xmin = 2500; xmax = 3500;}
	  if (Energy == "150keV") {xmin = 3500; xmax = 5000;}
	  if (Energy == "300keV") {xmin = 7500; xmax = 9500;}
	  if (Energy == "450keV") {xmin = 12000; xmax = 13500;}
	  if (Energy == "511keV") {xmin = 13800; xmax = 15000;}
	  if (Energy == "600keV") {xmin = 16000; xmax = 18000;}
	  if (Energy == "800keV") {xmin = 21000 ; xmax = 24000;}
	  if (Energy == "1000keV") {xmin = 26000; xmax = 30000;}
	}

	if(Width == "2.5cm") {
	  if (Energy == "100keV") {xmin = 2000; xmax = 3500;}
	  if (Energy == "150keV") {xmin = 3500; xmax = 5000;}
	  if (Energy == "300keV") {xmin = 7500; xmax = 9500;}
	  if (Energy == "450keV") {xmin = 12000; xmax = 14000;}
	  if (Energy == "511keV") {xmin = 13500; xmax = 15500;}
	  if (Energy == "600keV") {xmin = 15000; xmax = 18000;}
	  if (Energy == "800keV") {xmin = 21000 ; xmax = 24000;}
	  if (Energy == "1000keV") {xmin = 26000; xmax = 30000;}
	}

	if(Width == "1cm") {
	  if (Energy == "100keV") {xmin = 2000; xmax = 3500;}
	  if (Energy == "150keV") {xmin = 3500; xmax = 5000;}
	  if (Energy == "300keV") {xmin = 7500; xmax = 9500;}
	  if (Energy == "450keV") {xmin = 12000; xmax = 13500;}
	  if (Energy == "511keV") {xmin = 13500; xmax = 15500;}
	  if (Energy == "600keV") {xmin = 16000; xmax = 18000;}
	  if (Energy == "800keV") {xmin = 21000 ; xmax = 24000;}
	  if (Energy == "1000keV") {xmin = 26000; xmax = 30000;}
	}           	
      }

      //Fitting ranges for LaBr3
      if(Material == "LaBr3"){      
	if(Width == "5cm") {
	  if (Energy == "100keV") {xmin = 5000; xmax = 7000;}
	  if (Energy == "150keV") {xmin = 8000; xmax = 10000;}
	  if (Energy == "300keV") {xmin = 16000; xmax = 20000;}
	  if (Energy == "450keV") {xmin = 26000; xmax = 30000;}
	  if (Energy == "511keV") {xmin = 29000; xmax = 34000;}
	  if (Energy == "600keV") {xmin = 32000; xmax = 40000;}
	  if (Energy == "800keV") {xmin = 46000 ; xmax = 52000;}
	  if (Energy == "1000keV") {xmin = 58000; xmax = 65000;}
	}

	if(Width == "10cm") {
	  if (Energy == "100keV") {xmin = 5000; xmax = 7000;}
	  if (Energy == "150keV") {xmin = 8000; xmax = 10000;}
	  if (Energy == "300keV") {xmin = 16000; xmax = 20000;}
	  if (Energy == "450keV") {xmin = 26000; xmax = 30000;}
	  if (Energy == "511keV") {xmin = 28000; xmax = 34000;}
	  if (Energy == "600keV") {xmin = 32000; xmax = 40000;}
	  if (Energy == "800keV") {xmin = 46000 ; xmax = 52000;}
	  if (Energy == "1000keV") {xmin = 58000; xmax = 65000;}
	}

	if(Width == "2.5cm") {
	  if (Energy == "100keV") {xmin = 5000; xmax = 7000;}
	  if (Energy == "150keV") {xmin = 8000; xmax = 10000;}
	  if (Energy == "300keV") {xmin = 16000; xmax = 20000;}
	  if (Energy == "450keV") {xmin = 26000; xmax = 30000;}
	  if (Energy == "511keV") {xmin = 29000; xmax = 34000;}
	  if (Energy == "600keV") {xmin = 34000; xmax = 38000;}
	  if (Energy == "800keV") {xmin = 46000 ; xmax = 52000;}
	  if (Energy == "1000keV") {xmin = 58000; xmax = 65000;}
	}

	if(Width == "1cm") {
	  if (Energy == "100keV") {xmin = 5000; xmax = 7000;}
	  if (Energy == "150keV") {xmin = 8500; xmax = 10000;}
	  if (Energy == "300keV") {xmin = 16000; xmax = 20000;}
	  if (Energy == "450keV") {xmin = 26000; xmax = 30000;}
	  if (Energy == "511keV") {xmin = 29000; xmax = 34000;}
	  if (Energy == "600keV") {xmin = 34000; xmax = 38000;}
	  if (Energy == "800keV") {xmin = 46000 ; xmax = 52000;}
	  if (Energy == "1000keV") {xmin = 58000; xmax = 65000;}
	}           	
      }

      if(Material == "LYSO"){
	//defininig the initial range for fitting
	if(Width == "5cm") {
	  if (Energy == "100keV") {xmin = 2600; xmax = 3500;}
	  if (Energy == "150keV") {xmin = 3800; xmax = 5000;}
	  if (Energy == "300keV") {xmin = 8200; xmax = 9600;}
	  if (Energy == "450keV") {xmin = 12500; xmax = 14500;}
	  if (Energy == "511keV") {xmin = 14600; xmax = 16000;}
	  if (Energy == "600keV") {xmin = 17000; xmax = 19000;}
	  if (Energy == "800keV") {xmin = 23200 ; xmax = 25000;}
	  if (Energy == "1000keV") {xmin = 28500; xmax = 31000;}
	}

	if(Width == "10cm") {
	  if (Energy == "100keV") {xmin = 2500; xmax = 3500;}
	  if (Energy == "150keV") {xmin = 3500; xmax = 6000;}
	  if (Energy == "300keV") {xmin = 7500; xmax = 10500;}
	  if (Energy == "450keV") {xmin = 12500; xmax = 14400;}
	  if (Energy == "511keV") {xmin = 14000; xmax = 17500;}
	  if (Energy == "600keV") {xmin = 17000; xmax = 20000;}
	  if (Energy == "800keV") {xmin = 22500 ; xmax = 26000;}
	  if (Energy == "1000keV") {xmin = 28000; xmax = 32000;}
	}

	if(Width == "2.5cm") {
	  if (Energy == "100keV") {xmin = 2200; xmax = 3500;}
	  if (Energy == "150keV") {xmin = 3800; xmax = 5000;}
	  if (Energy == "300keV") {xmin = 8000; xmax = 9600;}
	  if (Energy == "450keV") {xmin = 12600; xmax = 14500;}
	  if (Energy == "511keV") {xmin = 14500; xmax = 16000;}
	  if (Energy == "600keV") {xmin = 17000; xmax = 19000;}
	  if (Energy == "800keV") {xmin = 23000 ; xmax = 25000;}
	  if (Energy == "1000keV") {xmin = 28500; xmax = 31100;}
	}                 	
      }


    
      sprintf(full_path,"%s%s/%s_%s",pathname.c_str(),folder.c_str(), "nOpticalPhotns_ItrFit", Energy.c_str());
    
      vector<double> params;
      params = generate_1Dplot(hist_list,full_path,xmax,xmin,leg_head,false,false,false,true,filetag.c_str(),xtitle.c_str(),ytitle.c_str(),rebin);
   
      //storing the parameters in the txt file
      ofstream outFile1(txtPath + "Resolution_vs_Einc_" + txtFileId + ".txt", ios::app);
      ofstream outFile2(txtPath + "Ngamma_vs_Einc_" + txtFileId + ".txt", ios::app);
      ofstream outFile5(txtPath + "Efficiency_vs_Einc_" + txtFileId + ".txt", ios::app);
    
      if (!outFile1 || !outFile2) {
	cerr << "Error opening file!" << endl;
	return 1;
      }

      double Mu, Mu_err, Sigma, Sigma_err, Resol, Resol_err, eff;
      string energy;

      int energy_val;

      Mu = params[0];
      Mu_err = params[1];
      Sigma = abs(params[2]);
      Sigma_err = params[3];
      eff = params[4];

      Resol = Sigma/Mu;
      Resol_err = Resol * sqrt((pow(Sigma_err/Sigma,2)) + (pow(Mu_err/Mu,2)));

      energy = Energy.erase(Energy.length() - 3);
      energy_val = stoi(energy);
  
      outFile1 << energy << fixed << setprecision(3)
	       << setw(12) << Resol << fixed << setprecision(5)
	       << setw(12) << Resol_err << endl;

      outFile2 << energy << fixed << setprecision(2)
	       << setw(12) << Mu << endl;

      outFile5 << energy << fixed << setprecision(4)
	       << setw(12) << eff << endl;

      outFile1.close();
      outFile2.close();
    }

    vector<string> TxtArr;
    double x, y, EMeasured;
    TxtArr = {"BGO_2.5cm", "BGO_5cm", "BGO_10cm", "PbWO4_2.5cm", "PbWO4_5cm", "PbWO4_10cm", "Plastic_5cm", "Plastic_10cm", "Plastic_20cm", "GAGG_1cm", "GAGG_2.5cm", "GAGG_5cm", "GAGG_10cm", "LaBr3_1cm", "LaBr3_2.5cm", "LaBr3_5cm", "LaBr3_10cm", "LYSO_2.5cm", "LYSO_5cm", "LYSO_10cm"};

    vector<string> data_files1;

    //storing the names of all txt files in the array:
    //string txtPath = "test_plots/Scintillation/Fits/txtFiles/";
    TSystemDirectory dir1(txtPath.c_str(), txtPath.c_str());
    TList* files1 = dir1.GetListOfFiles();
    if (!files1) return;

    TIter next1(files1);
    TObject* obj1;

    while ((obj1 = next1())) {
      TString fname = obj1->GetName();
      if (!obj1->IsA()->InheritsFrom("TSystemFile")) continue;
      if (fname.Contains("Ngamma") && fname.EndsWith(".txt")) {
	data_files1.push_back(txtPath + string(fname.Data()));
      }
    }

    for (int itxt = 0; itxt < TxtArr.size(); itxt++){
      for (int ifile = 0; ifile< data_files1.size(); ifile++){
	TString FileName;
	FileName = data_files1[ifile];
	if(FileName.Contains(TxtArr[itxt])){
	  ifstream ChkFile(txtPath + "Ngamma_vs_Einc_" + TxtArr[itxt] + ".txt");
	  ofstream outFile3(txtPath + "Emeasured_vs_Einc_" + TxtArr[itxt] + ".txt");

	  ofstream outFile4(txtPath + "Response_vs_Einc_" + TxtArr[itxt] + ".txt");
	  while (ChkFile >> x >> y){ 
	    if (x == 511) {E_calib = x; Mu_calib = y;}}

	  ChkFile.clear();                // Clear EOF flag
	  ChkFile.seekg(0, ios::beg);     // Go back to beginning
	
	  //now filling the values of the txt file for EMeaured and Response
	  while (ChkFile >> x >> y){
	    EMeasured = (E_calib/Mu_calib)*y;

	    outFile3 << x << setw(12) << EMeasured << endl;

	    outFile4 << x << fixed << setprecision(2)
		     << setw(12) << EMeasured/x << endl;}

	  ChkFile.close();
	  outFile3.close();
	  outFile4.close();			
	}
      }
    }
  

    // fitting the datapoints
    data_plot(Material, txtPath, PlotFolder);

    FitOverlay(Fit_array, Graph_array_Resol, Material, "Resol");
    FitOverlay(Fit_array, Graph_array_Eff, Material, "Eff");

    string command = "./sortTxt.sh " + txtPath + " " + Material;
    system(command.c_str());
  }
}














// struct MixedData {
//      vector<string> str1;
//      string str4;
//      string str2;
//      int intData;
//      double double1;
//      double double2;
//      double double3;
//      double double4;
//      vector<string> str3;
//  };
