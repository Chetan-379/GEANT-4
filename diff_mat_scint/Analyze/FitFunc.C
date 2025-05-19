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

TH1F* setMyRange(TH1F *h1,double xLow,double xHigh){
  //call it after setting last bin as overflow
  double err=0;
  if(xHigh > 13000) return h1;
  if(xLow < -13000) return h1;
  int nMax=h1->FindBin(xHigh);
  h1->SetBinContent(nMax,h1->IntegralAndError(nMax,h1->GetNbinsX(),err));
  h1->SetBinError(nMax,err);
  for(int i=nMax+1;i<=h1->GetNbinsX()+1;i++){
    h1->SetBinContent(i,0);
    h1->SetBinError(i,0);

  }
  h1->GetXaxis()->SetRangeUser(xLow,xHigh);
  cout<<xLow<<"\t"<<xHigh<<"\t"<<"set range"<<endl;
  return h1;
}

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

  TLegend *legend = new TLegend(0.71, 0.84, 0.85, 0.90);
  legend->SetTextSize(0.035);
  legend->SetLineColor(kWhite);
  legend->SetNColumns(2);
  legend->SetHeader(title);
  double MaxY=0;
  for(int i = 0; i < (int)hist.size(); i++) {
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
      //hist.at(i)->Scale(1.0 / hist.at(i)->Integral());
      hist.at(i)->GetYaxis()->SetTitle("Normalized");
    } else {
      hist.at(i)->GetYaxis()->SetTitle("Entries");
    }

    hist.at(i)->GetXaxis()->SetTitle("nOpticalPhotons");

    int maxBin = -1;
    double maxContent = -1;
    for (int j = 1; j <= hist.at(i)->GetNbinsX(); ++j) {
      double center = hist.at(i)->GetBinCenter(j);
      double content = hist.at(i)->GetBinContent(j);

      // Skip the peak at 0
      if (center >= -0.5 && center <= +0.5) continue;

      if (content > maxContent) {
	maxContent = content;
	maxBin = j;
      }
    }

    // Iterative fitting
    //int peakBin = hist.at(i)->GetMaximumBin();
    int peakBin = maxBin;
    double peakValue = hist.at(i)->GetXaxis()->GetBinCenter(peakBin);
    double totalIntegral = hist.at(i)->Integral();
    double targetIntegral = 0.95 * totalIntegral;

    int binLeft = peakBin, binRight = peakBin;
    double cumulativeSum = hist.at(i)->GetBinContent(peakBin);

    while (cumulativeSum < targetIntegral) {
      if (binLeft > 1) cumulativeSum += hist.at(i)->GetBinContent(--binLeft);
      if (binRight < hist.at(i)->GetNbinsX()) cumulativeSum += hist.at(i)->GetBinContent(++binRight);
      if (binLeft == 1 && binRight == hist.at(i)->GetNbinsX()) break;
    }

    double fitMin = hist.at(i)->GetXaxis()->GetBinCenter(binLeft);
    double fitMax = hist.at(i)->GetXaxis()->GetBinCenter(binRight);

    //        fitMin=xmin;
    fitMin=xmin;
    fitMax=xmax;

    double old_mean=hist.at(i)->GetMean();
    double old_sigma = hist.at(i)->GetRMS();
    // cout << "\nsigma is: " << old_sigma << endl;

    for (int itr = 0; itr < 400; ++itr) {
      //hist.at(i)->Fit(Fit_func[i], "RE");
      Fit_func[i] = new TF1("GaussFit", Gauss, fitMin, fitMax, 3);
      Fit_func[i]->SetParameters(peakValue, old_mean, 5.55);
      //else Fit_func[i]->SetParameters(peakValue, old_mean, old_sigma);
      hist.at(i)->Fit(Fit_func[i], "RQ");

      double mu = Fit_func[i]->GetParameter(1);
      double sigma = Fit_func[i]->GetParameter(2);
      double chi2byNDF = Fit_func[i]->GetChisquare()/Fit_func[i]->GetNDF();

      double rangeMin = peakValue - 2.0 * abs(sigma);
      double rangeMax = peakValue + 2.0 * abs(sigma);

      if (fabs(fitMin-rangeMin)/rangeMin <= 0.0002) break;

      fitMin = rangeMin;
      fitMax = rangeMax;

      //Fit_func[i]->SetRange(rangeMin, rangeMax);
      Fit_func[i]->SetRange(fitMin, fitMax);

      old_mean = mu;
      old_sigma = sigma;
    }

    //hist.at(i)->GetXaxis()->SetRangeUser(fitMin, fitMax);
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


    // sprintf(param1, "A: %0.2e #pm %0.4f", Fit_func[i]->GetParameter(0), Fit_func[i]->GetParError(0));
    // sprintf(param2, "\\mu: %0.3f \\pm %0.4f", Fit_func[i]->GetParameter(1), Fit_func[i]->GetParError(1));
    // sprintf(param3, "\\sigma: %0.4f \\pm %0.4f", Fit_func[i]->GetParameter(2), Fit_func[i]->GetParError(2));
    // sprintf(chisqrNDF, "\\chi^{2}/NDF: %.4f", Fit_func[i]->GetChisquare()/Fit_func[i]->GetNDF());

    sprintf(param1, "A: %0.2e #pm %0.4f", A, Fit_func[i]->GetParError(0));
    sprintf(param2, "\\mu: %0.3f \\pm %0.4f", Mu, Fit_func[i]->GetParError(1));
    sprintf(param3, "\\sigma: %0.4f \\pm %0.4f", Sigma, Fit_func[i]->GetParError(2));
    sprintf(chisqrNDF, "\\chi^{2}/NDF: %.4f", Fit_func[i]->GetChisquare()/Fit_func[i]->GetNDF());

    latex->DrawLatexNDC(0.70,0.80, param1);
    latex->DrawLatexNDC(0.70,0.75, param2);
    latex->DrawLatexNDC(0.70,0.70, param3);
    latex->DrawLatexNDC(0.70,0.65, chisqrNDF);

    results.push_back(Mu);
    results.push_back(Fit_func[i]->GetParError(1));
    results.push_back(Sigma);
    results.push_back(Fit_func[i]->GetParError(2));

  }

  legend->Draw();
  gPad->Modified();
  gPad->Update();

  if (save_canvas) {
    canvas_n1->SaveAs(Form("%s.png", tag_name));
  }

  return results;
}

// Defining the fit function
double sqrtE(double *x, double *par) {

  double a = par[0];   // coefficient
  double b = par[1];  // any constant
  double c = par[2];

  //return (sqrt((pow(a/sqrt(x[0]), 2)) + (pow(b,2))));
  return (sqrt((pow(a/sqrt(x[0]), 2)) + (pow(b/x[0],2)) + (pow(c,2))));

  // double a = par[0];   // coefficient
  // double n = par[1];  //power
  // double b = par[2];  // any constant
  // return (sqrt((pow(a/pow(x[0],n), 2)) + (pow(b,2))));
}


void data_plot() {
  vector<string> data_files;
  vector<string> xtitle, ytitle;

  //filling all the txt files in a loop

  // TSystemDirectory dir("plots/Scintillation/Fits/txtFiles", "plots/Scintillation/Fits/txtFiles");
  // TList* files = dir.GetListOfFiles();
  // if (!files) return;

  // TIter next(files);
  // TObject* obj;

  // while ((obj = next())) {
  //   TString fname = obj->GetName();
  //   if (!obj->IsA()->InheritsFrom("TSystemFile")) continue;
  //   if (fname.EndsWith(".txt")) {
  //     //f.push_back("out_root_files/" + string(fname.Data()));
  //     data_files.push_back("plots/Scintillation/Fits/txtFiles/" + string(fname.Data()));
  //   }
  // }

  // cout << "\n\nSize of the files vector is: " << data_files.size() << endl;

  for (int i=0; i<data_files.size(); i++){
    TString FileName = data_files[i];
    string ImgId;
    
    ifstream infile;

    vector<double> x_vals, y_vals, x_err_vals, y_err_vals;
    double x, y, y_err;


    // Create canvas
    data_files[i].resize(data_files[i].length() - 4);        //removing .txt part
    TCanvas *canvas = new TCanvas(data_files[i].c_str(), data_files[i].c_str(), 800, 600);
    canvas->SetLeftMargin(0.135);
    canvas->SetRightMargin(0.035);
    canvas->SetTopMargin(0.057);
    canvas->SetBottomMargin(0.11);


    // Fit with a straight line (pol1) and draw
    ;
    TGraph *graph;
    if (FileName.Contains("Emeasured_vs_Einc") || FileName.Contains("Ngamma_vs_Einc")){

      infile.open(FileName);
      while (infile >> x >> y){
        x_vals.push_back(x);
        y_vals.push_back(y);
      }
      cout << "EMEs/Ngama" << endl;
      if (x_vals.empty()) {
        std::cerr << "No data found!" << std::endl;
        return;
      }

      graph = new TGraph(x_vals.size(), &x_vals[0], &y_vals[0]);
      
      
      TF1 *fit = new TF1("fit", "pol1", *std::min_element(&x_vals[0], &x_vals[0]+x_vals.size()), *std::max_element(&x_vals[0], &x_vals[0]+x_vals.size()));
      graph->Fit(fit, "Q");  // "Q" = quiet, no print to stdout
      fit->SetLineColor(kRed);
      fit->Draw("same");

      // cout << "R2value: " << fit->GetR
      double a = fit->GetParameter(0);  // intercept
      double b = fit->GetParameter(1);  // slope

      // Create equation string
      TString eq;
      eq.Form("y = %.2fx + %.2f", b, a);

      // Draw equation on canvas using TLatex
      TLatex latex;
      latex.SetNDC();         // Use normalized coordinates
      latex.SetTextSize(0.05);
      latex.DrawLatex(0.25, 0.85, eq);

      canvas->Update();
      if (FileName.Contains("Ngamma_vs_Einc")) ImgId = "Ngamma_vs_Einc";
      if (FileName.Contains("Emeasured_vs_Einc")) ImgId = "Emeasured_vs_Einc";
    
    }

    if (FileName.Contains("Response_vs_Einc")){
      graph->GetYaxis()->SetRangeUser(0.5,1.2);
      TLine *line = new TLine(10, 1, 1090, 1);
      line->SetLineColor(kRed);
      line->SetLineWidth(2);
      line->SetLineStyle(2);
      line->Draw("same");
      canvas->Update();
      ImgId = "Response_vs_Einc";
    }

    if (FileName.Contains("Resolution_vs_Einc")){

      //infile = string(FileName).c_str();
      infile.open(FileName);
      while (infile >> x >> y >> y_err){
        x_vals.push_back(x);
        y_vals.push_back(y);
    	x_err_vals.push_back(0);
    	y_err_vals.push_back(y_err);
      }

      cout << "resolution" << endl;
      if (x_vals.empty()) {
        std::cerr << "No data found!" << std::endl;
        return;
      }

      auto EGraph = new TGraphErrors(x_vals.size(), &x_vals[0], &y_vals[0], &x_err_vals[0], &y_err_vals[0]);
    
      EGraph->SetTitle(" ");
      EGraph->SetMarkerStyle(20);
      EGraph->SetMarkerSize(1.0);
      EGraph->SetMarkerColor(kBlue);
      EGraph->Draw("AP");
      // EGraph->GetXaxis()->SetTitle("E_{True}");
      // EGraph->GetYaxis()->SetTitle("N_{\\gamma}");
    
      //EGraph->GetXaxis()->SetTitle(xtitle[i].c_str());
      EGraph->GetXaxis()->SetTitleSize(0.055);
      EGraph->GetXaxis()->SetTitleOffset(0.8);

      //EGraph->GetYaxis()->SetTitle(ytitle[i].c_str());
      EGraph->GetYaxis()->SetTitleSize(0.055);
      EGraph->GetYaxis()->SetTitleOffset(0.9);
    
    
      TF1 *fit_func = new TF1("fit_fcn", sqrtE, 100, 1000,3);
      fit_func->SetParameters(1,1/2,0);
    
      //TF1 *fit_func = new TF1("fit_fcn", sqrtE, 100, 1000,2);
      //fit_func->SetParameters(1,0);
      fit_func->SetRange(100, 1000);
    
      EGraph->Fit(fit_func, "EQ");
      // Filename.ReplaceAll(".txt", "");  
      fit_func->SetLineColor(kRed);
      fit_func->Draw("same");

      double a = fit_func->GetParameter(0);
      double b = fit_func->GetParameter(1);
      double c = fit_func->GetParameter(2);
    
      // double a = fit_func->GetParameter(0);
      // double n = fit_func->GetParameter(1);
      // double b = fit_func->GetParameter(2);
    
      // Create equation string
      TString eq;
    
      //eq.Form("R = %.2fE^{-\\frac{1}{2}} \\oplus %.2f", a, b);
      eq.Form("R = %.2fE^{-\\frac{1}{2}} \\oplus %.2fE^{-1} \\oplus %.2f", a, b, c);

      //eq.Form("R = %.2fE^{-%.2f} \\oplus %.2f", a, n, b);


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
      ImgId = "Resolution_vs_Einc";
    }

    infile.close();

    //Defining the folder names:

    string MainFolder, folder, ImgName;
    MainFolder = "plots/Scintillation/Fits/Resol_Resp";
    string Width, ImageId;

    if(FileName.Contains("BGO")){
      if (FileName.Contains("2.5cm")) {folder = MainFolder + "/BGO/2.5cm/"; Width = "2.5cm";}
      if (FileName.Contains("_5cm")) {folder = MainFolder + "/BGO/5cm/"; Width = "5cm";}
      if (FileName.Contains("10cm")) {folder = MainFolder + "/BGO/10cm/"; Width = "10cm";}
      
      ImgName = "BGO_" + Width + ImgId;
      
    }
    
    if(FileName.Contains("PbWO4")){
      if (FileName.Contains("2.5cm")) {folder = MainFolder + "/PbWO4/2.5cm/"; Width = "2.5cm";}
      if (FileName.Contains("_5cm")) {folder = MainFolder + "/PbWO4/5cm/"; Width = "5cm";}
      if (FileName.Contains("10cm")) {folder = MainFolder + "/PbWO4/10cm/"; Width = "10cm";}
      
      ImgName = "PbWO4_" + Width + "_" + ImgId;
    }
    
    if(FileName.Contains("Plastic")){
      if (FileName.Contains("20cm")) {folder = MainFolder + "/Plastic/20cm/"; Width = "20cm";}
      if (FileName.Contains("5cm")) {folder = MainFolder + "/Plastic/5cm/"; Width = "5cm";}
      if (FileName.Contains("10cm")) {folder = MainFolder + "/Plastic/10cm/"; Width = "10cm";}
      
      ImgName = "Plastic_" + Width + "_"  + ImgId;
    }
    
    //saving the canvas
    string image_name = folder + ImgName + ".png";
    cout << image_name << endl;
    canvas->SaveAs(image_name.c_str());
  }
}


void FitFunc(string pathname)
{
  char *hist_name = new char[100];
  char *full_path = new char[1000];
  char *leg_head = new char[100];
  vector<string> f;
  string filetag;
  string folder;

  double Mu_calib = 0, E_calib = 0;

  std::vector<std::string> rootFiles;

  //filling all the root files in an array f
  //TSystemDirectory dir("out_root_files", "out_root_files");
  TSystemDirectory dir("check_out_root_files", "check_out_root_files");
  TList* files = dir.GetListOfFiles();
  if (!files) return;

  TIter next(files);
  TObject* obj;

  while ((obj = next())) {
    TString fname = obj->GetName();
    if (!obj->IsA()->InheritsFrom("TSystemFile")) continue;
    if (fname.EndsWith(".root")) {
      //f.push_back("out_root_files/" + string(fname.Data()));
      f.push_back("check_out_root_files/" + string(fname.Data()));
    }
  }

  cout << "\n\nSize of the files vector is: " << f.size() << endl;


  vector<int >rebin = {1,1,1,1,1,1,1,1,1}; //keep it 1 if you don't want to change hist bins

  vector<int>xmax = {2000,2000,16,16,2000,2000,2000,2000};
  vector<int>xmin = {-200,0,0,0,0,0,0,0};
  
  string name = "nOptical_Photons";
  sprintf(hist_name,"%s",name.c_str());
  
  //cleaning the contents of the txt files before starting
  int result = system("./txtClean.sh");  // Ensure script is executable
  
  if (result != 0) {
    cerr << "Shell script failed with code: " << result << endl;
  } else {
    cout << "Old txt files removed!!." << endl;
  }

  //starting looping over the files
  for(int iFile=0; iFile < f.size(); iFile++){
    TFile *root_file = new TFile(f[iFile].c_str());
    vector<TH1F> HistList;
    TH1F* hist = (TH1F*)root_file->Get(hist_name);

    vector<TH1*> hist_list;

    hist_list.push_back(hist);

    int xrange=0.0;
    string xtitle, ytitle;

    xtitle = hist->GetXaxis()->GetTitle();
    ytitle = hist->GetYaxis()->GetTitle();


    //specifying the path of the folder etc
    string MainFolder;
    //if(histId.Contains("Edep") || histId.Contains("Compt")) MainFolder = "Pho_Interact";
    MainFolder = "Scintillation/Fits/ItrFit";
    TString FileId = f[iFile];
    string Energy, Width, txtFileId;
    
    //making identifier for energy
    if (FileId.Contains("511keV")) Energy = "511keV";
    if (FileId.Contains("100keV")) Energy = "100keV";
    if (FileId.Contains("150keV")) Energy = "150keV";
    if (FileId.Contains("300keV")) Energy = "300keV";
    if (FileId.Contains("450keV")) Energy = "450keV";
    if (FileId.Contains("600keV")) Energy = "600keV";
    if (FileId.Contains("800keV")) Energy = "800keV";
    if (FileId.Contains("1000keV")) Energy = "1000keV";


    if(FileId.Contains("BGO")){
      if (FileId.Contains("2.5cm")) {folder = "plots/" + MainFolder + "/BGO/2.5cm"; Width = "2.5cm";}
      if (FileId.Contains("_5cm")) {folder = "plots/" + MainFolder + "/BGO/5cm"; Width = "5cm";}
      if (FileId.Contains("10cm")) {folder = "plots/" + MainFolder + "/BGO/10cm"; Width = "10cm";}

      filetag = "BGO_" + Energy + "_" + Width;
      txtFileId = "BGO_" + Width;
	
    }

    if(FileId.Contains("PbWO4")){
      if (FileId.Contains("2.5cm")) {folder = "plots/" + MainFolder + "/PbWO4/2.5cm"; Width = "2.5cm";}
      if (FileId.Contains("_5cm")) {folder = "plots/" + MainFolder + "/PbWO4/5cm"; Width = "5cm";}
      if (FileId.Contains("10cm")) {folder = "plots/" + MainFolder + "/PbWO4/10cm"; Width = "10cm";}

      filetag = "PbWO4_" + Energy + "_" + Width;
      txtFileId = "PbWO4_" + Width;
    }

    if(FileId.Contains("Plastic")){
      filetag = "Plastic_" + Energy + "_" + Width;
      if (FileId.Contains("20cm")) {folder = "plots/" + MainFolder + "/Plastic/20cm"; Width = "20cm";}
      if (FileId.Contains("5cm")) {folder = "plots/" + MainFolder + "/Plastic/5cm"; Width = "5cm";}
      if (FileId.Contains("10cm")) {folder = "plots/" + MainFolder + "/Plastic/10cm"; Width = "10cm";}

      filetag = "Plastic_" + Energy + "_" + Width;
      txtFileId = "Plastic_" + Width;
    }

    sprintf(full_path,"%s%s/%s_%s",pathname.c_str(),folder.c_str(), "nOpticalPhotns_ItrFit", Energy.c_str());

    vector<double> params;
    params = generate_1Dplot(hist_list,full_path,xmax[0],xmin[0],leg_head,false,false,false,true,filetag.c_str(),xtitle.c_str(),ytitle.c_str(),rebin[0]);

    //storing the parameters in the txt file
    ofstream outFile1("plots/Scintillation/Fits/txtFiles/Resolution_vs_Einc_" + txtFileId + ".txt", ios::app), outFile2("plots/Scintillation/Fits/txtFiles/Ngamma_vs_Einc_" + txtFileId + ".txt", ios::app);
    
    if (!outFile1 || !outFile2) {
      std::cerr << "Error opening file!" << std::endl;
      return 1;
    }

    double Mu, Mu_err, Sigma, Sigma_err, Resol, Resol_err;
    string energy;

    int energy_val;

    Mu = params[0];
    Mu_err = params[1];
    Sigma = params[2];
    Sigma_err = params[3];

    Resol = Sigma/Mu;
    Resol_err = Resol * sqrt((pow(Sigma_err/Sigma,2)) + (pow(Mu_err/Mu,2)));

    energy = Energy.erase(Energy.length() - 3);
    energy_val = stoi(energy);

    //writing the txt files
    outFile1 << energy << fixed << setprecision(2)
	     << setw(12) << Resol << fixed << setprecision(5)
	     << setw(12) << Resol_err << endl;

    outFile2 << energy << fixed << setprecision(2)
	     << setw(12) << Mu << endl;

    outFile1.close();
    outFile2.close();
  }

  vector<string> TxtArr;
  double x, y, EMeasured;
  TxtArr = {"BGO_2.5cm", "BGO_5cm", "BGO_10cm", "PbWO4_2.5cm", "PbWO4_5cm", "PbWO4_10cm", "Plastic_5cm", "Plastic_10cm", "Plastic_20cm"};

  vector<string> data_files1;
  //storing the names of all txt files in the array:
   TSystemDirectory dir1("plots/Scintillation/Fits/txtFiles", "plots/Scintillation/Fits/txtFiles");
  TList* files1 = dir1.GetListOfFiles();
  if (!files1) return;

  TIter next1(files1);
  TObject* obj1;

  while ((obj1 = next1())) {
    TString fname = obj1->GetName();
    if (!obj1->IsA()->InheritsFrom("TSystemFile")) continue;
    if (fname.Contains("Ngamma") && fname.EndsWith(".txt")) {
      data_files1.push_back("plots/Scintillation/Fits/txtFiles/" + string(fname.Data()));
    }
  }

  for (int itxt = 0; itxt < TxtArr.size(); itxt++){
    for (int ifile = 0; ifile< data_files1.size(); ifile++){
      TString FileName;
      FileName = data_files1[ifile];
      if(FileName.Contains(TxtArr[itxt])){
	ifstream ChkFile("plots/Scintillation/Fits/txtFiles/Ngamma_vs_Einc_" + TxtArr[itxt] + ".txt");
	ofstream outFile3("plots/Scintillation/Fits/txtFiles/Emeasured_vs_Einc_" + TxtArr[itxt] + ".txt");

	ofstream outFile4("plots/Scintillation/Fits/txtFiles/Response_vs_Einc_" + TxtArr[itxt] + ".txt");
	while (ChkFile >> x >> y){ 
	  if (x == 511) {E_calib = x; Mu_calib = y;}   	  
	}

	ChkFile.clear();                // Clear EOF flag
	ChkFile.seekg(0, std::ios::beg); // Go back to beginning
	

	//now filling the values of the txt file for EMeaured and Response
	while (ChkFile >> x >> y){
	  EMeasured = (E_calib/Mu_calib)*y;

	  outFile3 << x << setw(12) << EMeasured << endl;

	  outFile4 << x << fixed << setprecision(2)
		   << setw(12) << EMeasured/x << endl;
	}

	ChkFile.close();
	outFile3.close();
	outFile4.close();
			
      }
    }
  }

  // fitting the datapoints
  data_plot();
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


