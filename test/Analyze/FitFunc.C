int line_width[12] = {2,2,2,2,2,2,2,2,2,2,2,2};
int line_style[12] = {1,1,1,1,1,1,1,1,1,1,1,1};                                                                   
//int marker_style[12] = {21,43,22,29,33,34,39,41,43,45,47,23};
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
    // double dx = x[0] - par[0];
    // double sigmaL = par[1];
    // double sigmaR = par[2];
    // double sigma = (sigmaL + sigmaR) / 2;
    // double alphaL = par[3];
    // double alphaR = par[4];

  double A = par[0];   // Amplitude
  double mu = par[1];  // Mean
  double sigma = par[2]; // Std. deviation
  double arg = ((x[0] - mu) * (x[0] - mu)) / (sigma*sigma);

  // return (A/(sqrt(2*3.14)*par[2])) * TMath::Exp(-0.5 * arg * arg);
   return A * TMath::Exp(-0.5 * arg);

    
    // if (dx < 0)
    //     return par[5] * TMath::Exp(-dx * dx / (2 * sigmaL * sigmaL + alphaL * dx * dx));
    // else
    //     return par[5] * TMath::Exp(-dx * dx / (2 * sigmaR * sigmaR + alphaR * dx * dx));
}

const int nfiles=100;                                                                                                                                                             
TFile *f[nfiles];

// void generate_1Dplot(vector<TH1*> hist, char const *tag_name="", char const *xlabel="", char const *ylabel="", int rebin=-1, double ymin=0, double ymax=0, double xmin=-1, double xmax=-1, char const *leg_head="", bool normalize=false, bool log_flag=false, bool DoRebin=false, bool save_canvas=true, char const *title="", vector<string> legend_texts={"nil"})
void generate_1Dplot(vector<TH1*> hist, char const *tag_name="", int xmax=-1,int xmin=-1,char const *leg_head="",
		     bool normalize=false, bool log_flag=true, bool DoRebin=false, bool save_canvas=true, char const *title="", const char* xtitile="", const char* ytitile="", int rebin=-1){  

  //{   //normalize=true;
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
 for(int i = 0; i < (int)hist.size(); i++) {

        if(normalize) {
            //hist.at(i)->Scale(1.0 / hist.at(i)->Integral());
            hist.at(i)->GetYaxis()->SetTitle("Normalized");
        } else {
            hist.at(i)->GetYaxis()->SetTitle("Entries");
        }

	hist.at(i)->GetXaxis()->SetTitle("nOpticalPhotons");

        //hist.at(i)->Rebin(rebin);
        //hist.at(i)->GetXaxis()->SetRangeUser(xmin, xmax);
        //hist.at(i)->Smooth(1);

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
	cout <<"Here2"<<endl;

        /*TF1* fitFunc = new TF1("CruijffFit", Cruijff, fitMin, fitMax, 6);
        fitFunc->SetParameters(peakValue, 0.01, 0.01, 0.1, 0.1, hist.at(i)->GetMaximum());
        fitFunc->SetParLimits(1, 0.00001, 0.1);
        fitFunc->SetParLimits(2, 0.00001, 0.1);*/

	//Fit_func[i] = new TF1("GaussFit", Gauss, fitMin, fitMax, 3);
	//Fit_func[i] = new TF1("GaussFit", Gauss, 0, 250, 3);
        //Fit_func[i]->SetParameters(peakValue, 0.01, 0.01, 0.1, 0.1, hist.at(i)->GetMaximum());
	//Fit_func[i]->SetParameters(peakValue, 90.5, 13.55); 

	// Fit_func[i]->SetParLimits(1, 0.00001, 0.1);
        // Fit_func[i]->SetParLimits(2, 0.00001, 0.1);
	
	double old_mean=hist.at(i)->GetMean();
        double old_sigma = hist.at(i)->GetRMS();
        //for (int iter = 0; iter < 400; ++iter) {
	for (int iter = 0; iter < 400; ++iter) {
	  //hist.at(i)->Fit(Fit_func[i], "RE");
	  Fit_func[i] = new TF1("GaussFit", Gauss, fitMin, fitMax, 3);
	  //Fit_func[i]->SetParameters(peakValue/(sqrt(2*3.14)*old_sigma), old_mean, old_sigma);
	  Fit_func[i]->SetParameters(peakValue, old_mean, old_sigma);
	  //Fit_func[i]->SetParLimits(2, 0.0000, 100000);
	  
	  hist.at(i)->Fit(Fit_func[i], "R");
	  
	  
	  // double mean = Fit_func[i]->GetParameter(0);
	  // double sigmaL = Fit_func[i]->GetParameter(1);
	  // double sigmaR = Fit_func[i]->GetParameter(2);
	  // double sigma = (sigmaL + sigmaR)/2;
	  // if ((abs(old_mean-mean)/mean) < 0.0002 && abs(old_sigma-sigma)/sigma < 0.0002) break;
	  // old_sigma=sigma;
	  // old_mean=mean;
	  // cout <<"Here3"<<"\t"<<iter<<endl;

	  double mu = Fit_func[i]->GetParameter(1);
	  double sigma = Fit_func[i]->GetParameter(2);
	  double chi2byNDF = Fit_func[i]->GetChisquare()/Fit_func[i]->GetNDF();

	  // hist.at(i)->Fit("gaus", "R");

	  // TF1 *fitFunc = hist.at(i)->GetFunction("gaus");
	  // double mu    = fitFunc->GetParameter(1); // mean
	  // double sigma = fitFunc->GetParameter(2); // standard deviation
	  // double chi2byNDF = fitFunc->GetChisquare()/fitFunc->GetNDF();
	  
	  cout << "mu: " << mu << endl;
	  cout << "sigma: " << sigma << endl;
	  cout << "chi2/NDF: " << chi2byNDF << "\n\n";


	  //if (sigma == 0) sigma = 10e-10;
	  if (sigma < 0) sigma = abs(sigma);

	  // double rangeMin = mu - 2.0 * sigma;
	  // double rangeMax = mu + 2.0 * sigma;

	  double rangeMin = peakValue - 2.0 * abs(sigma);
	  double rangeMax = peakValue + 2.0 * abs(sigma);

	  // if (chi2byNDF > 25.)  {
	  //    rangeMin = peakValue - 1.0 * sigma;
	  //    rangeMax = peakValue + 1.0 * sigma;

	  //}

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
        //hist.at(i)->GetXaxis()->SetRangeUser(0.96, 1.04);
        //fitFunc->SetRange(0.96, 1.04);
        // hist.at(i)->GetYaxis()->SetRangeUser(0, 1.5*MaxY);
        // hist.at(i)->SetMarkerColor(line_color[i]);
	// hist.at(i)->SetLineColor(line_color[i]);
        // hist.at(i)->SetMarkerStyle(marker_style[i]);
        // hist.at(i)->SetLineWidth(1);
	// hist.at(i)->SetMarkerSize(1);
        // Fit_func[i]->SetLineColor(line_color[i]);
        //hist.at(i)->Draw(i == 0 ? "" : "SAME");
	hist.at(i)->Draw("hist");
	Fit_func[i]->Draw("SAME");


	TLatex *latex = new TLatex();
	latex->SetTextSize(0.03);
    
	char* param1 = new char[50];
	char* param2 = new char[50];
	char* param3 = new char[50];
	char* chisqrNDF = new char[50];
    
    
	sprintf(param1, "A: %0.2e #pm %0.4f", Fit_func[i]->GetParameter(0), Fit_func[i]->GetParError(0));
	sprintf(param2, "\\mu: %0.3f \\pm %0.4f", Fit_func[i]->GetParameter(1), Fit_func[i]->GetParError(1));
	sprintf(param3, "\\sigma: %0.4f \\pm %0.4f", Fit_func[i]->GetParameter(2), Fit_func[i]->GetParError(2));
	sprintf(chisqrNDF, "\\chi^{2}/NDF: %.4f", Fit_func[i]->GetChisquare()/Fit_func[i]->GetNDF());
    
	latex->DrawLatexNDC(0.70,0.80, param1);
	latex->DrawLatexNDC(0.70,0.75, param2);
	latex->DrawLatexNDC(0.70,0.70, param3);
	latex->DrawLatexNDC(0.70,0.65, chisqrNDF);

	// latex->DrawLatexNDC(0.30,0.80, param1);
	// latex->DrawLatexNDC(0.30,0.75, param2);
	// latex->DrawLatexNDC(0.30,0.70, param3);
	// latex->DrawLatexNDC(0.30,0.65, chisqrNDF);
	
        //fitFunc->Draw("SAME");
        
       /* TLatex latex;
        latex.SetNDC();
        latex.SetTextColor(line_color[i]);
        latex.SetTextSize(0.035);
        latex.DrawLatex(0.8, 0.9, Form("#mu = %.4f", fitFunc->GetParameter(0)));
        latex.DrawLatex(0.8, 0.85, Form("#sigma_{L} = %.4f", fitFunc->GetParameter(1)));
        latex.DrawLatex(0.8, 0.8, Form("#sigma_{R} = %.4f", fitFunc->GetParameter(2)));
        legend->AddEntry(hist.at(i), legend_texts[i].c_str(), "p");*/

	// TLatex latex;
        // latex.SetNDC();
        // latex.SetTextColor(line_color[i]);
        // latex.SetTextSize(0.028);
        // latex.SetText(0.8, 0.9-(i*0.13), Form("#mu = %.4f", Fit_func[i]->GetParameter(0)));
        // latexVec[i].push_back(latex);
        // latex.SetText(0.8, 0.87-(i*0.13), Form("#sigma_{L} = %.4f", Fit_func[i]->GetParameter(1)));
        // latexVec[i].push_back(latex);
        // latex.SetText(0.8, 0.84-(i*0.13), Form("#sigma_{R} = %.4f", Fit_func[i]->GetParameter(2)));
        // latexVec[i].push_back(latex);
        // latex.SetText(0.8, 0.8-(i*0.13), Form("#chi^{2}/NDF = %.4f", Fit_func[i]->GetChisquare()/Fit_func[i]->GetNDF()));
        // latexVec[i].push_back(latex);
        // legend->AddEntry(hist.at(i), legend_texts[i].c_str(), "p");

        // cout<<"Legend entry " <<i<<endl;
    }
 //    for(int i = 0; i < (int)latexVec.size(); i++) {
//     for(int j = 0; j < (int)latexVec[i].size(); j++) {
//         latexVec[i][j].Draw();  // Draw each TLatex object
//     }
//         hist.at(i)->Draw("same PE");
//         Fit_func[i]->Draw("same");
//         cout<<"Chi2"<<"\t"<< Fit_func[i]->GetChisquare() << "\t"<< "NDF : "<<Fit_func[i]->GetNDF()<<endl;
// }

    

    legend->Draw();
    gPad->Modified();
    gPad->Update();


    
    if (save_canvas) {
        canvas_n1->SaveAs(Form("%s.png", tag_name));
        //canvas_n1->SaveAs(Form("%s.pdf", tag_name));
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

void FitFunc(string pathname, string Energy)
{
  char *hist_name = new char[100];
  char *full_path = new char[1000];
  char *leg_head = new char[100];
  vector<string> f;
  vector<string> filetag;
  vector<string> folder;

  // f = {"out_root_files/Si_FluON_10000evts_out.root", "out_root_files/PbWO4_FluON_10000evts_out.root", "out_root_files/BGO_FluON_10000evts_out.root"};
  // filetag={"Si", "PbWO4", "BGO"};
  // folder = {"plots/Si", "plots/PbWO4", "plots/BGO"};
  //f = {"out_root_files/gammaE_1000keV_water_RI_Birk_0_SY200_pram_reproduce_PbWO4_out.root"};
  f={"out_root_files/gammaE_" + Energy + "_water_RI_Birk_0_SY200_pram_reproduce_PbWO4_out.root"};
  filetag={"\\gamma_{inc} = " + Energy};
  folder = {"plots/ScintPlots/ItrFit"};

  vector<int >rebin = {1,1,1,1,1,1,1,1,1}; //keep it 1 if you don't want to change hist bins
      
  //x axis range
  vector<int>xmax = {2000,2000,16,16,2000,2000,2000,2000};
  vector<int>xmin = {-200,0,0,0,0,0,0,0};

  //char* hist_name  = new char[200];

  string name = "nOptical_Photons";
  sprintf(hist_name,"%s",name.c_str());
  
  for(int iFile=0; iFile < f.size(); iFile++){ //looping over each file	
    TFile *root_file = new TFile(f[iFile].c_str());

    // Loop over all objects in the file
    vector<TH1F> HistList;
    
    TH1F* hist = (TH1F*)root_file->Get(hist_name);
    
    vector<TH1*> hist_list;
    
    hist_list.push_back(hist);
    
    int xrange=0.0;	  
    
    //setting up the axis titles
     string xtitle, ytitle;
    // if ( histType == "TH2I" || histType == "TH2F" || histType == "TH2D") {
     	xtitle = hist->GetXaxis()->GetTitle();
     	ytitle = hist->GetYaxis()->GetTitle();//}
    
    //   else {
    // 	xtitle = hist->GetName();
    // 	ytitle = "Entries";}

      //path to save the files a jpg or pdf
      sprintf(full_path,"%s%s/%s_%s",pathname.c_str(),folder[0].c_str(), "nOpticalPhotns_ItrFit", Energy.c_str());

      //calling generate_1Dplot which will take this vector of histograms 
      generate_1Dplot(hist_list,full_path,xmax[0],xmin[0],leg_head,false,false,false,true,filetag[0].c_str(),xtitle.c_str(),ytitle.c_str(),rebin[0]);
    }
}



// void FitFunc(){

//     char* hname = new char[200];
  
//   char* hist_name = new char[200];
  
//   char* title= new char[2000];
 
//   char *leg_head = new char[200];
 
//   int n=0;
//   int n_files=1;
 
//     //f[0] = new TFile("EGM_DRN_plots.root");
//     f[0] = new TFile("plot.root");
    
//     vector<string> filetag=  {""};
//     //vector<vector<string>> varName;
//     //vector<vector<string>> legend_texts;
//     //vector<string> xLabel;
//     vector<string> loghist;
//     vector<string> norm;
// MixedData varName[]= {  // {{Array of names of plots},Title of plot, xlabel, rebin, ymin, ymax , xmin, xmax, {legend}}

// {{"DRN_ErecoByEgen","BDT_ErecoByEgen"},"ECorrByEGen","E_{Pred}/E_{Gen}",100,0,0.000001,0.95,1.05,{"DRN","BDT"}},
// //{"DRN_ErecoByEgen","ECorrByEGen","E_{Pred}/E_{Gen}",200,0,1,0,2,{"DRN"}},
// };
//  vector<string> GEN = {"Angle between gen photons"};

// loghist = {"Angle between gen photons","Total clustered rechits energy","Total unclustered rechit enrergy","Total energy of clustered rechits","Total energy of unclustered rechits","sigma_iE_iE","sigma_iphi_iphi"};
// norm ={""};
//   sprintf(hname,"temp.root");
//   TFile* fout = new TFile(hname,"RECREATE");
 
//     n_files=1;  
    
//   for(int i_file=0; i_file<n_files;i_file++)
//     {      
//     for(int i=0; i<size(varName); i++){
        
//         int rebin = varName[i].intData; 
//         string xLabel = varName[i].str2;
//         double ymin = varName[i].double1;
//         double ymax = varName[i].double2;
//         double xmin = varName[i].double3;
//         double xmax = varName[i].double4;
//         vector<string> legend_texts = varName[i].str3;
//         string Name = varName[i].str4;
//         vector<string> VarName = varName[i].str1;
//         cout << VarName.size() << endl;
//         vector<TH1F*> hist_list;
// 	cout << "XLABEL "<<xLabel<<endl;
//         for (int j=0; j<VarName.size();j++){
           
// 	  sprintf(hist_name,"%s",VarName[j].c_str());
// 	  cout<<hist_name<<"\t"<<i<<"\t"<<i_file<<"\t"<<f[i_file]->GetName()<<endl;
          
// 	  TH1F* h_resp2 = (TH1F*)f[i_file]->Get(hist_name); // SR
// 	  h_resp2->GetXaxis()->SetTitle(xLabel.c_str());
// 	  cout<<"resp2 "<<h_resp2->Integral()<<"\t"<<rebin<<"\t"<<xmin<<"\t"<<xmax<<endl;
	  
// 	  //h_resp2->Rebin(rebin);
	 
	  
// 	  //h_resp2= setMyRange(h_resp2,xmin,xmax);
// 	  //setLastBinAsOverFlow(h_resp2);
	  
	  
// 	  hist_list.push_back(h_resp2); 
//        } 
//       string  Savename;
//     int gen = count(GEN.begin(),GEN.end(),Name);
//     int LOG = count(loghist.begin(), loghist.end(),Name);
//     int NORM= count(norm.begin(), norm.end(),Name);
//     if(gen){Savename = "GEN_";}
//     else {Savename = "RECO_";}
//          string Savename2=Savename+to_string(1000 + i)+"_" +Name + "_Fit";

// generate_1Dplot(hist_list,Savename2.c_str(),xLabel.c_str(),"Normalized frequency",rebin,ymin,ymax,xmin,xmax,leg_head,true,false,false,true,filetag[i_file].c_str(),legend_texts);
// if(LOG && NORM){
//          string Savename1=Savename +to_string(1000 + i) +"_"+Name + "_Fit";
// generate_1Dplot(hist_list,Savename1.c_str(),xLabel.c_str(),"Normalized frequency",rebin,ymin,ymax,xmin,xmax,leg_head,true,true,false,true,filetag[i_file].c_str(),legend_texts);
// }

//          //  else if(LOG && !NORM){
// //           string Savename1=Savename+to_string(1000 + i)+"_" +Name + "_Fit";
// // generate_1Dplot(hist_list,Savename1.c_str(),xLabel.c_str(),"Normalized frequency",rebin,ymin,ymax,xmin,xmax,leg_head,false,true,false,true,filetag[i_file].c_str(),legend_texts);
// // 	  }
// //           else if(!LOG && NORM){
// //           string Savename1=Savename+to_string(1000 + i)+"_" +Name + "_Fit"; 
// // generate_1Dplot(hist_list,Savename1.c_str(),xLabel.c_str(),"Normalized frequency",rebin,ymin,ymax,xmin,xmax,leg_head,true,false,false,true,filetag[i_file].c_str(),legend_texts);
// // }
	 
          
        
//     }
//     }
//}
