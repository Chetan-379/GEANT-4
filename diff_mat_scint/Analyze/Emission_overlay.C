const int n_pl = 4;
bool logx = false;
//defining the legends for each plots
TString legend_text[15] = {"PbWO_{4}","BGO", "LYSO", "GAGG", "LaBr_{3}"};

int line_width[15] = {2,2,2,2,2,2,2,2,2,2,2,2,2,2};
int line_style[15] = {1,1,1,1,1,1,1,1,1,1,1,1,2,2};
int line_color[15] = {kBlue, kGreen+1, kPink+5, kOrange-3, kCyan+4, kAzure+7, kSpring+8, kPink+10, kGreen + 3,kRed,kGray, kViolet -7 };

TH1F* setLastBinAsOverFlow(TH1F*, int);
TH1F* setMyRange(TH1F*,double,double);
TH1F* DrawOverflow(TH1F*);

TH1F* DrawOverflow(TH1F* h,int xmin, int xrange){
  //function to paint the histogram h with an extra bin for overflows
  // This function paint the histogram h with an extra bin for overflows
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
  // Restore the number of entries
  htmp->SetEntries(h->GetEntries());
  // FillStyle and color
  // htmp->SetFillStyle(h->GetFillStyle());
  // htmp->SetFillColor(h->GetFillColor());
  return htmp;
}


TH1F* setLastBinAsOverFlow(TH1F* h_hist, int xrange){
  //     h_hist = setMyRange(h_hist,0,xrange);
  //  h_hist->GetXaxis()->SetRangeUser(0,xrange);
  double lastBinCt =h_hist->GetBinContent(h_hist->GetNbinsX()),overflCt =h_hist->GetBinContent(h_hist->GetNbinsX());
  //  cout<<h_hist->GetNbinsX()<<"\t"<<lastBinCt<<"\t"<<overflCt<<endl;

  double lastBinErr=h_hist->GetBinError(h_hist->GetNbinsX()),  overflErr=h_hist->GetBinError(h_hist->GetNbinsX()+1); 
  if(lastBinCt!=0 && overflCt!=0)
    lastBinErr = (lastBinCt+overflCt)* (sqrt( ((lastBinErr/lastBinCt)*(lastBinErr/lastBinCt)) + ((overflErr/overflCt)*(overflErr/overflCt)) ) );

  else if(lastBinCt==0 && overflCt!=0)
    lastBinErr = overflErr;
  else if(lastBinCt!=0 && overflCt==0)
    lastBinErr = lastBinErr;
  else lastBinErr=0;
  //h_temp->GetXaxis()->SetRangeUser(0,xrange);

  lastBinCt = lastBinCt+overflCt;
  //  cout<<lastBinCt<<endl;
  TH1F* h_temp = (TH1F*)h_hist->Clone();
  h_temp->SetBinContent(h_hist->GetNbinsX(),lastBinCt);
  h_temp->SetBinError(h_hist->GetNbinsX(),lastBinErr);
  //  h_temp->GetXaxis()->SetRangeUser(0,xrange);

  // h_hist = setMyRange(h_hist,0,xrange);
  //
  return h_temp;
}


TH1F* setMyRange(TH1F *h1,double xLow,double xHigh){
  //call it after setting last bin as overflow                                                                                                    
  double err=0;
  if(xHigh > 13000) return h1;
  if(xLow < -13000) return h1;

  // h1->Print("all");
  //  h1->GetXaxis()->SetRangeUser(xLow,xHigh);  
  int nMax=h1->FindBin(xHigh);
  h1->SetBinContent(nMax,h1->IntegralAndError(nMax,h1->GetNbinsX(),err));
  h1->SetBinError(nMax,err);

  //  cout<<nMax<<endl;
  //  h1->GetXaxis()->SetRangeUser(xLow,xHigh);
  for(int i=nMax+1;i<=h1->GetNbinsX()+1;i++){
    h1->SetBinContent(i,0);
    h1->SetBinError(i,0);
    //    cout<<":";
    //h1->GetXaxis()->SetRangeUser(xLow,xHigh); 
  }
  return h1;
}

//important function - change in this function if you want to change decoration of a plot
void generate_1Dplot(vector<TH1F*> hist, char const *tag_name="", double xmax=-1.,double xmin=-1.,char const *leg_head="",
		     bool normalize=false, bool log_flag=true, bool DoRebin=false, bool save_canvas=true, char const *title="", const char* xtitile="", int rebin=-1){  
  

  TCanvas *canvas_n1 = new TCanvas(tag_name, tag_name,900,750);//600,600,1200,1200);
  canvas_n1->Range(-60.25,-0.625,562.25,0.625);
  canvas_n1->SetFillColor(0);
  canvas_n1->SetBorderMode(0);
  canvas_n1->SetBorderSize(2);
  canvas_n1->SetRightMargin(0.045);
  canvas_n1->SetLeftMargin(0.12);
  canvas_n1->SetTopMargin(0.06);
  canvas_n1->SetBottomMargin(0.12);

  gStyle->SetOptStat(0);
  
  vector<TString> legName;
  std::string leg_head_str = leg_head;
  double x = 0.15;
  double y = 0.90;
  TLegend *legend; //legend to be drawn on the plot - shift x,ys if you want to move this on the canvas
  legend = new TLegend(0.3,0.7,0.5,0.9);  
  legend->SetTextSize(0.030);
  legend->SetLineColor(kWhite);
  char* lhead = new char[100];

  sprintf(lhead,"#bf{%s} ",title); // legend header -- for example you are plotting this for WGJets then title string can be "WGJets"
  legend->SetHeader(lhead);
  legend->SetLineColor(kWhite);

  TLegendEntry* leg_entry[15];
  TLegendEntry* leg_entry1[15];
  float x_label_size = 0.045;
  double ymin = 100000.0;
  double ymax = 0.0;
  cout<<" hist.size() = "<<hist.size()<<endl;


  for(int i = 0; i < (int)hist.size(); i++) {
    if(normalize) {
      hist.at(i)->Scale(1.0/hist.at(i)->Integral());
      hist.at(i)->GetYaxis()->SetTitle("Normalized");
    }
    else {
      hist.at(i)->Scale(1.0/hist.at(i)->GetMaximum());
      hist.at(i)->GetYaxis()->SetTitle("Entries");
    }

    hist.at(i)->Rebin(rebin);
    hist.at(i)->GetXaxis()->SetRangeUser(xmin, xmax); //to be given while calling this function
    hist.at(i)->SetLineWidth(line_width[i]); //these are defined on top of this script    
    hist.at(i)->SetLineStyle(line_style[i]);
    hist.at(i)->SetLineColor(line_color[i]);
    hist.at(i)->SetTitle(" "); //you can change it if you want
    //setLastBinAsOverFlow(hist.at(i),0);
    //
    hist.at(i)->GetXaxis()->SetTitleSize(0.05);
    hist.at(i)->GetXaxis()->SetLabelSize(x_label_size);
    hist.at(i)->GetXaxis()->SetLabelSize(0.0450);
    //    hist.at(i)->GetXaxis()->SetRange(0,1500);
    hist.at(i)->GetYaxis()->SetTitleSize(0.055);
    hist.at(i)->GetYaxis()->SetLabelSize(0.05);
    hist.at(i)->GetYaxis()->SetTitleOffset(1.2);
    hist.at(i)->GetYaxis()->SetLabelSize(x_label_size);
    hist.at(i)->GetXaxis()->SetTitle(xtitile); //setting the title of X axis

    hist.at(i)->GetXaxis()->SetTitleSize(0.05);
    hist.at(i)->GetXaxis()->SetLabelSize(0.04);
    hist.at(i)->GetYaxis()->SetLabelSize(0.04);
    hist.at(i)->GetYaxis()->SetTitleSize(0.05);
    hist.at(i)->GetYaxis()->SetTitleOffset(1.1);
    hist.at(i)->GetXaxis()->SetTitleOffset(1.1);

    TLatex* textOnTop = new TLatex();
    textOnTop->SetTextSize(0.04);

    textOnTop->DrawLatexNDC(0.13,0.96,"CMS it{#bf{Simulation Preliminary}}");
    sprintf(lhead,"%s ",title);
    textOnTop->DrawLatexNDC(0.79,0.96,lhead);
    char* en_lat = new char[500];
    textOnTop->SetTextSize(0.035);

    //setting up the legend style and all
    vector <string> Intg;
    char final_label[100];    
    float Integral = hist.at(i)->Integral();
    sprintf(final_label,"%s",legend_text[i].Data());
    //sprintf(final_label,"%.0f events", Integral);
   
    legName.push_back(hist.at(i)->GetName());
    leg_entry[i] = legend->AddEntry(hist.at(i), final_label, "l");
    leg_entry[i]->SetTextColor(hist.at(i)->GetLineColor());
    if(hist.at(i)->GetMaximum() > ymax) ymax = hist.at(i)->GetMaximum();
    if(hist.at(i)->GetMinimum() < ymin) ymin = hist.at(i)->GetMinimum();
    
  }
  if(ymin == 0.0) ymin = 1e-3;
  if(ymin<0.0) ymin = 1e-4;
  for(int i = 0; i < (int)hist.size(); i++) {
    if(!normalize && log_flag) hist.at(i)->GetYaxis()->SetRangeUser(1.0,ymax*10);
    else if(!normalize && !log_flag) hist.at(i)->GetYaxis()->SetRangeUser(0,ymax+(ymax/7));
    
    else
      hist.at(i)->GetYaxis()->SetRangeUser(0.00001,ymax*60.0);
    
    if(!i) hist.at(i)->Draw("hist");
    else   hist.at(i)->Draw("hist sames"); //overlaying the histograms
	
  }
  
  legend->Draw();

  if(log_flag) gPad->SetLogy();
  if(logx) gPad->SetLogx();
  
  gPad->Update();
  
  char* canvas_name = new char[1000];
  //saving the file
  if(save_canvas) {
    sprintf(canvas_name,"%s.png",tag_name);
    canvas_n1->SaveAs(canvas_name);       
  }
  
}

void overlay(string pathname)
{
  char* hname = new char[200];
  char* hist_name  = new char[200];
  char* full_path = new char[2000];
  char* path2 = new char[2000];
  char* title= new char[2000];
  char *leg_head = new char[200];

  vector<string> f;
 
  f= {"PbWO4_out_root_files/PbWO4_5cm_ScintON_511keV_Birk0_NoAbs_spline_out.root", "BGO_out_root_files/new_BGO_5cm_ScintON_511keV_Birk0_NoAbs_spline_out.root", "LYSO_out_root_files/LYSO_5cm_ScintON_511keV_Birk0_NoAbs_spline_out.root", "GAGG_out_root_files/GAGG_5cm_ScintON_511keV_Birk0_NoAbs_spline_out.root", "LaBr3_out_root_files/LaBr3_5cm_ScintON_511keV_Birk0_NoAbs_spline_out.root"};

  //f= {"BGO_out_root_files/BGO_5cm_ScintON_100keV_Birk0_NoAbs_NoCeren_SefDerEmiss_spline_out.root", "BGO_out_root_files/BGO_5cm_ScintON_511keV_Birk0_NoAbs_NoCeren_SefDerEmiss_spline_out.root", "BGO_out_root_files/BGO_5cm_ScintON_1000keV_Birk0_NoAbs_NoCeren_SefDerEmiss_spline_out.root"};
  
  //define your histograms to be read from here
  string histname1[100], histname2[100], histname3[100];

  vector<string> variables;
  //variables = {"Total_Edep"};//, "Compton_Edep", "Total_Edep","Photo_Edep"};   //put the names of variables here
  variables = {"OptPho_lmbda"};
  //variables = {"nOptical_Photons"};

 
  vector<string> filetag;
  filetag={" "};
 
  //rebin values
  vector<int >rebin = {1,1,1,1,1,1,1,1,1,1,1,1}; //keep it 1 if you don't want to change hist bins
      
  //x axis range
  vector<double>xmax = {750};
  vector<double>xmin = {250};

  for (int iVar =0; iVar < variables.size(); iVar++){
    vector<TH1F*> hist_list;
    for(int i_file=0; i_file < f.size(); i_file++) //looping over each file
      { TFile *root_file = new TFile(f[i_file].c_str()); 	    

	sprintf(hist_name,"%s", variables[iVar].c_str());
	cout<<"i_file "<<i_file<<"\t"<<iVar<<"\t"<<root_file->GetName()<< "\t" << hist_name<<endl;
	TH1F* resp = (TH1F*)root_file->Get(hist_name); //reading hist from the TFile
	hist_list.push_back(resp);
      }
	
    
    int xrange=0.0;

    //x axis title for all plots///
    vector<string>diff_title;

    //diff_title = { "Compton_Edep"};
    //diff_title = { "nOptPhotons"};
    //diff_title = { "OptPho_lmbda (nm)"};
    diff_title = { "\\lambda (nm)"};
    vector<string>xtitle;
    xtitle = {diff_title[iVar]};
	 
    //path to save the files a jpg or pdf
    vector<string> folder;   
    folder = {"Overlay_plots/", "Overlay_plots/", "Overlay_plots/"};
    sprintf(full_path,"%s/%simage",pathname.c_str(),folder[iVar].c_str());
	 	    
    //calling generate_1Dplot which will take this vector of histograms and 
    generate_1Dplot(hist_list,full_path,xmax[iVar],xmin[iVar],leg_head,false,false,false,true,filetag[iVar].c_str(),xtitle[iVar].c_str(),rebin[iVar]);
  }
}

















































