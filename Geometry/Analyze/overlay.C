const int n_pl = 4;
bool logx = false;
//defining the legends for each plots
//TString legend_text[15] = {"\\gamma_{E}: 511keV","511keV", "1000keV"};
TString legend_text[15] = {"Sim: 100k Evts","511keV", "1000keV"};

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
  

  //TCanvas *canvas_n1 = new TCanvas(tag_name, tag_name,900,750);//600,600,1200,1200);
  TCanvas *canvas_n1 = new TCanvas(tag_name, tag_name,1000,750);//600,600,1200,1200);
  canvas_n1->Range(-60.25,-0.625,562.25,0.625);
  //canvas_n1->Range(-70.25,-0.625,300.25,0.625);
  canvas_n1->SetFillColor(0);
  canvas_n1->SetBorderMode(0);
  canvas_n1->SetBorderSize(2);
  //canvas_n1->SetRightMargin(0.045);
  canvas_n1->SetRightMargin(0.095);
  canvas_n1->SetLeftMargin(0.12);
  canvas_n1->SetTopMargin(0.06);
  canvas_n1->SetBottomMargin(0.12);

  
  gStyle->SetOptStat(1);
  
  vector<TString> legName;
  std::string leg_head_str = leg_head;
  double x = 0.15;
  double y = 0.90;
  TLegend *legend; //legend to be drawn on the plot - shift x,ys if you want to move this on the canvas
  legend = new TLegend(0.40,0.8,0.60,0.9);
  //legend = new TLegend(0.15,0.8,0.35,0.9);  
  //legend->SetTextSize(0.030);
  legend->SetTextSize(0.040);
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
      hist.at(i)->GetYaxis()->SetTitle("Entries");
      //hist.at(i)->GetYaxis()->SetTitle("E_{\\gamma_{scattered}}");
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
    leg_entry[i] = legend->AddEntry(hist.at(i), final_label, " ");
    leg_entry[i]->SetTextColor(hist.at(i)->GetLineColor());
    if(hist.at(i)->GetMaximum() > ymax) ymax = hist.at(i)->GetMaximum();
    if(hist.at(i)->GetMinimum() < ymin) ymin = hist.at(i)->GetMinimum();
    
  }
  if(ymin == 0.0) ymin = 1e-3;
  if(ymin<0.0) ymin = 1e-4;
  for(int i = 0; i < (int)hist.size(); i++) {
    //if(!normalize) hist.at(i)->GetYaxis()->SetRangeUser(1.0,ymax*10);
    if(!normalize && log_flag) hist.at(i)->GetYaxis()->SetRangeUser(1, ymax*10);
    if(!normalize && !log_flag) hist.at(i)->GetYaxis()->SetRangeUser(0,ymax+(ymax/7));
    
    // else
    //   hist.at(i)->GetYaxis()->SetRangeUser(0.00001,ymax*60.0);
    
    //if(!i) hist.at(i)->Draw("colz");
    if(!i) hist.at(i)->Draw("hist");
    else   hist.at(i)->Draw("hist sames"); //overlaying the histograms
    //canvas_n1->Update();

    
    //TPaveStats *ptstats1 = new TPaveStats(0.8,0.83,0.955,0.94,"brNDC");
    TPaveStats *ptstats1 = new TPaveStats(0.71,0.78,0.905,0.94,"brNDC");
    ptstats1->SetBorderSize(1);
    ptstats1->SetFillColor(0);
    ptstats1->SetLineColor(hist.at(i)->GetLineColor());
    ptstats1->SetTextAlign(12);
    ptstats1->SetTextColor(hist.at(i)->GetLineColor());
    ptstats1->SetTextFont(42);
    ptstats1->SetOptStat(110011);
    ptstats1->SetOptFit(1);
    hist.at(i)->GetListOfFunctions()->Add(ptstats1);
    ptstats1->SetParent(hist.at(i));

  }
  
  legend->Draw();

  if(log_flag) gPad->SetLogy();
  if(logx) gPad->SetLogx();
  

  gPad->Update();
  canvas_n1->Update();

  
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
 
  // f= {"test_comp_egEta_out.root"};
  //f= {"test_comp_egEta_out.root"};
  //f= {"test_comp_egEta_nHit2plus_posTheta_100kEvts_z0_out.root"};
  //f= {"test_diffEta_out.root"};
  //f= {"test_diffEta_check_out_nHit_updte.root"};
  //f= {"test_BinWidth_change_out.root"};
  //f= {"random_mom_100kEvts_out.root"};
  //f= {"random_pol_random_mom_100kEvts_out.root"};
  f= {"mom_randomXYZ_pol_random_100kEvts_out.root"};

  
  //define your histograms to be read from here
  string histname1[100], histname2[100], histname3[100];

  vector<string> nHitCat;
  nHitCat = {"all_inc", "nHit0", "nHit1_S", "nHit1_B", "nHit2_S1S1", "nHit2_S1S2", "nHit2_SB", "nHit2_B1B1", "nHit2_B1B2", "nHit2_BS", "remain"};
  //nHitCat = {"all_inc"};//, "nHit0", "nHit1_S", "nHit1_B", "nHit2_S1S1", "nHit2_S1S2", "nHit2_SB", "nHit2_B1B1", "nHit2_B1B2", "nHit2_BS", "remain"};

  char hname_theta[100], hname_eta[100], hname_Eout[100], hname_ThetaVsEout[100]; 

  
  for(int i=0; i< nHitCat.size(); i++){
    sprintf(hname_theta, "h_theta_%s",nHitCat[i].c_str());
    sprintf(hname_eta, "h_eta_%s",nHitCat[i].c_str());
    sprintf(hname_Eout, "h_Eout_%s",nHitCat[i].c_str());
    //sprintf(hname_ThetaVsEout, "h_theta_Eout_%s",nHitCat[i].c_str());
    
    

    vector<string> variables;
    variables = {hname_theta, hname_eta, hname_Eout};//, hname_ThetaVsEout};
    //variables = {hname_ThetaVsEout};
  

  //vector<string> variables;  
  //variables = {"Eta_compSig", "nHits", "scat_theta", "HitTime", "HitPosX", "HitPosY", "HitPosZ"};
  //variables = {"pol_scat_eta", "nHits", "scat_theta", "HitTime", "HitPosX", "HitPosY", "HitPosZ", "ThetaVsEta"};
  //variables = {"h_theta_all_inc", "h_eta_all_inc", "h_Eout_all_inc"};

  
 
  //vector<string> filetag;
    vector<string> filetag;
    //filetag={"all inc", "nHit==0", "nHit==1_(Plastic)", "nHit==1_(BGO)", "nHit==2_Plastic_singleCrys", "nHit==2_Plastic_diffCryst", "nHit==2_Plastic_BGO" ,"nHit==2_BGO_singleCrys", "nHit==2_BGO_diffCryst", "nHit==2_BGO_Plastic", "rest (nHit>2)"};
    filetag={"all inc", "nHit==0", "nHit==1_(Plastic)", "nHit==1_(BGO)", "nHit>=2_Plastic_singleCrys", "nHit>=2_Plastic_diffCryst", "nHit>=2_Plastic_BGO" ,"nHit>=2_BGO_singleCrys", "nHit>=2_BGO_diffCryst", "nHit>=2_BGO_Plastic", "rest"};

    
  //filetag={"scatTheta inclusive", "nHits", "Compt_Scat_theta", "Hit_Time", "HitPosX", "HitPosY", "HitPosZ", "Theta vs Eta"};
  //filetag={"all_inc", "all_inc", "all_inc"};//, "Hit_Time", "HitPosX", "HitPosY", "HitPosZ", "Theta vs Eta"};
    //filetag = nHitCat[i]; 
 
  //rebin values
  vector<int >rebin = {1,1,1,1,1,1,1,1,1,1,1,1}; //keep it 1 if you don't want to change hist bins
      
  //x axis range
  vector<double>xmax = {200, 200, 200, 300, 700, 700, 700};
  vector<double>xmin = {0, -200, 0, 0, 0, -700, -700};

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
    //diff_title = {"eta(#circ)", "nHits", "theta(#circ)", "time(ns)", "posX(mm)", "posY(mm)", "posZ(mm)", "eta(#circ)"};
    diff_title = {"#theta_{scat}(#circ)", "eta(#circ)", "Escat_(MeV)", "test"};//, "posX(mm)", "posY(mm)", "posZ(mm)", "eta(#circ)"};
    // vector<string>xtitle;
    // xtitle = {diff_title[0], diff_title[1]};
	 
    //path to save the files a jpg or pdf
    vector<string> folder;   
    folder = {"plots/", "plots/", "plots/", "plots/", "plots/", "plots/", "plots/", "plots/1"};
    sprintf(full_path,"%s%s%s",pathname.c_str(),folder[iVar].c_str(),variables[iVar].c_str());

    cout << full_path << endl;
    bool logFlag = false;
    //if (filetag[iVar] == "nHits" || filetag[iVar] == "Hit_Time" || filetag[iVar] == "HitPosX" || filetag[iVar] == "HitPosY" || filetag[iVar] == "HitPosZ") logFlag = true;
    //if (filetag[iVar] == "Hit_Time" || filetag[iVar] == "HitPosX" || filetag[iVar] == "HitPosY" || filetag[iVar] == "HitPosZ") logFlag = true;
	 	    
    //calling generate_1Dplot which will take this vector of histograms and 
    generate_1Dplot(hist_list,full_path,xmax[iVar],xmin[iVar],leg_head,false,logFlag,false,true,filetag[i].c_str(),diff_title[iVar].c_str(),rebin[iVar]);
  }
  }
}


















































