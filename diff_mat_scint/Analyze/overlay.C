const int n_pl = 4;
bool logx = false;
//defining the legends for each plots
TString legend_text[15] = {"PbWO4","Si", "BGO"};//,"NHadJet>=2","ST>300","Trig_Eff","Event_clean","dPhi(MET,leadingJet)","remove overlap","MET>200","e/mu_veto","charge_trk_veto", "pMSSM_MCMC_106_19786","pMSSM_MCMC_473_54451"};
int line_width[15] = {2,2,2,2,2,2,2,2,2,2,2,2,2,2};
int line_style[15] = {1,1,1,1,1,1,1,1,1,1,1,1,2,2};
int line_color[15] = {kBlue, kGreen+1, kPink+5, kOrange-3, kCyan+4, kAzure+7, kSpring+8, kPink+10, kGreen + 3,kRed,kGray, kViolet -7 };
//int line_color[15] = {kBlack, kRed+1, kGreen, kMagenta+2, kOrange+1, kAzure+7, kCyan+4, kSpring+9, kPink+10, kGreen+3, kBlue, kGray+1, kViolet-7};
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
void generate_1Dplot(vector<TH1F*> hist, char const *tag_name="", float energy=-1, int xmax=-1,int xmin=-1,char const *leg_head="",
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
  // auto *pad_1 = new TPad("pad_1","pad_1",0.,0.0,1.,0.32); pad_1->Draw();
  // pad_1->SetTopMargin(0.06);
  // pad_1->SetBottomMargin(0.3);
  // pad_1->SetRightMargin(0.025);
  // pad_1->SetLeftMargin(0.14);
  
  // auto *p1 = new TPad("p1","p1",0.,0.32,1.,1.);  p1->Draw();
  // p1->SetBottomMargin(0.04);
  // p1->SetRightMargin(0.025);
  // p1->SetLeftMargin(0.14);
  // p1->SetTopMargin(0.05);
  // p1->cd();
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

  // TLegend *legend1; //legend to be drawn on the plot - shift x,ys if you want to move this on the canvas
  // legend1 = new TLegend(0.65,0.55,0.79,0.8);  
  // legend1->SetTextSize(0.030);
  // legend1->SetLineColor(kWhite);
  //char* lhead = new char[100];
  
  sprintf(lhead,"#bf{%s} ",title); // legend header -- for example you are plotting this for WGJets then title string can be "WGJets"
  legend->SetHeader(lhead);
  legend->SetLineColor(kWhite);

  //legend1->SetHeader(lhead);
  //legend1->SetLineColor(kWhite);

  TLegendEntry* leg_entry[15];
  TLegendEntry* leg_entry1[15];
  float x_label_size = 0.045;
  double ymin = 100000.0;
  double ymax = 0.0;
  cout<<" hist.size() = "<<hist.size()<<endl;


  for(int i = 0; i < (int)hist.size(); i++) {
    // if(DoRebin) {
    //  hist.at(i)->Rebin(2);

    // }
    //hist.at(i)= setLastBinAsOverFlow(hist.at(i),xrange);
     

    // normalize = true;

    // Int_t NEvnts;
    // float Evnt_frac;
    // Int_t no_of_events;
    
    // no_of_events = hist.at(i)->GetEntries();
    // NEvnts = hist.at(i)->Integral();
    // Evnt_frac = (no_of_events/(hist.at(0)->GetEntries()))*100;
    // int overflCt= (hist.at(i)->GetBinContent(hist.at(i)->GetNbinsX()+1));
    
    // cout << no_of_events << "    (frac: " << setprecision(3) << Evnt_frac << "%)" << overflCt << endl;
    if(normalize) {
      	hist.at(i)->Scale(1.0/hist.at(i)->Integral());
      	hist.at(i)->GetYaxis()->SetTitle("Normalized");
    }
    else {
      hist.at(i)->GetYaxis()->SetTitle("Entries");
       }
    hist.at(i)->GetXaxis()->SetRangeUser(xmin,xmax); //to be given while calling this function
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
    hist.at(i)->GetYaxis()->SetTitleSize(0.05);
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
    if(DoRebin) { //if rebin flag is on - this will reduce the bin size by half
     hist.at(i)->Rebin(rebin);
    }
    //setting up the legend style and all
       vector <string> Intg;
    char final_label[100];    
    float Integral = hist.at(i)->Integral();
    sprintf(final_label,"%s (%.1f)",legend_text[i].Data(), Integral);
   
    legName.push_back(hist.at(i)->GetName());
    leg_entry[i] = legend->AddEntry(hist.at(i), final_label, "l");
    leg_entry[i]->SetTextColor(hist.at(i)->GetLineColor());
    if(hist.at(i)->GetMaximum() > ymax) ymax = hist.at(i)->GetMaximum();
    if(hist.at(i)->GetMinimum() < ymin) ymin = hist.at(i)->GetMinimum();
    
  }
  if(ymin == 0.0) ymin = 1e-3;
  if(ymin<0.0) ymin = 1e-4;
  for(int i = 0; i < (int)hist.size(); i++) {
    if(!normalize) hist.at(i)->GetYaxis()->SetRangeUser(1.0,ymax*10);
    else
      hist.at(i)->GetYaxis()->SetRangeUser(0.00001,ymax*60.0);
    //    p1->SetGrid();
    
    if(!i) hist.at(i)->Draw("hist");
    else   hist.at(i)->Draw("hist sames"); //overlaying the histograms
	
  }
  
  legend->Draw();
  //legend1->Draw();


  if(log_flag) 
    gPad->SetLogy();
  
  if(logx)
    gPad->SetLogx();
  
  gPad->Update();

 
  TLatex* textOnTop = new TLatex();
  textOnTop->SetTextSize(0.04);
  //textOnTop->DrawLatexNDC(0.12,0.96,"CMS #it{#bf{Preliminary}}");
  textOnTop->DrawLatexNDC(0.12,0.96,"#bf{Work in Progress}");
  
  char* en_lat = new char[500];
  textOnTop->SetTextSize(0.04);
  float inlumi=energy;
  sprintf(en_lat,"#bf{%0.2f fb^{-1} (13 TeV)}",inlumi);
  //textOnTop->DrawLatexNDC(0.7,0.96,en_lat);

 

  char* canvas_name = new char[1000];
  //c->Print(canvas_name);
  //saving the file
  if(save_canvas) {
    sprintf(canvas_name,"%s.png",tag_name);
    canvas_n1->SaveAs(canvas_name);   
    // sprintf(canvas_name,"%s.pdf",tag_name);
    // canvas_n1->SaveAs(canvas_name);
    
  }
  
}
const int nfiles=100,nBG=6;                                                                                                                                                              
//TFile *f[nfiles];


void overlay(string pathname)
{
  char* hname = new char[200];
  char* hist_name  = new char[200];
  char* hist_name1 = new char[200];
  char* hist_name2 = new char[200];
  char* hist_name3 = new char[200];
  char* hist_name4 = new char[200];
  char* hist_name5 = new char[200];
  char* hist_name6 = new char[200];
  char* hist_name7 = new char[200];
  char* full_path = new char[2000];
  char* full_path1 = new char[2000];
  char* full_path2 = new char[2000];
  char* path2 = new char[2000];
  char* title= new char[2000];
  char* full_path3 = new char[2000];
  char* full_path4 = new char[2000];
  char* full_path5 = new char[2000];
  char* full_path6 = new char[2000];
  char* full_path7 = new char[2000];
  char* full_path8 = new char[2000];
  char* full_path9 = new char[2000];
  char* full_path10 = new char[2000];
  char* full_path11= new char[2000];
  char *leg_head = new char[200];
  //define your files here
  
  vector<string> f;
  f = {"out_root_files/Gamma_Interaction_PbWO4_PPT_PbWO4_out.root", "out_root_files/Gamma_Interaction_Si_PPT_Si_out.root", "out_root_files/Gamma_Interaction_BGO_PPT_BGO_out.root"};
  //define your histograms to be read from here
  int n_files=f.size(); //you have n files in this example
 
  string histname1[100], histname2[100], histname3[100];
  char hname_NHadJets[100], hname_Jet_Pt[100], hname_Jet_Eta[100], hname_Jet_Phi[100], hname_Met[100], hname_PhoPt[100], hname_PhoEta[100], hname_PhoPhi[100];
  // Book your histograms & summary counters here
  // vector<string> selection = {"Total_Edep"};

  // for (int i=0; i<selection.size();i++)
  //   {
  //     sprintf(hname_NHadJets,"h_NHadJets_%s",selection[i].c_str());
  //     histname3[i] = hname_NHadJets;
  //     sprintf(hname_Jet_Pt, "h_Jet_Pt_%s",selection[i].c_str());
  //     sprintf(hname_Jet_Eta, "h_Jet_Eta_%s",selection[i].c_str());
  //     sprintf(hname_Jet_Phi, "h_Jet_Phi_%s",selection[i].c_str());
  //     sprintf(hname_Met, "h_MET_%s",selection[i].c_str());
  //     histname2[i] = hname_Met; 
  //     sprintf(hname_PhoPt, "h_Pho_Pt_%s",selection[i].c_str());
  //     histname1[i] = hname_PhoPt; 


  //     sprintf(hname_PhoEta, "h_Pho_Eta_%s",selection[i].c_str());
  //     sprintf(hname_PhoPhi, "h_Pho_Phi_%s",selection[i].c_str());
  //  }

  
  vector<string> variables, baseline2, baseline3;
  // for (int i=0; i < selection.size(); i++){
  // baseline1.push_back(histname1[i]);
  // baseline2.push_back(histname2[i]);
  // baseline3.push_back(histname3[i]); 
  //}

  variables = {"Total_Edep", "Compton_Edep", "Photo_Edep"};   //put the names of variables here
  //baseline2 = 
  //string to be added to output file name - useful when you have different files and reading the same histograms from these
  vector<vector<string>> bigbaseline;
  //bigbaseline = {baseline1}; //, baseline2, baseline3};
 
  // for (int bigi=0; bigi<bigbaseline.size(); bigi++)
  //   {
      vector<string> filetag;
      filetag={"PBWO4","Si", "BGO"};//, "WGJets_2018", "ZJetsToNuNu_2018", "ZNuNuGJets_2018", "QCD_2018","GJets_2018"};

      //luminosity for each year - depends if you want to use it or - generate1Dplot uses this number and add it on the top
      vector<float>energyy;
      energyy={59.74,41.529};//,16.5,137.19,59.74,41.529,16.5,137.19,59.74,41.529,16.5,137.19,59.74,41.529,16.5,137.19,59.74,41.529,16.5,137.19,59.74,41.529,16.5,137.19,59.74,41.529,16.5,137.19,19.5,19.5,19.5,19.5,19.5,19.5,19.5,36.0,36.0,36.0,36.0,36.0,36.0,36.0};

      //rebin values
      vector<int >rebin = {2,2,2,2,2,2,2,2,2,2,2,2}; //keep it 1 if you don't want to change hist bins
      
      //x axis title for all plots///
      // vector<string>xtitle = {title[bigi], title[bigi]};

      //x axis range
      vector<int>xmax = {2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000};
      vector<int>xmin = {0,0,0,0,0,0,0,0,0,0,0};

      //looping over each files///
      for (int iVar =0; iVar < variables.size(); iVar++){
	vector<TH1F*> hist_list;
	for(int i_file=0; i_file < f.size(); i_file++) //looping over each file
	  { TFile *root_file = new TFile(f[i_file].c_str()); 
	    
	    //for(int i_cut=0; i_cut<bigbaseline[bigi].size();i_cut++) //looping over different histograms which should be overlayed on the same canvas and these histograms are saved in the same file
	    //{
	    //sprintf(hist_name,"%s",baseline[i_cut].c_str());
	    sprintf(hist_name,"%s", variables[iVar].c_str());
	       cout<<"i_file "<<i_file<<"\t"<<iVar<<"\t"<<root_file->GetName()<< "\t" << hist_name<<endl;
	      TH1F* resp = (TH1F*)root_file->Get(hist_name); //reading hist from the TFile
	      hist_list.push_back(resp);
	  }
	
	  float energy=energyy[0];
	  int xrange=0.0;

	   //x axis title for all plots///
	  vector<string>diff_title;
	  diff_title = { "Total_Edep", "Compton_Edep", "Photo_Edep"}; // , "Pt_Miss", "NHadJets"};
	  vector<string>xtitle;
	  xtitle = {diff_title[iVar], diff_title[iVar], diff_title[iVar]}; //, diff_title[bigi], diff_title[bigi], diff_title[bigi], diff_title[bigi], diff_title[bigi],diff_title[bigi],diff_title[bigi],diff_title[bigi]};
	 
	  //path to save the files a jpg or pdf
	  vector<string> folder;
	  folder = {"plots/PPT/", "plots/PPT/", "plots/PPT/"}; //, "plots/WG_/", "plots/ZNuNu_/", "plots/ZGNuNu_/", "plots/QCD_/","plots/GJets_/"};
	  sprintf(full_path,"%s/%s%s_Overlay_%s",pathname.c_str(),folder[iVar].c_str(),diff_title[iVar].c_str(),filetag[iVar].c_str());
	 	    
	  //calling generate_1Dplot which will take this vector of histograms and 
	  generate_1Dplot(hist_list,full_path,energy,xmax[iVar],xmin[iVar],leg_head,false,true,false,true,filetag[iVar].c_str(),xtitle[iVar].c_str(),rebin[iVar]);
      
	  //}

      }
    }



















































