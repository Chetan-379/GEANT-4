const int n_pl = 4;
bool logx = false;
//defining the legends for each plots
TString legend_text[10] = {"Integral","PbWO4","BGO", "CZT", "pMSSM_MCMC_106_19786","pMSSM_MCMC_473_54451"};
int line_width[12] = {2,2,2,2,2,2,2,2,2,2,2,2};
int line_style[12] = {1,1,1,1,1,1,1,1,1,1,1,1};
int line_color[9] = {kBlue, kBlack, kRed, kGreen+2, kMagenta, kAzure + 7 , kCyan + 1 , kGreen + 3 };
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
//void generate_1Dplot(vector<TH1F*> hist, TH1* hist_ratio,char const *tag_name="", float energy=-1, int xmax=-1,int xmin=-1,char const *leg_head="",
//		     bool normalize=false, bool log_flag=true, bool DoRebin=false, bool save_canvas=true, char const *title="", const char* xtitile="", int rebin=-1){  
void generate_1Dplot(vector<TH1*> hist, char const *tag_name="", int xmax=-1,int xmin=-1,char const *leg_head="",
		     bool normalize=false, bool log_flag=true, bool DoRebin=false, bool save_canvas=true, char const *title="", const char* xtitile="", const char* ytitile="", int rebin=-1, bool statFlag = false){  

  TCanvas *canvas_n1 = new TCanvas(tag_name, tag_name,900,750);//600,600,1200,1200);
  canvas_n1->Range(-60.25,-0.625,562.25,0.625);
  canvas_n1->SetFillColor(0);
  canvas_n1->SetBorderMode(0);
  canvas_n1->SetBorderSize(2);
  gStyle->SetOptStat(statFlag);
  

  auto *p1 = new TPad("p1","p1",0.,0.01,1.,1.);  p1->Draw();
  p1->SetBottomMargin(0.12);  
  p1->SetRightMargin(0.035);
  p1->SetLeftMargin(0.13);
  p1->SetTopMargin(1);

  p1->cd();
  
  vector<TString> legName;
  std::string leg_head_str = leg_head;
  double x = 0.15;
  double y = 0.90;
  TLegend *legend; //legend to be drawn on the plot - shift x,ys if you want to move this on the canvas
  legend = new TLegend(0.37,0.76,0.45,0.89);  
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

  TLegendEntry* leg_entry[11];
  TLegendEntry* leg_entry1[11];
  float x_label_size = 0.045;
  double ymin = 100000.0;
  double ymax = 0.0;
  double xrange = xmax;

  cout<<" hist.size() = "<<hist.size()<<endl;

  for(int i = 0; i < (int)hist.size(); i++) {

    char final_label[100];
    float Integral = hist.at(i)->Integral();     
    sprintf(final_label,"%s (%0.2f)",legend_text[i].Data(), Integral); 

    if(normalize) {
      hist.at(i)->Scale(1.0/hist.at(i)->Integral());
      hist.at(i)->GetYaxis()->SetTitle("Normalized");
    }
    else {
      hist.at(i)->GetYaxis()->SetTitle(ytitile);
    }

    //getting the value of xmax:
    int nbins = hist.at(i)->GetNbinsX();
    int lastNonZeroBin = 0;
    
    for (int ibin = nbins; ibin >= 1; --ibin) {
      //if (hist.at(i)->GetBinContent(ibin) > 0) {
      if (hist.at(i)->GetBinContent(ibin) >= (hist.at(i)->GetMaximum())*0.001) {
        lastNonZeroBin = ibin;
        break;
      }
    }

    //double effectiveXmax = hist.at(i)->GetBinLowEdge(lastNonZeroBin + 1);

    hist.at(i)->GetXaxis()->SetRangeUser(xmin, xmax); //to be given while calling this function
    //hist.at(i)->GetXaxis()->SetRangeUser(xmin, effectiveXmax);
    hist.at(i)->SetLineWidth(line_width[i]); //these are defined on top of this script    
    hist.at(i)->SetLineStyle(line_style[i]);
 hist.at(i)->SetLineColor(line_color[i]);
 hist.at(i)->SetTitle(xtitile); //you can change it if you want
 //setLastBinAsOverFlow(hist.at(i),0);
    //
    hist.at(i)->GetXaxis()->SetTitleSize(0.05);
    hist.at(i)->GetXaxis()->SetLabelSize(x_label_size);
    hist.at(i)->GetXaxis()->SetLabelSize(0.00450);
    //    hist.at(i)->GetXaxis()->SetRange(0,1500);
    hist.at(i)->GetYaxis()->SetTitleSize(0.05);
    hist.at(i)->GetYaxis()->SetLabelSize(0.05);
    hist.at(i)->GetYaxis()->SetTitleOffset(1.0);
    hist.at(i)->GetYaxis()->SetLabelSize(x_label_size);
    //hist.at(i)->GetYaxis()->SetLabelFormat("6.3E");
    hist.at(i)->GetXaxis()->SetTitle(xtitile); //setting the title of X axis

    hist.at(i)->GetXaxis()->SetTitleSize(0.05);
    hist.at(i)->GetXaxis()->SetLabelSize(0.04);
    hist.at(i)->GetYaxis()->SetLabelSize(0.04);
    hist.at(i)->GetYaxis()->SetTitleSize(0.05);
    hist.at(i)->GetYaxis()->SetTitleOffset(1.1);
    hist.at(i)->GetXaxis()->SetTitleOffset(1.1);

    //gStyle->SetTitleFontSize(0.0);
    
    TLatex* textOnTop = new TLatex();
    textOnTop->SetTextSize(0.04);

    //textOnTop->DrawLatexNDC(0.13,0.96,"CMS it{#bf{Simulation Preliminary}}");
    //textOnTop->DrawLatexNDC(0.13,0.96,"it{#bf{GEANT4 Simulation}}");
  
    sprintf(lhead,"%s ",title);
    //textOnTop->DrawLatexNDC(0.79,0.96,lhead);
    textOnTop->DrawLatexNDC(0.79,0.4,lhead);
    char* en_lat = new char[500];
    textOnTop->SetTextSize(0.035);
    if(DoRebin) { //if rebin flag is on - this will reduce the bin size by half
      hist.at(i)->Rebin(rebin);
    }
    //setting up the legend style and all

    
    //leg_entry[i] = legend->AddEntry(hist.at(i),(hist.at(i)->Integral().c_str()),"xsl");

    legName.push_back(hist.at(i)->GetName());
    leg_entry[i] = legend->AddEntry(hist.at(i), final_label, "l");
    
    //leg_entry1[i] = legend1->AddEntry(hist.at(i),final_label,"l");
    //leg_entry[i] = legend->AddEntry(hist.at(i),(hist.at(i)->Integral().c_str()),"xsl");
    leg_entry[i]->SetTextColor(hist.at(i)->GetLineColor());
    //leg_entry1[i]->SetTextColor(hist.at(i)->GetLineColor());
    if(hist.at(i)->GetMaximum() > ymax) ymax = hist.at(i)->GetMaximum();
    if(hist.at(i)->GetMinimum() < ymin) ymin = hist.at(i)->GetMinimum();
    //setLastBinAsOverFlow(hist.at(i),0);
    
  }
  

  if(ymin == 0.0) ymin = 1e-3;
  if(ymin<0.0) ymin = 1e-4;
  //  if(ymax<=10) ymax=10;
  for(int i = 0; i < (int)hist.size(); i++) {
    if(!normalize && !log_flag) hist.at(i)->GetYaxis()->SetRangeUser(ymin*100,ymax+(ymax/7));
    if(!normalize && log_flag) hist.at(i)->GetYaxis()->SetRangeUser(1,ymax*10);
    else
      //hist.at(i)->GetYaxis()->SetRangeUser(ymin/10,ymax*5.0);
      hist.at(i)->GetYaxis()->SetRangeUser(ymin/10,ymax+(ymax/7));
    //    p1->SetGrid();
    
    // if(!i) hist.at(i)->Draw("hist colz");
    // else   hist.at(i)->Draw("hist sames colz"); //overlaying the histograms

    if(!i) hist.at(i)->Draw("hist colz");
    else   hist.at(i)->Draw("hist colz"); //overlaying the histograms
    
    TPaveStats *ptstats1 = new TPaveStats(0.80,0.75,0.965,0.9,"brNDC");
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
    //ptstats1->SetTitle("fsdfa");
    //resp1->Scale(1./resp1->Integral());

    //TGaxis::SetMaxDigits(0);
    // TPaletteAxis *palette = (TPaletteAxis*)hist.at(i)->GetListOfFunctions()->FindObject("palette");
  //   if (palette) {
  //       // Set new position in NDC (normalized device coordinates: 0 to 1)
  //       palette->SetX1NDC(0.55); // left side
  //       palette->SetX2NDC(0.58); // right side
  //       palette->SetY1NDC(0.2);  // bottom
  //       palette->SetY2NDC(0.5);  // top
  //   }

  //   canvas_n1->Modified(); // Mark canvas as modified
  //   canvas_n1->Update();   // Redraw with new palette position
    
  }
  
  //legend->Draw();
  //legend1->Draw();


  if(log_flag) 
    gPad->SetLogy();
  
  if(logx)
    gPad->SetLogx();

  gPad->Update();
  canvas_n1->Update();

   
 
  TLatex* textOnTop = new TLatex();
  textOnTop->SetTextSize(0.04);
  //textOnTop->DrawLatexNDC(0.12,0.91," #it{#bf{GEANT4 Simulation}}");

  TLatex* MatName = new TLatex();
  MatName->SetTextSize(0.035);
  MatName->DrawLatexNDC(0.50,0.83, title);
 
  // char* en_lat = new char[500];
  // textOnTop->SetTextSize(0.04);
  // float inlumi=energy;
  // sprintf(en_lat,"#bf{%0.2f fb^{-1} (13 TeV)}",inlumi);
  // textOnTop->DrawLatexNDC(0.7,0.91,en_lat);

 

  
  char* canvas_name = new char[1000];
  //c->Print(canvas_name);
  
  if(save_canvas) {
    sprintf(canvas_name,"%s.png",tag_name);//.png",tag_name);//_wnormalize.png",tag_name);
    canvas_n1->SaveAs(canvas_name);   
    // sprintf(canvas_name,"%s.pdf",tag_name);
    // canvas_n1->SaveAs(canvas_name);
    
  }
}

const int nfiles=100,nBG=6;                                                                                                                                                              


void plot_all_hists(string pathname)
{
  char *hist_name = new char[100];
  char *full_path = new char[1000];
  char *leg_head = new char[100];
  vector<string> f;
  string filetag, folder, Energy, Width;

  // f = {"out_root_files/Si_FluON_10000evts_out.root", "out_root_files/PbWO4_FluON_10000evts_out.root", "out_root_files/BGO_FluON_10000evts_out.root"};
  // filetag={"Si", "PbWO4", "BGO"};
  // folder = {"plots/Si", "plots/PbWO4", "plots/BGO"};

  
  
  // string path = "./out_root_files";  // current directory
  // vector<string> rootFile;
  //   for (const auto& entry : fs::directory_iterator(path)) {
  //       if (entry.is_regular_file()) {
  //           std::string filename = entry.path().filename().string();
  //           if (filename.find(".root") != std::string::npos) {
  //               rootFile.push_back(filename);
  //           }
  //       }
  //   }

   std::vector<std::string> rootFiles;

   string FileFolder;
   FileFolder = "BGO_out_root_files";
   TSystemDirectory dir(FileFolder.c_str(), FileFolder.c_str());
   TList* files = dir.GetListOfFiles();
   if (!files) return;

    TIter next(files);
    TObject* obj;

    while ((obj = next())) {
      TString fname = obj->GetName();
        if (!obj->IsA()->InheritsFrom("TSystemFile")) continue;
        if (fname.EndsWith(".root")) {
	  f.push_back(FileFolder + "/" + string(fname.Data()));
        }
    }

    cout << "\n\nSize of the files vector is: " << f.size() << endl;

  
    // filetag={material};

  // vector<TString> Energy, Material, Widths;
  // vector<double> hist_Max;
  // vector<int> Energy_val;
  
  // Material = {"BGO", "PbWO4", "Plastic"};
  // Widths = {"2.5cm", "5cm", "10cm", "20cm"};

  
  // Energy = {"511keV", "100keV", "150keV", "300keV", "450keV", "600keV", "800keV", "1000keV"};
  // Energy_val = {511, 100, 150, 300, 450, 600, 800, 1000};

  //folder = {"plots/BGO/PhoInteract"};
  //folder = {"plots"};

  vector<int >rebin = {1,1,1,1,1,1,1,1,1}; //keep it 1 if you don't want to change hist bins
      
  //x axis range
  vector<int>xmax = {700,2000,16,16,2000,2000,2000,2000};
  vector<int>xmin = {300,0,0,0,0,0,0,0};

      
  for(int iFile=0; iFile < f.size(); iFile++){ //looping over each file    		    
    TFile *root_file = new TFile(f[iFile].c_str());

    // Loop over all objects in the file
    vector<TH1F> HistList;
    TIter keyList(root_file ->GetListOfKeys());
    TKey* key;

    while ((key = (TKey*)keyList())) {      
      TClass *cl = gROOT->GetClass(key->GetClassName());
      if (!cl->InheritsFrom("TH1")) continue;
      TH1 *hist = (TH1*)key->ReadObj();

      string histType = key->GetClassName();
      vector<TH1*> hist_list;

      TString histId = hist->GetName();

      //if(!(histId.Contains("Edep") || histId.Contains("Compt"))) continue;
      hist_list.push_back(hist);
      //if (hist_list.size() > 0) cout << hist_list[0]->GetName() << endl;
      int xrange=0.0;	  
	
      //setting up the axis titles
      string xtitle, ytitle;
      bool hist_2d = false;
      if ( histType == "TH2I" || histType == "TH2F" || histType == "TH2D") {
	xtitle = hist->GetXaxis()->GetTitle();
	ytitle = hist->GetYaxis()->GetTitle();
	hist_2d = true;}

      else {
	xtitle = hist->GetName();
	ytitle = "Entries";}

      //specifying the folder to save the plots
      string MainFolder;
      bool logFlag = false;
      bool statFlag = false;
      if(histId.Contains("Edep") || histId.Contains("Compt")) {MainFolder = "Pho_Interact"; if (!hist_2d) logFlag = true;}
      else {MainFolder = "Scintillation/Opt_Pho_Ana"; statFlag = true;}
      TString FileId = f[iFile];

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
	
      }
      
      if(FileId.Contains("PbWO4")){
	if (FileId.Contains("2.5cm")) {folder = "plots/" + MainFolder + "/PbWO4/2.5cm"; Width = "2.5cm";}
	if (FileId.Contains("_5cm")) {folder = "plots/" + MainFolder + "/PbWO4/5cm"; Width = "5cm";}
	if (FileId.Contains("10cm")) {folder = "plots/" + MainFolder + "/PbWO4/10cm"; Width = "10cm";}

	filetag = "PbWO4_" + Energy + "_" + Width;

	if (histId.Contains("nOptical_Photons")) xmax[0] = 750;
      }
      
      if(FileId.Contains("Plastic")){
	filetag = "Plastic_" + Energy + "_" + Width;
	if (FileId.Contains("20cm")) {folder = "plots/" + MainFolder + "/Plastic/20cm"; Width = "20cm";}
	if (FileId.Contains("5cm")) {folder = "plots/" + MainFolder + "/Plastic/5cm"; Width = "5cm";}
	if (FileId.Contains("10cm")) {folder = "plots/" + MainFolder + "/Plastic/10cm"; Width = "10cm";}

	filetag = "Plastic_" + Energy + "_" + Width;
      }
      
      //path to save the files a jpg or pdf
      sprintf(full_path,"%s%s/%s_%s",pathname.c_str(),folder.c_str(), xtitle.c_str(),filetag.c_str());

      //calling generate_1Dplot which will take this vector of histograms 
      generate_1Dplot(hist_list,full_path,xmax[0],xmin[0],leg_head,false,logFlag,false,true,filetag.c_str(),xtitle.c_str(),ytitle.c_str(),rebin[0], statFlag);
    }
  }
}




