void Overlay_txt() {

  vector<string> Variables = {"Resolution", "Efficiency"};
  vector<string> Depths = {"2.5cm", "5cm", "10cm"};

  vector<float> ymax = {0.2, 1.1};
  vector<float> ymin = {0.0, 0.2};

  vector<float> xLegMin = {0.7, 0.15};   
  vector<float> xLegMax = {0.9, 0.35};
  
  vector<float> yLegMin = {0.7, 0.2};
  vector<float> yLegMax = {0.9, 0.4};

   
  string folder;

  for (int iVar = 0; iVar < Variables.size(); iVar++){
    for (int iDepth = 0; iDepth < Depths.size(); iDepth++){
      string FileName1 = "Chetan_Thesis_Plots/Trans_20cmX20cm/Scintillation/Fits/txtFiles/PbWO4/"+ Variables[iVar] +"_vs_Einc_PbWO4_"+ Depths[iDepth] +".txt";
      string FileName2 = "Chetan_Thesis_Plots/Trans_20cmX20cm/Scintillation/Fits/txtFiles/BGO/"+ Variables[iVar] +"_vs_Einc_BGO_"+ Depths[iDepth] +".txt";
      string FileName3 = "Chetan_Thesis_Plots/Trans_20cmX20cm/Scintillation/Fits/txtFiles/LYSO/"+ Variables[iVar] +"_vs_Einc_LYSO_"+ Depths[iDepth] +".txt";
      string FileName4 = "Chetan_Thesis_Plots/Trans_20cmX20cm/Scintillation/Fits/txtFiles/GAGG/"+ Variables[iVar] +"_vs_Einc_GAGG_"+ Depths[iDepth] +".txt";
      string FileName5 = "Chetan_Thesis_Plots/Trans_20cmX20cm/Scintillation/Fits/txtFiles/LaBr3/"+ Variables[iVar] +"_vs_Einc_LaBr3_"+ Depths[iDepth] +".txt";
 
      vector<string> filenames = {FileName1, FileName2, FileName3, FileName4, FileName5};

      //vector<string> filenames = {"PbWO4_emission.txt", "BGO_emission.txt", "LYSO_emission.txt", "GAGG_emission.txt", "LaBr3_emission.txt"};
      for (int i = 0; i< filenames.size(); i++){
	cout << "plotting the files: " << filenames[i] << endl;}

      vector<TGraph*> graphs;
      int colors[] = {kRed, kBlue, kGreen + 2, kBlack, kMagenta + 3};
      //int colors[] = {kRed, kBlue, kGreen + 2};
  
      for (size_t i = 0; i < filenames.size(); ++i) {
	ifstream infile(filenames[i]);
	if (!infile.is_open()) {
	  cerr << "Error opening file: " << filenames[i] << endl;
	  continue;
	}

	vector<double> x_vals, y_vals;
	double x, y, er;
        //
	if (iVar == 0) {
	  while (infile >> x >> y >> er) {
	  x_vals.push_back(x);
	  y_vals.push_back(y);
	  }
	}

	else if (iVar == 1) {
	  while (infile >> x >> y) {
	  x_vals.push_back(x);
	  y_vals.push_back(y);
	  }
	}

	infile.close();

	TGraph* g = new TGraph(x_vals.size(), &x_vals[0], &y_vals[0]);
	// TGraph* g = new TGraph(filenames[i], "%lg %lg");
	g->SetTitle(filenames[i].c_str());
	g->SetLineColor(colors[i % 5]);
	g->SetMarkerColor(colors[i % 5]);

	graphs.push_back(g);
      }

      // Create canvas and draw graphs
      string canvas_name = Variables[iVar] + Depths[iDepth];
      TCanvas* c = new TCanvas(canvas_name.c_str(), canvas_name.c_str(), 800, 600);

      c->cd();

      graphs[0]->SetMarkerStyle(22);
      graphs[0]->SetLineWidth(10);
      graphs[0]->SetMarkerSize(1.8);

      graphs[1]->SetMarkerStyle(24);
      graphs[1]->SetLineWidth(10);
      graphs[1]->SetMarkerSize(1.8);

      graphs[2]->SetMarkerStyle(23);
      graphs[2]->SetLineWidth(15);
      graphs[2]->SetMarkerSize(1.8);

      if (graphs[3]){
      graphs[3]->SetMarkerStyle(25);
      graphs[3]->SetLineWidth(15);
      graphs[3]->SetMarkerSize(1.8);
      }

      if (graphs[4]){
      graphs[4]->SetMarkerStyle(47);
      graphs[4]->SetLineWidth(10);
      graphs[4]->SetMarkerSize(1.8);
      }



      for (int i = 0; i < graphs.size(); ++i) {
	if (i==0){
	  graphs[i]->SetTitle(" ");
	  graphs[i]->SetMinimum(ymin[iVar]);
	  graphs[i]->SetMaximum(ymax[iVar]);

	  // graphs[i]->SetMinimum(0.00);
	  // graphs[i]->SetMaximum(0.02);

	  graphs[i]->GetXaxis()->SetTitle("E_{inc} (keV)");
	  graphs[i]->GetXaxis()->SetTitleSize(0.05);
	  graphs[i]->GetXaxis()->SetTitleOffset(0.9);
	  graphs[i]->GetXaxis()->SetTitleFont(42);

	  if (iVar == 0) graphs[i]->GetYaxis()->SetTitle((Variables[iVar] + "(#sigma/#mu)").c_str());
	  else graphs[i]->GetYaxis()->SetTitle(Variables[iVar].c_str());
	  graphs[i]->GetYaxis()->SetTitleSize(0.05);
	  graphs[i]->GetYaxis()->SetTitleOffset(1);    
	  graphs[i]->GetYaxis()->SetTitleFont(42);
     
	  graphs[i]->Draw("AP");
	  //graphs[i]->GetXaxis()->SetLimits(300, 700);;
	}

	else {
	  graphs[i]->Draw("P same");	
	}
	gPad->SetLeftMargin(0.12);
      }

      vector<string> legtxt = {"PbWO4", "BGO", "LYSO", "GAGG", "LaBr3"};
    
      TLegend* legend = new TLegend(xLegMin[iVar], yLegMin[iVar], xLegMax[iVar], yLegMax[iVar]);
      for (size_t i = 0; i < graphs.size(); ++i) {
	legend->AddEntry(graphs[i], legtxt[i].c_str(), "p");
      }
    
      legend->Draw();

      string MainFolder = "Chetan_Thesis_Plots/Trans_20cmX20cm/Scintillation/Fits/Resol_Resp/";
      //string MainFolder = "test_plots/";	
	string FolderName = MainFolder + "Summary_Plots/";

      vector<string> diff_title = {"Resolution", "Efficiency"};

      string ImgName;
      ImgName = FolderName + diff_title[iVar] + "_Material_Overlay_" + Depths[iDepth] + ".png";
      //cout << ImgName << endl;
      c->SaveAs(ImgName.c_str());

      ImgName = FolderName + diff_title[iVar] + "_Material_Overlay_" + Depths[iDepth] + ".pdf";
      //cout << ImgName << endl;
      c->SaveAs(ImgName.c_str());
    }
  }
}
