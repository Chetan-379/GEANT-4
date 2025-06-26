void Overlay_txt() {
  string folder;
  
  //  vector<string> filenames = {"plots/Scintillation/Fits/txtFiles/PbWO4/Efficiency_vs_Einc_PbWO4_5cm.txt", "plots/Scintillation/Fits/txtFiles/LaBr3/Efficiency_vs_Einc_LaBr3_5cm.txt", "plots/Scintillation/Fits/txtFiles/GAGG/Efficiency_vs_Einc_GAGG_5cm.txt", "plots/Scintillation/Fits/txtFiles/LaBr3/Efficiency_vs_Einc_LaBr3_5cm.txt"};

  vector<string> filenames = {"plots/Scintillation/Fits/txtFiles/LaBr3/Efficiency_vs_Einc_LaBr3_1cm.txt", "plots/Scintillation/Fits/txtFiles/LaBr3/Efficiency_vs_Einc_LaBr3_2.5cm.txt", "plots/Scintillation/Fits/txtFiles/LaBr3/Efficiency_vs_Einc_LaBr3_5cm.txt", "plots/Scintillation/Fits/txtFiles/LaBr3/Efficiency_vs_Einc_LaBr3_10cm.txt"};

  for (int i = 0; i< filenames.size(); i++){
    cout << "plotting the files: " << filenames[i] << endl;}

  vector<TGraph*> graphs;
  int colors[] = {kRed, kBlue, kGreen + 2, kBlack};
  //int colors[] = {kRed, kBlue, kGreen + 2};

  for (size_t i = 0; i < filenames.size(); ++i) {
        ifstream infile(filenames[i]);
        if (!infile.is_open()) {
            cerr << "Error opening file: " << filenames[i] << endl;
            continue;
        }

        vector<double> x_vals, y_vals;
        double x, y;
        while (infile >> x >> y) {
            x_vals.push_back(x);
            y_vals.push_back(y);
        }

        infile.close();

        TGraph* g = new TGraph(x_vals.size(), &x_vals[0], &y_vals[0]);
        g->SetTitle(filenames[i].c_str());
        g->SetLineColor(colors[i % 4]);
        g->SetMarkerColor(colors[i % 4]);

        graphs.push_back(g);
    }

    // Create canvas and draw graphs
    TCanvas* c = new TCanvas("c", "Overlay Graphs", 800, 600);

    c->cd();

    graphs[0]->SetMarkerStyle(22);
    graphs[0]->SetLineWidth(10);
    graphs[0]->SetMarkerSize(1.5);

    graphs[1]->SetMarkerStyle(24);
    graphs[1]->SetLineWidth(10);
    graphs[1]->SetMarkerSize(1.5);

    graphs[2]->SetMarkerStyle(23);
    graphs[2]->SetLineWidth(10);
    graphs[2]->SetMarkerSize(1.5);

    graphs[3]->SetMarkerStyle(25);
    graphs[3]->SetLineWidth(10);
    graphs[3]->SetMarkerSize(1.8);



    for (int i = 0; i < graphs.size(); ++i) {
      if (i==0){
	graphs[i]->SetTitle("Efficiency_LaBr3");
	graphs[i]->SetMinimum(0.);
	graphs[i]->SetMaximum(1.1);                 

	graphs[i]->GetXaxis()->SetTitle("E_{inc}");
	graphs[i]->GetXaxis()->SetTitleSize(0.05);
	graphs[i]->GetXaxis()->SetTitleOffset(0.9);
	graphs[i]->GetXaxis()->SetTitleFont(42);

	graphs[i]->GetYaxis()->SetTitle("Efficiency");
	graphs[i]->GetYaxis()->SetTitleSize(0.05);
	graphs[i]->GetYaxis()->SetTitleOffset(0.9);    
	graphs[i]->GetYaxis()->SetTitleFont(42);
     
	graphs[i]->Draw("AP");
      }

      else {
	graphs[i]->Draw("P same");
      }
    }

    vector<string> legtxt = {"1cm", "2.5cm", "5cm", "10cm"};
    //vector<string> legtxt = {"2.5cm", "5cm", "10cm"};
    //vector<string> legtxt = {"PbWO4", "LaBr3", "GAGG", "LaBr3"};
    
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    for (size_t i = 0; i < graphs.size(); ++i) {
      legend->AddEntry(graphs[i], legtxt[i].c_str(), "p");
    }
    
    legend->Draw();

    c->SaveAs(string("Overlay_plots/image.png").c_str());
}
