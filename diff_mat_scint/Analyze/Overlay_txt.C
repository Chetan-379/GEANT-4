void Overlay_txt() {
  string folder;
  
  //string material = "PbWO4";
  //string material = "PbWO4";

  // folder = "plots/Scintillation/Fits/txtFiles/" + material + "/";
  // cout << "plotting txtfiles in folder: "<< folder << endl;
  
 
  //std::vector<std::string> filenames = {folder + "Efficiency_vs_Einc_" + material + "_2.5cm.txt", folder + "Efficiency_vs_Einc_" + material + "_5cm.txt", folder + "Efficiency_vs_Einc_" + material +"_10cm.txt"};

  //std::vector<std::string> filenames = {folder + "Efficiency_vs_Einc_" + material + "_2.5cm.txt", folder + "Efficiency_vs_Einc_" + material + "_5cm.txt", folder + "Efficiency_vs_Einc_" + material +"_10cm.txt"};
  std::vector<std::string> filenames = {"plots/Scintillation/Fits/txtFiles/PbWO4/Efficiency_vs_Einc_PbWO4_2.5cm.txt", "plots/Scintillation/Fits/txtFiles/PbWO4/Efficiency_vs_Einc_PbWO4_5cm.txt", "plots/Scintillation/Fits/txtFiles/PbWO4/Efficiency_vs_Einc_PbWO4_10cm.txt"};
  //std::vector<std::string> filenames = {"plots/Scintillation/Fits/txtFiles/PbWO4/Efficiency_vs_Einc_PbWO4_5cm.txt", "plots/Scintillation/Fits/txtFiles/PbWO4/Efficiency_vs_Einc_PbWO4_5cm.txt"};

  for (int i = 0; i< filenames.size(); i++){cout << "plotting the files: " << filenames[i] << endl;}
  std::vector<TGraph*> graphs;

    // Colors for each graph
    int colors[] = {kRed, kBlue, kGreen + 2};

    // Loop over each file
    for (size_t i = 0; i < filenames.size(); ++i) {
        std::ifstream infile(filenames[i]);
        if (!infile.is_open()) {
            std::cerr << "Error opening file: " << filenames[i] << std::endl;
            continue;
        }

        std::vector<double> x_vals, y_vals;
        double x, y;
        while (infile >> x >> y) {
            x_vals.push_back(x);
            y_vals.push_back(y);
        }

        infile.close();

        TGraph* g = new TGraph(x_vals.size(), &x_vals[0], &y_vals[0]);
        g->SetTitle(filenames[i].c_str());
        g->SetLineColor(colors[i % 3]);
        g->SetMarkerColor(colors[i % 3]);
        // g->SetMarkerStyle(4);
        // g->SetLineWidth(2);
	// g->SetMarkerSize(1.5);

        graphs.push_back(g);
    }

    // Create canvas and draw graphs
    TCanvas* c = new TCanvas("c", "Overlay Graphs", 800, 600);
    c->cd();

    bool first = true;
    
    //for (auto g : graphs) {
    //if (first) {

    graphs[0]->SetMarkerStyle(22);
    graphs[0]->SetLineWidth(10);
    graphs[0]->SetMarkerSize(1.5);

    graphs[1]->SetMarkerStyle(24);
    graphs[1]->SetLineWidth(10);
    graphs[1]->SetMarkerSize(1.5);

    graphs[2]->SetMarkerStyle(23);
    graphs[2]->SetLineWidth(10);
    graphs[2]->SetMarkerSize(1.5);



    for (int i = 0; i < graphs.size(); ++i) {
      if (i==0){
	//graphs[i]->SetTitle("Efficiency_PbWO4_vs_PbWO4_5cm");
	graphs[i]->SetTitle("Efficiency_PbWO4");
	// graphs[i]->SetMinimum(0.5);                 // pb efficiency
	// graphs[i]->SetMaximum(1.1);                
	//graphs[i]->SetMinimum(0.075);                 // pb resolution
	//graphs[i]->SetMaximum(0.25);                 
	graphs[i]->SetMinimum(0.55);                 // bgo resolution
	graphs[i]->SetMaximum(1.1);                 

	graphs[i]->GetXaxis()->SetTitle("E_{inc}");
	graphs[i]->GetXaxis()->SetTitleSize(0.05);
	graphs[i]->GetXaxis()->SetTitleOffset(0.9);
	graphs[i]->GetXaxis()->SetTitleFont(42);

	//graphs[i]->GetYaxis()->SetTitle("Efficiency");
	graphs[i]->GetYaxis()->SetTitle("Efficiency");
	graphs[i]->GetYaxis()->SetTitleSize(0.05);
	graphs[i]->GetYaxis()->SetTitleOffset(0.9);    
	graphs[i]->GetYaxis()->SetTitleFont(42);

	

	graphs[i]->Draw("AP"); // A=axis, P=points, L=line
	first = false;
      } else {
	//graphs[i]->SetTitle("");
	graphs[i]->Draw("P same");
      }
    }

    vector<string> legtxt = {"2.5cm", "5cm", "10cm"};
    //vector<string> legtxt = {"PbWO4", "PbWO4"};
    // Add a legend
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    for (size_t i = 0; i < graphs.size(); ++i) {
      legend->AddEntry(graphs[i], legtxt[i].c_str(), "p");
    }
    legend->Draw();

    // Save the canvas
    //c->SaveAs(string("Overlay_plots/Efficiency_overlay_width_"+material+".png").c_str());      //efficiency
    c->SaveAs(string("Overlay_plots/image.png").c_str());      //resolution
}
