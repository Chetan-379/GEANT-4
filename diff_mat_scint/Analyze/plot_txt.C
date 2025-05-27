void plot_txt() {

  ifstream infile("plots/Scintillation/Fits/txtFiles/PbWO4/Efficiency_vs_Einc_PbWO4_5cm.txt");
  if (!infile.is_open()) {
    cerr << "Error: Could not open data.txt" << endl;
    return;
  }

  vector<double> x_vals, y_vals;
  double x, y;
  while (infile >> x >> y) {
    x_vals.push_back(x);
    y_vals.push_back(y);
  }

  infile.close();

  int nPoints = x_vals.size();
  if (nPoints == 0) {
    cerr << "Error: No data read from file." << endl;
    return;
  }

  TGraph* graph = new TGraph(nPoints, &x_vals[0], &y_vals[0]);
  graph->SetTitle("PbWO4_5cm");
  graph->SetMarkerStyle(20);
  graph->SetMarkerColor(kBlue);
  graph->SetLineColor(kRed);
  graph->SetMarkerSize(1.0);

  TCanvas* canvas = new TCanvas("canvas", "Efficiency", 800, 600);
  graph->Draw("AP");
    
  graph->GetXaxis()->SetTitle("E_{inc}");
  graph->GetXaxis()->SetTitleOffset(0.9);
  graph->GetXaxis()->SetTitleFont(42);

  graph->GetYaxis()->SetTitle("Efficiency");
  graph->GetYaxis()->SetTitleOffset(0.9);
    
  graph->GetYaxis()->SetTitleFont(42);
  graph->GetXaxis()->SetTitleSize(0.05);
  graph->GetYaxis()->SetTitleSize(0.05);

  graph->SetMinimum(0.8); 
  graph->SetMaximum(1.1);                 
        
  canvas->Update();
    
  canvas->SaveAs("Overlay_plots/image.png"); 
}
