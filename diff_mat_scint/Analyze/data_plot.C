
//defining fit function for resolution
double sqrtE(double *x, double *par) {
 
  double a = par[0];   // coefficient
  double n = par[1];  //power
  double b = par[2];  // any constant
 
  return (sqrt((pow(a/pow(x[0],n), 2)) + (pow(b,2))));
}

void data_plot() {
  vector<string> data_files;
  vector<string> xtitle, ytitle;

  data_files = {"Emeasured_vs_Einc.txt", "Ngamma_vs_Einc.txt", "Response_vs_Einc.txt", "Resolution_vs_Einc.txt"};
  xtitle = {"E_{incident}(keV)", "E_{incident}(keV)", "E_{incident}(keV)", "E_{incident}(keV)"};
  ytitle = {"E_{Measured}(keV)", "N_{\\gamma}", "\\frac{E_{Measured}}{E_{incident}}", "\\frac{\\sigma_{N_{\\gamma}}}{\\mu_{N_{\\gamma}}}"};

  for (int i=0; i<4; i++){
    std::ifstream infile(data_files[i]);
    if (!infile.is_open()) {
        std::cerr << "Could not open data.txt!" << std::endl;
        return;
    }

    std::vector<double> x_vals, y_vals, x_err_vals, y_err_vals;
    double x, y, y_err;

    while (infile >> x >> y >> y_err) {
        x_vals.push_back(x);
        y_vals.push_back(y);
	x_err_vals.push_back(0);
	y_err_vals.push_back(y_err);
    }

    if (x_vals.empty()) {
        std::cerr << "No data found!" << std::endl;
        return;
    }

    // Create canvas
     data_files[i].resize(data_files[i].length() - 4);        //removing .txt part 
     TCanvas *canvas = new TCanvas(data_files[i].c_str(), data_files[i].c_str(), 800, 600);
     canvas->SetLeftMargin(0.135);
     canvas->SetRightMargin(0.035);
     canvas->SetTopMargin(0.057);
     canvas->SetBottomMargin(0.11);


    // char *xtitle = new char[50];
    // char *ytitle = new char[50];

    //sprintf(xtitle, %s, "");

    // Create and draw the graph
     //TGraph *graph = new TGraph(x_vals.size(), &x_vals[0], &y_vals[0]);

    auto graph = new TGraphErrors(x_vals.size(), &x_vals[0], &y_vals[0], &x_err_vals[0], &y_err_vals[0]);

    graph->SetTitle(" ");
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.0);
    graph->SetMarkerColor(kBlue);
    graph->Draw("AP");
    // graph->GetXaxis()->SetTitle("E_{True}");
    // graph->GetYaxis()->SetTitle("N_{\\gamma}");

    graph->GetXaxis()->SetTitle(xtitle[i].c_str());
    graph->GetXaxis()->SetTitleSize(0.055);
    graph->GetXaxis()->SetTitleOffset(0.8);
    
    graph->GetYaxis()->SetTitle(ytitle[i].c_str());
    graph->GetYaxis()->SetTitleSize(0.055);
    graph->GetYaxis()->SetTitleOffset(0.9);

    // Fit with a straight line (pol1) and draw
    if (data_files[i] == "Emeasured_vs_Einc" || data_files[i] == "Ngamma_vs_Einc"){
    TF1 *fit = new TF1("fit", "pol1", *std::min_element(&x_vals[0], &x_vals[0]+x_vals.size()), *std::max_element(&x_vals[0], &x_vals[0]+x_vals.size()));
    graph->Fit(fit);  // "Q" = quiet, no print to stdout
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
    }

    if (data_files[i] == "Response_vs_Einc"){
      graph->GetYaxis()->SetRangeUser(0.5,1.2);
      TLine *line = new TLine(10, 1, 1090, 1);
      line->SetLineColor(kRed);
      line->SetLineWidth(2);
      line->SetLineStyle(2);
      line->Draw("same");
      canvas->Update();  
    }

    if (data_files[i] == "Resolution_vs_Einc"){      
      TF1 *fit_func = new TF1("fit_fcn", sqrtE, 100, 1000,3);
      fit_func->SetParameters(1,1/2,0);
      
      //TF1 *fit_func = new TF1("fit_fcn", sqrtE, 100, 1000,2);
      //fit_func->SetParameters(1,0);
      fit_func->SetRange(100, 1000);

      graph->Fit(fit_func, "E");  

      fit_func->SetLineColor(kRed);
      fit_func->Draw("same");

      double a = fit_func->GetParameter(0); 
      double n = fit_func->GetParameter(1); 
      double b = fit_func->GetParameter(2); 
     
    // Create equation string
    TString eq;
    eq.Form("y = %.2fx^{-%.2f} \\oplus %.2f", a, n, b);
    //eq.Form("y = %.2fx^{-\\frac{1}{2}} \\oplus %.2f", a, b);

    char* chisqrNDF = new char[50];
    sprintf(chisqrNDF, "\\chi^{2}/NDF: %e", fit_func->GetChisquare()/fit_func->GetNDF());

    TLatex latex;
    latex.SetNDC();   
    latex.SetTextSize(0.045);
    latex.DrawLatex(0.25, 0.85, eq);
    latex.DrawLatexNDC(0.25,0.78, chisqrNDF);
    canvas->Update();    
    }

    //saving the canvas
    string image_name = "./plots/ScintPlots/" +  data_files[i] + ".png";
    canvas->SaveAs(image_name.c_str());
  }
}
