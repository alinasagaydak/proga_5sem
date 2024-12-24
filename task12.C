#include <fstream>
#include <vector>
#include <TH1D.h>

void FillDataFromFile (std::string name_file, std::vector<double>& data) {
    std::ifstream ifile_data(name_file);
    if (ifile_data.is_open()) {
        std::string line;
        double val;
        while (std::getline(ifile_data, line)) {
            ifile_data >> val;
            data.push_back(val);
        }
    }
}

double fitfunc(double* x, double* par) {
    return par[0] * sin(x[0] + par[1]) + par[2] + par[3] * TMath::Gamma(0.5, x[0]);
}

void task12() {
    gStyle->SetOptFit(1111);
    
    std::vector<double> data;
    FillDataFromFile("task10Nov.dat", data);
    TH1D* h = new TH1D("h", ";mm;yields / 0.1 mm", 100, 0., 10.);
    h->SetLineWidth(2);
    for (double val : data) h->Fill(val);

    auto fitf = new TF1("fitf", fitfunc, 0., 10., 4);
    fitf->SetParameter(0, 5.);
    fitf->SetParameter(1, TMath::Pi()); 
    fitf->SetParameter(2, 10.);
    fitf->SetParameter(3, 1.); 
    fitf->SetParNames("Ampl sine", "sine phase", "const", "Ampl Gamma");
    h->Fit(fitf, "PS");
    h->Draw();
}
