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
    //auto f1 = new TF1("f1", "TMath::Gamma(par[3], x[0])", 0., 10., 1);
    //double C = 1000 + par[0] * (cos(10+par[1]) - cos(par[1])) - par[2] * f1->Integral(0., 10.);
    //return par[0] * sin(x[0] + par[1]) + C + par[2] * TMath::Gamma(par[3], x[0]);
    return par[0] + par[1] * sin(par[2] * (x[0] + par[3]));
}

void task12() {
    gStyle->SetOptFit(1111);
    
    std::vector<double> data;
    FillDataFromFile("task10Nov.dat", data);
    TH1D* h = new TH1D("h", ";mm;yields / 0.1 mm", 100, 0., 10.);
    h->SetLineWidth(2);
    for (double val : data) h->Fill(val);

    auto fitf = new TF1("fitf", fitfunc, 0., 10., 4);
    fitf->SetParameter(0, 10.);
    fitf->SetParameter(1, 1.); 
    fitf->SetParameter(2, 1.);
    fitf->SetParameter(3, TMath::Pi()/2.); 
    h->Fit(fitf, "PS");
    h->Draw();
}
