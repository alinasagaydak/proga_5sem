#include <vector>
#include <TH1D.h>
#include <TCanvas.h>
#include <TF1Convolution.h>

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

double function2(double* x, double* par) {
    return TMath::Gaus(x[0], par[0], par[1], kTRUE) + TMath::Gaus(x[0], par[2], par[3], kTRUE); 
}

void task13() {
    std::vector<double> data;
    FillDataFromFile("m3piJPSI_cut.dat", data);
    auto hist = new TH1D("hist", "", 400, 3., 3.2);
    for (double val : data) hist->Fill(val);
    auto c1 = new TCanvas();
    c1->SetLogy();
    hist->SetLineWidth(2);
    hist->Draw();

    TF1 *f2 = new TF1("f2", function2, 3., 3.2, 4);
    f2->FixParameter(2, 0.6);
    TF1Convolution* Conv = new TF1Convolution(f2, f2, 3., 3.2, true);
    Conv->SetNofPointsFFT(1000);
    TF1 *f = new TF1("NFFT = 1000", *Conv, 3., 3.2, Conv->GetNpar());
    f->SetNpx(10000);f->Draw("same");
}
