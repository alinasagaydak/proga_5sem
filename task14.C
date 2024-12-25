#include <fstream>
#include <TH1D.h>
#include <TCanvas.h>
#include <TVirtualFFT.h>
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

void task14() {
    std::vector<double> data;
    FillDataFromFile("dataFFT.dat", data);
    
    auto myc = new TCanvas("myc", "", 1200, 800);
    myc->Divide(2, 2);

    myc->cd(1);
    auto hsig = new TH1D("hsig","", 100, 0., 1.);
    for (double val : data) hsig->Fill(val);
    hsig->SetLineWidth(2);
    hsig->Draw();

    myc->cd(2);
    TH1* hfft = 0;
    hfft = hsig->FFT(hfft, "MAG");
    hfft->Draw();
}
