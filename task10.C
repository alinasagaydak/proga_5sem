#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>

#include <TH1D.h>
#include <TCanvas.h>
#include <TMinuit.h>
#include <TF1.h>

//Global variables
double hist1[100];
double hist2[100];
double error1[100];
double error2[100];
double x[100];
double chi2;

double func(double x, double *par) {
    // ampl - кол-во сигнала
    return par[0] * TMath::Gaus(x, par[1], par[2]);
}

void fcn(int &npar, double *gin, double &f, double *par, int iflag) {
    const int nbins = 100;

//calculate chisquare
    double chisq = 0;
    chi2 = 0;
    for (int i = 0; i < nbins; i++) {
        double delta1  = (hist1[i] - func(x[i], par) - par[3]) / error1[i];
        double delta2 = (hist2[i] - par[3]) / error2[i];
        chisq += delta1 * delta1 + delta2 * delta2;
        chi2 += chisq;
    }
    chi2 /= 100.;
    f = chisq;
}


void FillDataFromFile (std::string name_file1, std::string name_file2, std::vector<double>& data1, std::vector<double>& data2) {
    std::ifstream ifile_data1(name_file1);
    if (ifile_data1.is_open()) {
        std::string line;
        double val;
        while (std::getline(ifile_data1, line)) {
            ifile_data1 >> val;
            data1.push_back(val);
        }
    }

    std::ifstream ifile_data2(name_file2);
    if (ifile_data2.is_open()) {
        std::string line;
        double val;
        while (std::getline(ifile_data2, line)) {
            ifile_data2 >> val;
            data2.push_back(val);
        }
    }
}

void task10() {
    // Отрисовка гистограммы
    std::string name_file1 ("data_1.dat");
    std::string name_file2("data_2.dat");
    std::vector<double> data1;
    std::vector<double> data2;
    FillDataFromFile(name_file1, name_file2, data1, data2);


    TH1D* h1 = new TH1D("h1", "Data 1", 100, 500, 600);
    TH1D* h2 = new TH1D("h2", "Data 2", 100, 500, 600);
    for (double val : data1) h1->Fill(val); 
    for (double val : data2) h2->Fill(val);

// Минимизация хи-квадрат
    // Init values
    for (int i = 0; i < 100; i++) {
        hist1[i] = h1->GetBinContent(i);
        hist2[i] = h2->GetBinContent(i);
    }
    
    for (int i = 0; i < 100; i++) x[i] = i + 500;
    
    //The errors values
    for (int i = 0; i < 100; i++) {
        if (hist1[i] == 0) error1[i] = sqrt(3.09);
        else error1[i] = sqrt(hist1[i]);
    }

    for (int i = 0; i < 100; i++) {
        if (hist2[i] == 0) error2[i] = sqrt(3.09);
        else error2[i] = sqrt(hist2[i]);
    }
    for (int i = 0; i < 100; i++) {
        std::cout << "error1 = \t" << error1[i] << "\terror2 = \t" << error2[i] << std::endl;
    }
    TMinuit *gMinuit = new TMinuit(4);
    gMinuit->SetFCN(fcn);

    double arglist[10];
    int ierflg = 0;

    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

    // Set starting values and step sizes for parameters
    static double vstart[4] = {1., 530., 10., 1.};
    static double step[4] = {0.1, 0.1 , 0.1 , 0.1};
    gMinuit->mnparm(0, "ampl", vstart[0], step[0], 0, 0, ierflg);
    gMinuit->mnparm(1, "mean gaus", vstart[1], step[1], 0, 0, ierflg);
    gMinuit->mnparm(2, "sigma gaus", vstart[2], step[2], 0, 0, ierflg);
    gMinuit->mnparm(3, "const", vstart[3], step[3], 0, 0, ierflg);

    // Now ready for minimization step
    arglist[0] = 500;
    arglist[1] = 0.1;
    gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

    // Print results
    double amin, edm, errdef;
    int nvpar, nparx, icstat;
    gMinuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
    gMinuit->mnprin(3, amin);

    //Инициализация функций подгонки
    double ampl, mean, sigma, Const;
    double ampl_err, mean_err, sigma_err, Const_err;
    gMinuit->GetParameter(0, ampl, ampl_err);
    gMinuit->GetParameter(1, mean, mean_err);
    gMinuit->GetParameter(2, sigma, sigma_err);
    gMinuit->GetParameter(3, Const, Const_err);

    auto f1 = new TF1("f1", "[0] * TMath::Gaus(x, [1], [2]) + [3]", 500, 600);
    f1->SetParameter(0, ampl);
    f1->SetParameter(1, mean);
    f1->SetParameter(2, sigma);
    f1->SetParameter(3, Const);

    auto f2 = new TF1("f2", "[0]", 500, 600);
    f2->SetParameter(0, Const);

    std::cout << "chi2 = \t" << chi2 << std::endl;

    TCanvas* c = new TCanvas("c", " ", 800, 600);
    c->Divide(1, 2);
    c->cd(1);
    h1->SetLineWidth(2);
    h1->Draw("e");
    f1->Draw("same");
    c->cd(2);
    h2->SetLineWidth(2);
    h2->Draw("e");
    f2->Draw("same");
    c->SaveAs("resluts_hists.png");    
}
