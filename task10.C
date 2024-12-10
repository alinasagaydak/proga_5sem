#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>

#include <TH1D.h>
#include <TCanvas.h>
#include <TMinuit.h>
#include <TFitter.h>
#include <TVirtualFitter.h>

//Global variables
double hist1[100];
double hist2[100];
double error1[100];
double error2[100];
double x[100];
double ampl;


double func(double ampl, double x, double *par) {
    // ampl - кол-во сигнала
    return ampl * TMath::Gaus(x, par[0], par[1]);
}

void fcn(int &npar, double *gin, double &f, double *par, int iflag) {
    const int nbins = 100;

//calculate chisquare
    double chisq = 0;
    for (int i = 0; i < nbins; i++) {
        double delta1  = (hist1[i] - func(ampl, x[i], par) - par[2]) / error1[i];
        double delta2 = (hist2[i] - par[2]) / error2[i];
        chisq += delta1 * delta1 + delta2 * delta2;
    }
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

    std::vector<double> error1(data1.size());
    std::vector<double> error2(data2.size());
    std::for_each(error1.begin(), error1.end(), [&data1](double err){for (double val : data1) err = sqrt(val);});
    std::for_each(error2.begin(), error2.end(), [&data2](double err){for (double val : data2) err = sqrt(val);});

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
    
    ampl = 0.;
    for (double val : hist1) ampl += val;
    
    for (int i = 0; i < 100; i++) x[i] = i + 500;
    
    //The errors values
    //for (int i = 0; i < 100; i++) {
    //    if (hist1[i] != 0) { error1[i] = sqrt(hist1[i]); continue; }
    //    if (hist2[i] != 0) { error2[i] = sqrt(hist2[i]); continue; }
    //    error1[i] = sqrt(3.09); // уровень достоверности 95%
    //    error2[i] = sqrt(3.09); 
   // }
    for (int i = 0; i < 100; i++) {
        if (hist1[i] == 0) error1[i] = sqrt(3.09);
        else error1[i] = sqrt(hist1[i]);
    }

    for (int i = 0; i < 100; i++) {
        if (hist2[i] == 0) error2[i] = sqrt(3.09);
        else error2[i] = sqrt(hist2[i]);
    }

    TMinuit *gMinuit = new TMinuit(3);
    gMinuit->SetFCN(fcn);

    double arglist[10];
    int ierflg = 0;

    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

    // Set starting values and step sizes for parameters
    static double vstart[3] = {550., 20., 1.};
    static double step[3] = {0.1 , 0.1 , 0.1};
    gMinuit->mnparm(0, "mean gaus", vstart[0], step[0], 0, 0, ierflg);
    gMinuit->mnparm(1, "sigma gaus", vstart[1], step[1], 0, 0, ierflg);
    gMinuit->mnparm(2, "const", vstart[2], step[2], 0, 0, ierflg);

    // Now ready for minimization step
    arglist[0] = 500;
    arglist[1] = 0.1;
    gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

    // Print results
    //double amin, edm, errdef;
    //int nvpar, nparx, icstat;
    //gMinuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
    //gMinuit->mnprin(3, amin);


    TCanvas* c = new TCanvas();
    c->Divide(1, 2);
    c->cd(1);
    h1->SetLineWidth(2);
    h1->Draw("e");
    c->cd(2);
    h2->SetLineWidth(2);
    h2->Draw("e");
        
}
