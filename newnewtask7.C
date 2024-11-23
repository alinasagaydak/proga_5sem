#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <vector>
#include <cmath>

double fitf(double* x, double* par) {
   return par[0] * pow(x[0], par[1]) * exp(-par[2] * x[0]) + par[3] * pow(x[0], par[1]) * exp(-par[4] * x[0]);
}

void newnewtask7() {
    TFile ifile("m3pimc.root", "read");
    TTree* h10 = (TTree*)ifile.Get("h10");
    int nph;
    float eph[45];
    float thetaph[45];
    float phiph[45];
    unsigned char Isrfilter;
    float chi2_3p;

    h10->SetBranchAddress("nph", &nph);
    h10->SetBranchAddress("eph", eph);
    h10->SetBranchAddress("thetaph", thetaph);
    h10->SetBranchAddress("phiph", phiph);
    h10->SetBranchAddress("Isrfilter", &Isrfilter);
    h10->SetBranchAddress("chi2_3p", &chi2_3p);

    TFile* ofile = new TFile("newroot.root", "RECREATE");
    TTree* MyTree = new TTree("MyTree", "tree with data");
    int my_nph;
    std::vector<float> my_eph;
    std::vector<float> my_thetaph;
    std::vector<float> my_phiph;
  
    MyTree->Branch("nph", &my_nph);
    MyTree->Branch("eph", &my_eph);
    MyTree->Branch("thetaph", &my_thetaph);
    MyTree->Branch("phiph", &my_phiph);
    
    TH1D* h_eph = new TH1D("h_eph", "eph;Mev;Events", 50, 0., 9.);

    for (Long64_t i = 0; i < h10->GetEntries(); i++) {
        h10->GetEntry(i);
        if (nph == 0) continue;
        else if (Isrfilter == 1 & chi2_3p < 30) {
            my_nph = nph;
            for (int j = 0; j < nph; j++) {
                my_eph.push_back(eph[j]);
                my_thetaph.push_back(thetaph[j]);
                my_phiph.push_back(phiph[j]);

                h_eph->Fill(eph[j]);
            }
            MyTree->Fill();
            
            my_eph.erase(my_eph.begin(), my_eph.end());
            my_thetaph.erase(my_thetaph.begin(), my_thetaph.end());
            my_phiph.erase(my_phiph.begin(), my_phiph.end());
        } 
    }

    TF1* fitfunc = new TF1("fitfunc", fitf, 0, 9, 5);
    fitfunc->SetParameters(1e4, -1, 0.5, 1e3, 0.1);
    fitfunc->SetParNames("A1", "n", "k1", "A2", "k2");
    fitfunc->SetParLimits(1, -5, 0);
    fitfunc->SetParLimits(2, 0, 10);
    fitfunc->SetParLimits(4, 0, 10);

    auto c = new TCanvas("c", "c");
    c->cd();
    c->SetLogy();
    h_eph->Fit(fitfunc, "R");
    h_eph->Draw();
    c->Draw();


    ofile->cd();
    MyTree->Write();
    c->Write();
    ofile->Close();
    ifile.Close();
}
