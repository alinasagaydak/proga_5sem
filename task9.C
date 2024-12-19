#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TGraphPainter.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TDirectory.h>
#include <iostream>

double fitfunc_th(double* x, double* par) {
    return par[0] * TMath::Gaus(x[0], par[1], par[2]) + par[3] * TMath::Gaus(x[0], par[4], par[5]);
}

double fitfunc_phi(double* x, double* par) {
    return par[0];
}

void task9() {
    TFile file("newroot.root","read"); // открываем файл
    TTree* tree = (TTree*)file.Get("MyTree"); // достаем дерево
    int nph;
    float eph[45];
    float thetaph[45];
    float phiph[45];

    tree->SetBranchAddress("nph", &nph);
    tree->SetBranchAddress("eph", eph); 
    tree->SetBranchAddress("thetaph", thetaph); 
    tree->SetBranchAddress("phiph", phiph);

    const double critical_phi[] = {0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180.,
                          190., 200., 210., 220., 230., 240., 250., 260., 270., 280., 290., 300., 310., 320., 330., 340., 350., 360.};
    const double critical_theta[] = {0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180.};

    const int size_crtical_phi = sizeof(critical_phi) / sizeof(double);
    const int size_crtical_theta = sizeof(critical_theta) / sizeof(double);

    double candidates_phi[size_crtical_phi] = {0};
    double candidates_theta[size_crtical_theta] = {0};


    for (size_t i = 0, N = tree->GetEntries(); i < N; i++) {
        tree->GetEntry(i); // читаем событие
        for (int j = 0; j < nph - 1; j++) {
            for (int k = j + 1; k < nph; k++) {
                TLorentzVector vecph1(eph[j], 0., 0., eph[j]);
                vecph1.SetTheta(thetaph[j]);
                vecph1.SetPhi(phiph[j]);

                TLorentzVector vecph2(eph[k], 0., 0., eph[k]);
                vecph2.SetTheta(thetaph[k]);
                vecph2.SetPhi(phiph[k]);

                TLorentzVector vecpi = vecph1 + vecph2;
                double m_inv = vecpi.M();
                if (0.1 <= m_inv && m_inv < 0.2) {
                    double pitheta = vecpi.Theta() * 180. / TMath::Pi();
                    double piphi = (vecpi.Phi() + TMath::Pi()) * 180. / TMath::Pi();
                
                    for (int phi = 0; phi < size_crtical_phi; phi++) {
                        if (piphi > critical_phi[phi]) continue;
                        else {
                            candidates_phi[phi-1] += 1;
                            break;
                        }
                    }
                    for (int th = 0; th < size_crtical_theta; th++) {
                        if (pitheta > critical_theta[th]) continue;
                        else {
                            candidates_theta[th-1] += 1;
                            break;
                        }
                    }
                }
            }
        }
    }

// Azimuthal
    TF1* fitphi = new TF1("fitphi", fitfunc_phi, 0., 360., 1);
    fitphi->SetLineColor(38);
    double err_phi[size_crtical_phi] = {0.};
    double err_bins_phi[size_crtical_phi] = {0.};
    for (int i = 0; i < size_crtical_phi; i++) err_phi[i] = sqrt(candidates_phi[i]);
    
    auto gr_phi = new TGraphErrors(size_crtical_phi, critical_phi, candidates_phi, err_bins_phi, err_phi);
    gr_phi->SetMarkerColor(4);
    gr_phi->SetMarkerStyle(21);
    gr_phi->Fit(fitphi, "R");


// Polar
    TF1* fitth = new TF1("fitth", fitfunc_th, 0., 180., 6);
    fitth->SetLineColor(46);
    double err_theta[size_crtical_theta] = {0.};
    double err_bins_th[size_crtical_theta] = {0.};
    for (int i = 0; i < size_crtical_theta; i++) err_theta[i] = sqrt(candidates_theta[i]);
    
    auto gr_th = new TGraphErrors(size_crtical_theta, critical_theta, candidates_theta, err_bins_th, err_theta);
    gr_th->SetMarkerColor(2);
    gr_th->SetMarkerStyle(21);
    gr_th->Fit(fitth, "R");
    //
    
// Multigraph    
    auto c1 = new TCanvas("c1"," ");
    c1->SetGrid();
    auto multigr = new TMultiGraph();
    multigr->SetTitle("Distribution of the number of #pi^{0} candidates by angles; angle, degree; number of candidates");
    multigr->Add(gr_phi);
    multigr->Add(gr_th);
    multigr->Draw("AP");
    
    auto* legend = new TLegend(0.7, 0.7, .85, .8);
    legend->AddEntry(gr_phi, "By azimuthal angle", "1ep");
    legend->AddEntry(gr_th, "By polar angle", "1ep");
    legend->AddEntry(fitphi, "fit function", "l");
    legend->AddEntry(fitth, "fit function", "l");
    legend->Draw();

    file.Close();
// Write to dir
    TFile newfile("newroot.root","update");
    TDirectory* newsubdir = newfile.mkdir("newsubdir","new_dir_task9");
    newsubdir->cd();
    c1->Write();
}
