#include <vector>
#include <iostream>

void newnewtask8() {
//////////////////////////////////////////////////////////
//   This file has been automatically generated 
//     (Sun Nov 24 00:45:49 2024 by ROOT version6.32.06)
//   from TTree MyTree/tree with data
//   found on file: newroot.root
//////////////////////////////////////////////////////////


//Reset ROOT and connect tree file
   gROOT->Reset();
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("newroot.root");
   if (!f) {
      f = new TFile("newroot.root");
   }
   TTree* MyTree = (TTree*)f->Get("MyTree");

//Declaration of leaves types
   Int_t           nph;
   Float_t         eph[45];
   Float_t         thetaph[45];
   Float_t         phiph[45];

   // Set branch addresses.
   MyTree->SetBranchAddress("nph",&nph);
   MyTree->SetBranchAddress("eph",eph);
   MyTree->SetBranchAddress("thetaph",thetaph);
   MyTree->SetBranchAddress("phiph",phiph);

//     This is the loop skeleton
//       To read only selected branches, Insert statements like:
// MyTree->SetBranchStatus("*",0);  // disable all branches
// TTreePlayer->SetBranchStatus("branchname",1);  // activate branchname

    Long64_t nentries = MyTree->GetEntries();
    std::vector<int> all_candidates;
    TH1D* h_invm = new TH1D("h_invm","Inv Mass;Gev", 100, 0., 1.);
    TH1D* h_ang = new TH1D("h_ang", "Angle between photons;rad", 100, 0., TMath::Pi());

    Long64_t nbytes = 0;
    for (Long64_t i=0; i<nentries;i++) {
        nbytes += MyTree->GetEntry(i);
        int p0_candidate = 0;
        std::vector<double> inv_mass_vec;
        for (int j = 0; j < nph; j++) {
            for (int k = j; k < nph; k++) {
                double ang = TMath::ACos(sin(thetaph[j]) * sin(thetaph[k]) * cos(phiph[j] - phiph[k]) + cos(phiph[j]) * cos(phiph[k]));
                h_ang->Fill(ang);
                double m_inv = sqrt(2. * eph[j] * eph[k] * (1. - (sin(thetaph[j]) * sin(thetaph[k]) * cos(phiph[j] - phiph[k]) + cos(thetaph[j]) * cos(thetaph[k]))));
                inv_mass_vec.push_back(m_inv);
                if (0.1 <= m_inv && m_inv < 0.2) p0_candidate++;
            }
        }
        all_candidates.push_back(p0_candidate);
        if (p0_candidate == 2) for (double val : inv_mass_vec) h_invm->Fill(val);
    }
    
    TFile* ofile = new TFile("ofile.root", "RECREATE");
    h_invm->Write();
    h_ang->Write();
//    ofile->Close();
}
