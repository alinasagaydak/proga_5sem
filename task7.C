#include <TTree.h>
#include <TFile.h>

void task7() {
    TFile ifile("m3pimc.root", "read");
    TTree* h10 = (TTree*)ifile.Get("h10");
    int nph; Float_t eph;
    Float_t thetaph[nph];
    Float_t phiph[nph];
    unsigned char Isrfilter;
    float chi2_3p;

    h10->SetBranchAddress("nph", &nph);
    h10->SetBranchAddress("eph", &eph);
    h10->SetBranchAddress("thetaph", thetaph);
    h10->SetBranchAddress("phiph", phiph);
    h10->SetBranchAddress("Isrfilter", &Isrfilter);
    h10->SetBranchAddress("chi2_3p", &chi2_3p);


    TFile* ofile = new TFile("newroot.root", "RECREATE");
    TTree* MyTree = new TTree("MyTree", "tree with data");
    int my_nph; Float_t my_eph;
    Float_t my_thetaph[my_nph];
    Float_t my_phiph[my_nph];

    MyTree->Branch("nph", &my_nph);
    MyTree->Branch("eph", &my_eph);
    MyTree->Branch("thetaph", my_thetaph, "my_thetaph[my_nph]/F");
    MyTree->Branch("phiph", my_phiph, "my_phiph[my_nph]/F");

    for (int i = 0; i < h10->GetEntries(); i++) {
        h10->GetEntry(i);
        my_nph = nph;
        my_eph = eph;
        if (Isrfilter == 1 && chi2_3p < 30) {
            for (int j; j < nph; j++) {
                my_thetaph[j] = thetaph[j];
                my_phiph[j] = phiph[j];
            }
        } 
        MyTree->Fill();
    }

    ifile.Close();
    
    ofile->cd();
    MyTree->Write();
    ofile->Close();
}
