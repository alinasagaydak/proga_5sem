#include <iostream>
#include <cmath>
#include <TRandom3.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TF1.h>

void task6() {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1); 

    TRandom3* rand = new TRandom3(0);

    double E_CM = 1020.0; // Полная энергия в CЦМ МэВ
    double m_Ks = 497.611; // Масса Ks M'D
    double m_pi = 139.570; // Масса pi МэВ
    double c = 29.9792458; // Скорость света в см/нс
    double tau_Ks = 0.0895; // Время жизни Ks в нс

    double rad_det = 30.0; // см
    double len_det = 50.0; // см
    double min_pion_momentum = 40.0; // MэВ/c

    auto h_theta_Ks_lab = new TH1D("h_theta_Ks_lab", "#theta_{Ks} in lab ref;angle, rad; num of events", 100, 0, TMath::Pi());
    auto h_phi_Ks_lab = new TH1D("h_phi_Ks_lab", "#phi_{Ks} in lab ref;angle, rad; num of events", 100, -TMath::Pi(), TMath::Pi());
    auto h_theta_pi_lab = new TH1D("h_theta_pi_lab", "#theta_{#pi} in lab ref;angle, rad; num of events", 100, 0, TMath::Pi());
    auto h_phi_pi_lab = new TH1D("h_phi_pi_lab", "#phi_{#pi} in lab ref;angle, rad; num of events", 100, -TMath::Pi(), TMath::Pi());
    auto h_theta_pi_Ks = new TH1D("h_theta_pi_Ks", "#theta_{#pi} in Ks ref;angle, rad; num of events", 100, 0, TMath::Pi());
    auto h_phi_pi_Ks = new TH1D("h_phi_pi_Ks", "#phi_{#pi} in Ks ref;angle, rad; num of events", 100, -TMath::Pi(), TMath::Pi());

    int N_total = pow(10, 5);
    int N_registered = 0;

    TF1* theta_dist_Ks = new TF1("theta_dist_Ks", "[0] * pow(sin(x), 3)", 0, TMath::Pi());
    TF1* theta_dist_pi = new TF1("theta_dist_pi", "[0] * sin(x)", 0, TMath::Pi());

    theta_dist_Ks->SetParameter(0, 1.);
    theta_dist_pi->SetParameter(0, 1.);
    
    for (Int_t i = 0; i < N_total; ++i) {
        double theta_Ks = theta_dist_Ks->GetRandom();
        double phi_Ks = 2 * TMath::Pi() * rand->Rndm() - TMath::Pi(); // От -pi до pi
        h_theta_Ks_lab->Fill(theta_Ks);
        h_phi_Ks_lab->Fill(phi_Ks);

        double theta_pi = theta_dist_pi->GetRandom();
        double phi_pi = 2 * TMath::Pi() * rand->Rndm() - TMath::Pi();
        h_theta_pi_Ks->Fill(theta_pi);
        h_phi_pi_Ks->Fill(phi_pi);
        
        double E_Ks = E_CM / 2.0;
        double p_Ks = sqrt(E_Ks * E_Ks - m_Ks * m_Ks);

        // Вектор импульса K_S
        double p_Ksx = p_Ks * sin(theta_Ks) * cos(phi_Ks);
        double p_Ksy = p_Ks * sin(theta_Ks) * sin(phi_Ks);
        double p_Ksz = p_Ks * cos(theta_Ks);
        TLorentzVector Ks_lab(p_Ksx, p_Ksy, p_Ksz, E_Ks);

        double E_pi = m_Ks / 2.0; 
        double p_pi = sqrt(E_pi * E_pi - m_pi * m_pi); 

        // Вектора импульсов пионов в системе K_S
        double p_pion_x = p_pi * sin(theta_pi) * cos(phi_pi);
        double p_pion_y = p_pi * sin(theta_pi) * sin(phi_pi);
        double p_pion_z = p_pi * cos(theta_pi);
        TLorentzVector pi_plus(p_pion_x, p_pion_y, p_pion_z, E_pi);
        TLorentzVector pi_minus(-p_pion_x, -p_pion_y, -p_pion_z, E_pi);

        // Направление движения K_S в лабораторной системе
        TVector3 Ks_direction = Ks_lab.Vect().Unit();

        double old_angle = pi_plus.Z() / pi_plus.P();
        // Поворачиваем импульсы пионов так, чтобы ось z совпала с направлением Ks
        pi_plus.RotateUz(Ks_direction);
        pi_minus.RotateUz(Ks_direction);
        double new_angle = pi_plus.Vect().Dot(Ks_direction) / pi_plus.P();

        // Скорость K_S в лабораторной системе
        TVector3 beta_Ks = Ks_lab.BoostVector();

        // Выполняем буст пионов в лабораторную систему
        pi_plus.Boost(beta_Ks);
        pi_minus.Boost(beta_Ks);

        // Заполняем гистограммы углов для пионов в лабораторной системе
        h_theta_pi_lab->Fill(pi_plus.Theta());
        h_theta_pi_lab->Fill(pi_minus.Theta());
        h_phi_pi_lab->Fill(pi_plus.Phi());
        h_phi_pi_lab->Fill(pi_minus.Phi());

        double beta_Ks_mag = beta_Ks.Mag();
        double gamma_Ks = Ks_lab.Gamma();

        // Средняя длина пробега K_S в см
        double lambda = beta_Ks_mag * gamma_Ks * c * tau_Ks;

        // Генерируем длину пробега K_S
        double r_decay = rand->Rndm();
        double L_decay = -lambda * log(r_decay);

        // Положение распада K_S
        TVector3 decay_vertex = Ks_direction * L_decay;


        // Траектории пионов
        TVector3 pi_plus_dir = pi_plus.Vect().Unit();
        TVector3 pi_minus_dir = pi_minus.Vect().Unit();

        // Расстояние до пересечения с цилиндрической поверхностью детектора
        double t_plus = (rad_det - decay_vertex.Perp()) / (pi_plus_dir.Perp());
        double t_minus = (rad_det - decay_vertex.Perp()) / (pi_minus_dir.Perp());

        // Положение пересечения с цилиндром
        TVector3 pi_plus_pos = decay_vertex + pi_plus_dir * t_plus;
        TVector3 pi_minus_pos = decay_vertex + pi_minus_dir * t_minus;

        // Проверяем, находятся ли точки пересечения внутри детектора по оси z
        bool registered_plus = (fabs(pi_plus_pos.Z()) <= len_det / 2.0) && (pi_plus.P() >= 40);
        bool registered_minus = (fabs(pi_minus_pos.Z()) <= len_det / 2.0) && (pi_minus.P() >= 40);

        if (registered_plus && registered_minus) {
            N_registered++;
        }
    }

    double efficiency = (double)N_registered / (double)N_total;
    double error = sqrt(efficiency * (1. - efficiency) / N_total);

    std::cout << "Total events: " << N_total << std::endl;
    std::cout << "Registered events: " << N_registered << std::endl;
    std::cout << "Reconstruction efficiency: " << efficiency * 100 << " ± " << error * 100 << " %" << std::endl;


    TCanvas* c1 = new TCanvas("c1");
    c1->Divide(2, 3);
    
    c1->cd(1);
    theta_dist_Ks->SetParameter(0, h_theta_Ks_lab->GetMaximum());
    h_theta_Ks_lab->Fit(theta_dist_Ks);
    h_theta_Ks_lab->SetLineWidth(2);
    h_theta_Ks_lab->Draw();

    c1->cd(2);
    h_phi_Ks_lab->SetLineWidth(2);
    h_phi_Ks_lab->Fit("pol0");
    h_phi_Ks_lab->Draw();

    c1->cd(3);
    theta_dist_pi->SetParameter(0, h_theta_pi_lab->GetMaximum());
    h_theta_pi_lab->Fit(theta_dist_pi);
    h_theta_pi_lab->SetLineWidth(2);
    h_theta_pi_lab->Draw();

    c1->cd(4);
    h_phi_pi_lab->SetLineWidth(2);
    h_phi_pi_lab->Fit("pol0");
    h_phi_pi_lab->Draw();

    c1->cd(5);
    theta_dist_pi->SetParameter(0, h_theta_pi_Ks->GetMaximum());
    h_theta_pi_Ks->Fit(theta_dist_pi);
    h_theta_pi_Ks->SetLineWidth(2);
    h_theta_pi_Ks->Draw();

    c1->cd(6);
    h_phi_pi_Ks->SetLineWidth(2);
    h_phi_pi_Ks->Fit("pol0");
    h_phi_pi_Ks->Draw();

}
