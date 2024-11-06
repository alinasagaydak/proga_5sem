#include <iostream>
#include <math.h>
#include <TF1.h>
#include <TMath.h>
#include <TCanvas.h>

//Определим функцию для в.ф.:
double psi(double* x, double* par) {
    // return exp(-fabs(x[0])/par[0]);
    return pow(2/TMath::Pi(), 0.25) * (1 / sqrt(par[0])) * exp(-pow(x[0]/par[0], 2));
}

double potential(double* x, double* par) {
    if (fabs(x[0]) < 10) return -0.5;
    else return 0;
}

//Определим подынтегральное выражение скалярного произведения:
double scalar_expr(double* x, double* par) {
    double A = par[0];
    TF1 Psi("psi", psi, -10. * A, 10. * A, 1);
    TF1 U("U", potential, -10. * A, 10. * A, 0);
    Psi.SetParameter(0, A);
    // psi* (-ihd/dx)(-ihd/dx) psi 
    // постоянная Планка Math::Hbarcgs()
    double m = 9.1 * pow(10, -28); // масса электронав граммах
    double result = - pow(TMath::Hbarcgs(), 2) / (2 * m) * Psi.Eval(x[0]) * Psi.Derivative2(x[0]) + Psi.Eval(x[0]) * Psi.Derivative2(x[0]) * U.Eval(x[0]);
    //double result = - pow(TMath::Hbarcgs(), 2) / (2 * m) * Psi.Eval(x[0]) * Psi.Derivative2(x[0]);
    return result;
}

//Скалярное произведение как функция от a:
double scal_prod(double* x,double* par) {
    TF1 scal_exp("sc", scalar_expr, -10. * x[0], 10. * x[0], 1);
    scal_exp.SetParameter(0, x[0]);
    return scal_exp.Integral(-10. * x[0], 10. * x[0]);
}

void scalar_prod(double amin, double amax) {
    gStyle->SetLabelSize(0.05, "xy");
    gStyle->SetTitleSize(0.05, "xy");
    TCanvas *s = new TCanvas();
    TF1 *scalar_pr = new TF1("sp", scal_prod, amin, amax, 0);
    scalar_pr->Draw();  
    scalar_pr->SetTitle("<p2> as a function of a");
    scalar_pr->GetXaxis()->SetTitle("a");
    scalar_pr->GetYaxis()->SetTitle("<p2>");
    s->SaveAs("sc_pr.png");
}

void task3() {
    double amin = 0;
    double amax = 1;

    scalar_prod(amin, amax);
}

