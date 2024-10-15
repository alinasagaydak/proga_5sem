// ROOT libraries
#include <TRandom3.h>
#include <TH1.h>
#include <TF1.h>
#include <TFile.h>

// C++ libraries
#include <cmath>
#include <iostream>
#include <vector>

Double_t inv_func(Double_t x){ return sqrt(-log(1-x)); }

Double_t myfunc(Double_t x) { return exp(-x) * pow(x, 3); }

Double_t AverageValMethod(Double_t (*myfunc)(Double_t x), Double_t N_sample) {
    Double_t sum = 0;
    for (size_t i = 0; i < N_sample; i++) {
        Double_t r = gRandom->Rndm();
        sum += (*myfunc)(r);
    }
    return 1 / N_sample * sum;
}

Double_t GeometricAlgorithm(Double_t (*myfunc)(Double_t x), Double_t N_sample) {
    Double_t M = (*myfunc)(1);
    Double_t N_selected = 0.;
    for (size_t i = 0; i < N_sample; i++) {
        Double_t eta = gRandom->Rndm(); // x
        Double_t xi = M * gRandom->Rndm(); // y
        if (abs(xi - (*myfunc)(eta)) <= pow(10, -3)) N_selected += 1.;
    }
    return (N_selected / N_sample) * M;
}

Double_t MainPartIntegration(Double_t N_sample) {
    //Double_t N_selected = 0.;
    Double_t sum = 0.;
    for (size_t i = 0; i < N_sample; i++) {
        Double_t r = gRandom->Rndm();
        Double_t xi = pow(r, 0.25);
        //Double_t x = gRandom->Rndm();
        //Double_t y = gRandom->Rndm();
        //if (xi < exp(r)) N_selected += 1;
        sum += exp(xi);
    }
    return 1 / N_sample * sum;
    //return N_selected / N_sample;
}

std::vector<Double_t> MethodsError(Double_t N_sample) {
    Double_t I = 6 - 16 / std::exp(1.0);
    Double_t In_AverVal = AverageValMethod((*myfunc), N_sample);
    Double_t In_GeomAlg = GeometricAlgorithm((*myfunc), N_sample);
    Double_t In_MainPartInt = MainPartIntegration(N_sample);
    std::vector<Double_t> errors = {abs(In_AverVal - I), abs(In_GeomAlg - I), abs(In_MainPartInt - I)};
    return errors;
}

void PrintResults(Double_t& resAverVal, Double_t& resGeomAlg, Double_t& resMainPartInt, std::vector<Double_t>& errors) {
    std::cout << "Average Value Method:\t" << resAverVal << "\terror:\t" << errors[0] << std::endl;
    std::cout << "Geometric Algorithm:\t" << resGeomAlg << "\terror:\t" << errors[1] << std::endl;
    std::cout << "Main Part Integration:\t" << resMainPartInt << "\terror:\t" << errors[2] << std::endl;
}

void task5() {
    // метод обратных функций f(x) = 2*x*exp(-x^2)
    TH1D* h_invfunc = new TH1D("h_invfunc", "The method of inverse functions", 150, 0., 3.5);

    for (size_t i = 0; i < 1000000; i++) {
        Double_t r = gRandom->Rndm();
        Double_t xi = inv_func(r);
        h_invfunc->Fill(xi);
    }
    
    // вычисление интеграла от 0 до 1 функции x^3 * exp(-x)
    Double_t resAverVal = AverageValMethod((*myfunc), 1000000.);
    Double_t resGeomAlg = GeometricAlgorithm((*myfunc), 1000000.);
    Double_t resMainPartInt = MainPartIntegration(1000000.);

    std::vector<Double_t> errors = MethodsError(1000000.);

    PrintResults(resAverVal, resGeomAlg, resMainPartInt, errors);

    TFile* ofile = new TFile("ofile.root", "RECREATE");
    h_invfunc->Write();    
}
