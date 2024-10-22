// ROOT libraries
#include <TRandom3.h>
#include <TH1.h>
#include <TF1.h>
#include <TFile.h>
#include <TString.h>

// C++ libraries
#include <cmath>
#include <iostream>
#include <vector>
#include <thread>

Double_t inv_func(Double_t x){ return sqrt(-log(1-x)); }

Double_t myfunc(Double_t x) { return exp(-x) * pow(x, 3); }

Double_t AverageValMethod(Double_t (*myfunc)(Double_t x), Double_t N_sample) {
    TRandom3* rand = new TRandom3(0);
    Double_t sum = 0;
    for (size_t i = 0; i < N_sample; i++) {
        Double_t r = rand->Rndm();
        sum += (*myfunc)(r);
    }
    return 1 / N_sample * sum;
}

Double_t GeometricAlgorithm(Double_t (*myfunc)(Double_t x), Double_t N_sample) {
    TRandom3* rand = new TRandom3(0);
    Double_t M = (*myfunc)(1);
    Double_t N_selected = 0.;
    for (size_t i = 0; i < N_sample; i++) {
        Double_t eta = rand->Rndm(); // x
        Double_t xi = M * rand->Rndm(); // y
        if (xi <=  (*myfunc)(eta)) N_selected++;
    }
    return (N_selected / N_sample) * M;
}

Double_t MainPartIntegration(Double_t N_sample) {
    //Double_t N_selected = 0.;
    TRandom3* rand = new TRandom3(0);
    Double_t sum = 0.;
    for (size_t i = 0; i < N_sample; i++) {
        Double_t r = rand->Rndm();
        Double_t xi = pow(r, 0.25);
        //Double_t x = gRandom->Rndm();
        //Double_t y = gRandom->Rndm();
        //if (xi < exp(r)) N_selected += 1;
        sum += exp(xi)/4;
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
    Double_t N_sample = pow(10, 5);
    Int_t num_itrs = 1000;
    std::vector<Double_t> resultsVec(num_itrs * 3);

    // Сразу запускаем потоки, чтобы программа быстрее завершилась. Три вида интегрирования, набираем данные    
    std::thread th1([&resultsVec, &N_sample, &num_itrs]()
    {
      for (int i = 0; i < num_itrs; i++) {
        Double_t res = AverageValMethod((*myfunc), N_sample);
        resultsVec[i] = res;  
      }             
    });

    std::thread th2([&resultsVec, &N_sample, &num_itrs]()
    {
        for (size_t i = num_itrs; i < 2 * num_itrs; i++) {
            Double_t res = GeometricAlgorithm((*myfunc), N_sample);
            resultsVec[i] = res;
        }     
    });

    std::thread th3([&resultsVec, &N_sample, &num_itrs]()
    {
        for (size_t i = 2 * num_itrs; i < 3 * num_itrs; i++) {
            Double_t res = MainPartIntegration(N_sample);
            resultsVec[i] = res;
        }        
    });

    // метод обратных функций f(x) = 2*x*exp(-x^2)
    TH1D* h_invfunc = new TH1D("h_invfunc", "The method of inverse functions", 150, 0., 3.5);
    for (size_t i = 0; i < N_sample; i++) {
         Double_t r = gRandom->Rndm();
         Double_t xi = inv_func(r);
         h_invfunc->Fill(xi);
    }

    // вычисление интеграла от 0 до 1 функции x^3 * exp(-x)
    Double_t resAverVal = AverageValMethod((*myfunc), N_sample);
    Double_t resGeomAlg = GeometricAlgorithm((*myfunc), N_sample);
    Double_t resMainPartInt = MainPartIntegration(N_sample);

    std::vector<Double_t> errors = MethodsError(N_sample);
    PrintResults(resAverVal, resGeomAlg, resMainPartInt, errors);


    th1.join();
    th2.join();
    th3.join();

    TH1D* h_resAverVal = new TH1D("h_resAverVal", Form("Average Value Method Results, N_sample = %0.f", N_sample), 100, 0.11, 0.12);
    std::for_each(resultsVec.begin(), resultsVec.begin() + num_itrs, [&h_resAverVal](Double_t& res){h_resAverVal->Fill(res);});

    TH1D* h_resGeomAlg = new TH1D("h_resGeomAlg", Form("Geometric Algorithm Results, N_sample = %0.f", N_sample), 100, 0.11, 0.12);
    std::for_each(resultsVec.begin() + num_itrs, resultsVec.begin() + 2*num_itrs, [&h_resGeomAlg](Double_t& res){h_resGeomAlg->Fill(res);});
    
    TH1D* h_resMainPartInt = new TH1D("h_resMainPartInt", Form("Main Part Integration Method Results, N_sample = %0.f", N_sample), 100, 0.56, 0.57);
    std::for_each(resultsVec.begin() + 2*num_itrs, resultsVec.end(), [&h_resMainPartInt](Double_t& res){h_resMainPartInt->Fill(res);});


    TFile* ofile = new TFile("ofile.root", "RECREATE");
    h_invfunc->Write();
    h_resAverVal->Write();
    h_resGeomAlg->Write();
    h_resMainPartInt->Write();
}
