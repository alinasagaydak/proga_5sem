#include <TRandom.h>
#include <cmath>

#define _USE_MATH_DEFINES

void task4() {
    TRandom3* rand = new TRandom3(0);
    double L = 4; // расстояние между линиями
    double l = 3; // длина иголки
    double N = pow(10, 4); // сколько раз кинули иголку
    double n = 0; // сколько раз было пересечение

    for (int i = 0; i < N; i++) {
        double phi = rand->Rndm() * M_PI / 2.;
        double x = rand->Rndm() * l;
        if (x <= l * cos(phi)) n++;
    }
    
    double my_pi = 4 * N * l / (L * n);     
    std::cout << my_pi << std::endl;
}
