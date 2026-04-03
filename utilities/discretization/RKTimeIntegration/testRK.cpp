#include <iostream>
#include "rkTableau.H"

int main() {
    auto rk = rkTableau::CarpenterLSRK33();
    std::cout << "Stages: " << rk.stages() << std::endl;
    for (auto a : rk.alpha())
        std::cout << "alpha: " << a << std::endl;
    for (auto b : rk.beta())
        std::cout << "beta: " << b << std::endl;
    return 0;
}