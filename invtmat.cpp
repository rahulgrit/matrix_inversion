/** @file invtmat.cpp 
*
*
*   @date 21-3-2018
*/
#include <iostream>
#include "matinv.h"
#include "xtensor/xrandom.hpp"
#include <chrono>

int main(){

    const int N = 5;
    xt::xtensorf<double, xt::xshape<N,N>> rdn = xt::random::rand<double>({N,N},-5,5);

////  to test different number of matrices    
//    auto rdn = xt::random::rand<double>({N,N},-5,5);
//    double time_c=0;
//    for(int k = 0; k < 1e3; ++k){
//    const auto t0 = std::chrono::high_resolution_clock::now();
    std::cout << "Input Matrix :\n" << rdn << std::endl;
    xt::xarray<double> op = matinv<N>::doIt(rdn);
//    const auto t1 = std::chrono::high_resolution_clock::now();
//    const std::chrono::duration<double> dt = t1 - t0;
//    time_c += dt.count();}
    std::cout << "Output Matrix :\n" << op << std::endl;

return 0;}
