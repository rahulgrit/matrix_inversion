
/** @file matinv.h 
*
*   @brief Gauss-Jordan method for matrix inversion using xtensor
*
*   matinv.h
*   @date 20-3-2018
*/
#ifndef MATINV
#define MATINV
#include <stdexcept>
#include "xtensor/xio.hpp"
#include "xtensor/xarray.hpp"
#include "xtensor/xfixed.hpp"
#include "xtensor/xrandom.hpp"
#include "xtensor/xsort.hpp"

/// structure template that inverts a given N x N matrix
//  Gauss-Jordan elimination with pivoting

template <unsigned N> struct matinv{
	template <typename I1TYPE>
	 static auto doIt(I1TYPE& F) noexcept{

	    //appending identity matrix

            xt::xtensorf<double, xt::xshape<N, N>> y = xt::eye<double>({N,N});
	    xt::xarray<double> b = xt::concatenate(xt::xtuple(F, y),1); 
    	 
            for (int k = 0; k < N; ++k){
    		//estimating position of largest pivot
                size_t p = k;
        	auto res0 = xt::view(b,xt::range(k,N),k);
        	xt::xarray<double> valmax = xt::amax(xt::fabs(res0));
        	xt::xarray<int> pos = xt::argmax(xt::fabs(res0));
        	p = pos(0) + k;
	
                // check if non-singular
		if(valmax(0) < 0) throw std::runtime_error("Matrix not invertible");
      		
                // exchange current row and pivot row
                if( p != k){
       	    	    auto swpa = xt::view(b, p, xt::range(0,2*N));
	    	    auto swpb = xt::view(b, k, xt::range(0,2*N));
	    	    swpa = swpa + swpb;
	    	    swpb = swpa - swpb;
	    	    swpa = swpa - swpb;
		}
                // divide pivot row by pivot value
	        double lval = b(k,k);
	        auto res1 = xt::view(b, k, xt::range(0,2*N));
	        res1 = res1 * (1/lval);
	
	        for (int i = 0; i < N; ++i){
	            if(i != k){
	            // subtract row k with elements
		    // multiplied by kth element of current row
		    double lval = b(i,k);
 	            auto rwslice = xt::view(b, i, xt::range(0,2*N));
	            rwslice = rwslice - lval*(res1);
	            }
	        }
            }
            // inverse matrix extracted from the appended matrix
	    xt::xtensorf<double, xt::xshape<N,N>> invmat = xt::view(b,xt::range(0,N),xt::range(N,2*N));
            return invmat;
	}
};
#endif
