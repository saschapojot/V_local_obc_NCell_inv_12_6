//
// Created by polya on 8/17/24.
//


#include "mc_read_load_compute.hpp"



///
/// @param x
/// @param leftEnd
/// @param rightEnd
/// @param eps
/// @return return a value within distance eps from x, on the open interval (leftEnd, rightEnd)
double mc_computation::generate_uni_open_interval(const double &x, const double &leftEnd, const double &rightEnd, const double &eps){


    double xMinusEps=x-eps;
    double xPlusEps=x+eps;

    double unif_left_end=xMinusEps<leftEnd?leftEnd:xMinusEps;
    double unif_right_end=xPlusEps>rightEnd?rightEnd:xPlusEps;

//    std::random_device rd;
//    std::ranlux24_base e2(rd());

    double unif_left_end_double_on_the_right=std::nextafter(unif_left_end, std::numeric_limits<double>::infinity());



    std::uniform_real_distribution<> distUnif(unif_left_end_double_on_the_right,unif_right_end); //[unif_left_end_double_on_the_right, unif_right_end)

    double xNext=distUnif(e2);
    return xNext;



}

///
/// @param x proposed value
/// @param y current value
/// @param a left end of interval
/// @param b right end of interval
/// @param epsilon half length
/// @return proposal probability S(x|y)
double mc_computation::S_uni(const double &x, const double &y,const double &a, const double &b, const double &epsilon){

    if (a<y and y<a+epsilon){
        return 1.0/(y-a+epsilon);
    } else if( a+epsilon<=y and y<b+epsilon){
        return 1.0/(2.0*epsilon);
    }else if(b-epsilon<=y and y<b){
        return 1.0/(b-y+epsilon);
    } else{

        std::cerr<<"value out of range."<<std::endl;
        std::exit(10);


    }


}


///
/// @param xVecCurr
/// @param j
/// @param xVecNext
void mc_computation::proposal_uni(const std::shared_ptr<double[]> & xVecCurr,const int&j,std::shared_ptr<double[]>&xVecNext){
    double lm=potFuncPtr->getLm();


}
