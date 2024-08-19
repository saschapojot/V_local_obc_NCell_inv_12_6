//
// Created by polya on 8/16/24.
//

#ifndef V_LOCAL_OBC_NCELL_INV_12_6_POTENTIALFUNCTIONPROTOTYPE_HPP
#define V_LOCAL_OBC_NCELL_INV_12_6_POTENTIALFUNCTIONPROTOTYPE_HPP
#include <cstring> // For memcpy
#include <fstream>
#include <immintrin.h>
#include <iostream>
#include <math.h>
#include <memory>
#include <regex>
#include <stdexcept>
#include <string>


class potentialFunction {
//base class for potential function
public:
    virtual double operator()(const double * xVec, const int& j) = 0;
//    virtual double plain_for(const double&L,const double *d0Vec, const double *d1Vec)=0;
    virtual void json2Coefs(const std::string &coefsStr)=0;
    virtual  void init()=0;
    virtual double potentialFull(const double * xVec)=0;
    virtual double getLm() const = 0; //  method to get lm
    virtual double get_eps() const = 0; //  method to get eps
    virtual ~ potentialFunction() {};
};

std::shared_ptr<potentialFunction>  createPotentialFunction(const std::string& funcName, const std::string &row) ;


#endif //V_LOCAL_OBC_NCELL_INV_12_6_POTENTIALFUNCTIONPROTOTYPE_HPP
