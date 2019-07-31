/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   tFunc.cpp
 * Author: thuydt
 * 
 * Created on July 8, 2019, 3:31 PM
 */

#include "tFunc.h"
//#include <cmath>
//#include <iostream>
#include </opt/local/include/boost/math/special_functions/bessel.hpp> 

tFunc::tFunc() {
}

tFunc::tFunc(const tFunc& orig) {
}

tFunc::~tFunc() {
}

double tFunc::fL(double x){
    double I0 = boost::math::cyl_bessel_i(0, -0.5 * x);//cyl_bessel_i(0, x);
    double I1 = boost::math::cyl_bessel_i(1, -0.5 * x);
    double term2 = (1 - x) * I0 - x * I1;
    return exp(0.5 * x) * term2;
}
