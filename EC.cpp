/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   EC.cpp
 * Author: thuydt
 * 
 * Created on July 30, 2019, 4:36 PM
 */

#include "EC.h"
#include <math.h>
#include <iostream>
#include "./lib/declare.h"

//empirical cost
using namespace std;
EC::EC(const EC& orig) {
}

EC::~EC() {
}
EC::EC(int n, int m, std::vector<tPair>* edge,  double al, double ld, std::vector<double>* xp, 
        std::vector<double>* yp, double* uu, double *vv, int* zz){
    N = n;
    M = m;
    E = edge;
    alp = al;
    lambda = ld;
    X = xp;
    Y = yp;
    u = uu;
    v = vv;
    z = zz;
    
    
}

double EC::getEC(){
    //std::cout <<"\nin EC::getEC()";
    int count = 0;
    
    double exp[N];
    
    for (int i = 0; i < N; i++){
        exp[i] = 0;
        int num_exp = X[i].size();
        int m = z[i];
        for (int j = 0; j < num_exp; j++){
            double dist = pow(X[i][j] - u[m], 2) + pow(Y[i][j] - v[m], 2);
            
            exp[i] += pow(dist, alp/2.0);
        }
        exp[i] = exp[i] / num_exp;

    }
    
    double term1 = 0;
    for (int i = 0; i < N; i++){
        term1 += exp[i];
    }
    term1 = lambda * term1;
    
    double term2 = 0;
    for (int i = 0; i < N; i++){
        for (tPair f: E[i]){
            int j = f.key;
            if (z[i] != z[j]){ // i and j are assigned to different events
                term2 += f.val;
            }
        }
    }
    term2 = term2 * (1-lambda);
    return term1 + term2;
    
}

