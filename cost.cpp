/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   cost.cpp
 * Author: thuydt
 * 
 * Created on February 19, 2019, 4:17 PM
 */
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <stdexcept>
#include <math.h>


#include "cost.h"
#include "algorithms.h"
//#include "declares.h"

using namespace std;
cost::cost() {
}

cost::cost(const cost& orig) {
}

cost::~cost() {
    /*if (z != NULL) delete[] z;
    if (mu != NULL) delete[] mu;
    if (sigma != NULL) delete[] sigma;
    if (gamma != NULL) delete[] gamma;
    if (kappa != NULL) delete[] kappa;
    if (X != NULL){ 
        for (int i = 0; i < N; i++)
            delete[] X[i];
            delete[] X;
    }
    if (K != NULL){
        for (int i = 0; i < N; i++)
            delete[] K[i];
        delete[] K;
    }
     */
}


double cost::EmpiricalCost(){
    //according to eq 6
    //in this case num_experiments = 192
   
    double x[N];
    double kappa[M]; 
    double sum = 0;
    for (int t = 0; t < num_experiments; t++){
        for (int i = 0; i < N; i++){
            x[i] = X[i][t];
        }
        for (int s = 0; s < M; s++){
            kappa[s] = K[s][t];
        }
        
        double max = 0;
        for (int s = 0; s < M; s++){
            double term = 0;
        
            for (int i = 0; i < N; i++){
                if (z[i] == s)
                    term += x[i];
                //term += z[i][s] * x[i];
            }
            if (term > kappa[s]){
                max += term-kappa[s];
            }
        }
        sum += max;
    }
    return sum / (double)num_experiments;
   
}


void cost::setN(int n_vertex){
    if (n_vertex < 1){
        throw invalid_argument("Num of vertices must be >= 1!");
    }
    N = n_vertex;
}
void cost::setM(int n_cluster){
    if (n_cluster < 1){
        throw invalid_argument("Num of clusters must be >= 1!");
    }
    M = n_cluster;
}

void cost::setz(int* z_arr){
    
    z = z_arr;
}