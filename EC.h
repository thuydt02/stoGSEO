/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   EC.h
 * Author: thuydt
 *
 * Created on July 30, 2019, 4:36 PM
 */

#ifndef EC_H
#define EC_H

#include "./lib/declare.h"
class EC {
public:
    EC();
    EC(int n, int m, std::vector<tPair>* edge,  double al, double ld, std::vector<double>* xp, 
        std::vector<double>* yp, double* uu, double *vv, int* zz);
    EC(const EC& orig);
    virtual ~EC();
    
    int N;
    int M;
    double alp;
    double lambda; // = alpha in GSEO pp
    std::vector<tPair> *E;   // adj list of edges with weights : friendship measurement of every 2 users w[i] 
                             //for each user, we have an adj list of edges.
                             //E = N vectors in which each element = pair(neighbor, weight)
    
    double* u;
    double* v;
    std::vector<double>* X;
    std::vector<double>* Y;
    int* z;
    //void setup();
    double getEC();
    
private:

};

#endif /* EC_H */

