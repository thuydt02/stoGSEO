/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   cost.h
 * Author: thuydt
 *
 * Created on February 19, 2019, 4:17 PM
 */
//to compute the cost of StoMEC
#ifndef COST_H
#define COST_H
#include <iostream>


class cost {
public:
    cost();
    cost(const cost& orig);
    virtual ~cost();
    double EmpiricalCost();
    void setN(int n_vertex);
    void setM(int n_cluster);
    void setz(int* z_arr);
    
private:
    int N;
    int M;
    int* z; //
};

#endif /* COST_H */

