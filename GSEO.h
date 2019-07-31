/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   GSEO.h
 * Author: thuydt
 *
 * Created on April 12, 2019, 3:27 PM
 * https://www.cse.ust.hk/~entaflos/SIGSPATIAL17-GAME-long_version.pdf (GGSEO pp)
 * https://www.cs.ucsb.edu/~klee/papers/On_Social_Event_Organization_papers.pdf (GSEO pp)
 */
#include <queue>

#include "./lib/declare.h"
#
#ifndef GSEO_H
#define GSEO_H


using namespace std;

class GSEO {
public:
    
    GSEO();
    GSEO(const GSEO& orig);
    virtual ~GSEO();
    
    
    void setN(int num_users) {N = num_users;}
    void setM(int num_events) {M = num_events;}
    void setnum_iteration(int n) {num_iteration = n;};
    void setlambda(double l) {lambda = l;};
    
    void setc(double **dist_cost){c = dist_cost;}
    void setE(std::vector<tPair> *edges) {E = edges;}
    void setminp(int *min_cap_event) {minp = min_cap_event;}
    void setmaxp(int *max_cap_event) {maxp = max_cap_event;}

    void execute();             // = combined dynamics in GSEO pp
    void initialAssignment();
    int shifting();            // = UNI function in GSEO pp. return num of moves
    int swap();                // = BI function in GSEO pp, return num of moves
    
    void update_tcost();        // to update the tcost matrix when an assignment is made
                                // do not call this func when the z = NULL
    double getCost();
    int* getSolution(){return z;}
    
private:
    
    //given inputs area
    
    int N; //num of users
    int M; //num of events
    int num_iteration;

    double lambda; // = alpha in GSEO pp
    
    double **c ;            // size = N x M , distance from user u to event p
    std::vector<tPair> *E;   // adj list of edges with weights : friendship measurement of every 2 users w[i] 
                             //for each user, we have an adj list of edges.
                             //E = N vectors in which each element = pair(neighbor, weight)
    int *minp;              //min of capacities of events, size = M;
    int *maxp;              //max of capacities of events, size = M;
    
    //end of given inputs
    
    //need-to-be-computed area
    
    double **tCost;         //size = N x M, total cost  = assignment cost = distance cost + weight cost, need to be computed
    //MinHeap* H;             //size = N heaps, each user has its heap to store total cost of assignments => each heap has M elements
                            //the heap is to quickly find min cost, it is equivalent to tCost, but in a different order
    
    //end of need-to-be-computed
    
     //solution area
    
    
    int *z;                 //size of N, z[i] = p means user i is assigned to event p
    //end of solution area
    
    //for convenience
    std::priority_queue<tPair, std::vector<tPair>, CompareVal> getQueue(int p1, int p2);
    
    //end of for convenience
};

#endif /* GSEO_H */

