/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   computeEC.cpp
 * Author: thuydt
 *
 * Created on July 30, 2019, 4:43 PM
 */

#include <cstdlib>
#include <vector>
#include "cost.h"
#include "./lib/tFile.h"
#include "EC.h"

using namespace std;

/*
 * 
 */
double lambda = 0.5;
std::string base_folder = "./data_run/configuration1/";
std::string input_folder = base_folder + "input/";
//std::string checkins_file = input_folder + "gowalla_dl_at_normxy_gt9_checkins.csv";//gowalla_ny_normxy_checkins.csv";
//std::string edges_file = input_folder + "gowalla_dl_at_gt9_edges.csv";//gowalla_ny_edges.csv";
//std::string eventloc_file = input_folder + "gowalla_dl_at_normxy_eventloc_M200.csv";//gowalla_ny_normxy_eventloc_M20.1";
//std::string vertex_mapping_file = input_folder + "gowalla_dl_at_gt9_uid_mapping.csv";//gowalla_ny_vertex_mapping.csv";

std::string checkins_file = input_folder + "gowalla_ny_normxy_checkins.csv";
std::string edges_file = input_folder + "gowalla_ny_edges.csv";
std::string eventloc_file = input_folder + "gowalla_ny_normxy_eventloc_M20.1";
std::string vertex_mapping_file = input_folder + "gowalla_ny_vertex_mapping.csv";

std::string out_folder = base_folder + "output/";

std::string eval_folder = base_folder + "evaluation/";

double alp[] = {1, 2};
//for being ease to use

std::vector<std::vector<double> > eventlocvec;
std::vector<std::vector<double> > checkinvec;
std::vector<std::vector<int> > edgevec;
std::vector<std::vector<int> > vertex_mappingvec;

int N; //num of users
int M; //num of events
int num_checkins;
int num_edges;
std::string sN;
std::string sM;



int main(int argc, char** argv) {

    checkinvec = tFile::readDoubleFileByRow(checkins_file, ',');
    edgevec = tFile::readIntFileByRow(edges_file, ',');
    eventlocvec = tFile::readDoubleFileByRow(eventloc_file, ',');
    vertex_mappingvec = tFile::readIntFileByRow(vertex_mapping_file, ',');
    
    N = vertex_mappingvec.size(); vertex_mappingvec.clear();
    M = eventlocvec.size();
    
    std::ostringstream streamN; streamN << N;           sN = streamN.str();
    std::ostringstream streamM; streamM << M;           sM = streamM.str();
    
    int num_checkins = checkinvec.size();
    int num_edges = edgevec.size();
    
    //------------------------    
    
    //computing X[N] of vectors and Y[N] of vectors
    //u[M], v[M] for event's locations
    
    
    std::vector<double>* X = new std::vector<double>[N];
    std::vector<double>* Y = new std::vector<double>[N];
    
    double* u = new double [M]; double* v = new double[M];
    
    for (int i = 0; i < num_checkins; i++){
        int uid = int(checkinvec[i][0] + 0.1); 
        
        X[uid].push_back(checkinvec[i][2]);
        Y[uid].push_back(checkinvec[i][3]);
    }
    
    
    for (int i = 0; i < M; i++){
                                      //eventlocvec[i].pop_back();
        v[i] = eventlocvec[i][1];//.back(); eventlocvec[i].pop_back();
        u[i] = eventlocvec[i][0];//.back(); eventlocvec[i].pop_back();
        
    }
    
    // for calling gseo
    
    std::vector<tPair> *E = new std::vector<tPair>[N];
    
    
    for (int i = 0; i < num_edges; i++){
        E[edgevec[i][0]].push_back(tPair(edgevec[i][1], 1)); //weight = 1
        E[edgevec[i][1]].push_back(tPair(edgevec[i][0], 1));
    }
    
    //void EC::EC(int n, int m, std::vector<tPair>* edge,  double al, double ld, std::vector<double>* xp, 
    //    std::vector<double>* yp, double* uu, double *vv, int* zz){
    
    
    std::string alg_name[] = {"z_RndPartition", "z_BRMean", "z_BR"};
    int num_alg = 3;
    std::vector<int> z_rnd_vec = tFile::readIntFile1Column(out_folder + alg_name[0] + "_N" + sN + "_M" + sM + ".csv");
    
    std::vector<int> zvec;
    
    EC eCost(N, M, E, 0, lambda, X, Y, u, v, NULL);
    
    std::vector<std::vector<double> > costsave;
    std::cout <<"\n-----------------------------------------";
    std::cout <<"Computing Empirical Cost";
    std::cout <<"\nz_rnd_vec.size() " << z_rnd_vec.size();
    for (double al: alp){
        eCost.alp = al;
        std::ostringstream streamk; streamk << al;           
        string sk = streamk.str();
    
        std::vector<double> costvec;
        
        costvec.push_back(al);
        
        for (int i = 0; i < num_alg; i++ ){
            std::string z_file = out_folder + alg_name[i] + "_k" + sk + "_N" + sN + "_M" + sM + ".csv";
            if (i > 0){
                zvec = tFile::readIntFile1Column(z_file);
                std::cout <<"\nReading " <<z_file;
            }else{
                std::cout <<"\nReading " <<out_folder + alg_name[0] + "_N" + sN + "_M" + sM + ".csv";
                zvec = z_rnd_vec;
            }
            int* z = new int[N];
            double ec = -1;
            int u = 0;
            for (int s: zvec){
                if (s > -1){
                    z[u] = s; 
                    u++;     
                }else{
                    break;
                }
            }
            if (u == N){
                eCost.z = z;
                ec = eCost.getEC();
                eCost.z = NULL;
            }
            if (z != NULL){ 
               delete[] z;
               z = NULL;
            }

            costvec.push_back(ec);
        } //algorithm
        costsave.push_back(costvec);
    }//alpha
    std::string fout = eval_folder + "EC.csv";
    std::cout <<"\nsaving to " <<fout;
    tFile::saveDoubleFile(costsave, fout, ',');
    
}

