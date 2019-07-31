/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: thuydt
 *
 * Created on April 12, 2019, 2:21 PM
 */
//test github

#include <cstdlib>

#include "tFunc.h"

#include "GSEO.h"
#include "./lib/tFile.h"
#include "./lib/declare.h"

#include <cmath>
#include <iostream>
#include <math.h>
#include <cfloat>
#include <limits.h>


using namespace std;

/*
 * running GSEO on gowalla_ny dataset, including: gowalla_ny_checkins.csv, gowalla_ny_edges.csv
 * we set the weight for each edge = 1
 * 
 *
 * 
 *  
 */
const int num_iteration = 20;
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

int num_padding_users = 100;
double max_unfriend_weight = 1;

double r_user[] = {1e-5}; //Standard derivation for user locations 
int k_euclidean[] = {1, 2};

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


void read_dataset(){
    checkinvec = tFile::readDoubleFileByRow(checkins_file, ',');
    edgevec = tFile::readIntFileByRow(edges_file, ',');
    eventlocvec = tFile::readDoubleFileByRow(eventloc_file, ',');
    vertex_mappingvec = tFile::readIntFileByRow(vertex_mapping_file, ',');
    
    N = vertex_mappingvec.size(); vertex_mappingvec.clear();
    
    M = eventlocvec.size();
  
    std::ostringstream streamN; streamN << N;           sN = streamN.str();
    std::ostringstream streamM; streamM << M;           sM = streamM.str();
    num_checkins = checkinvec.size();
    num_edges = edgevec.size();
    
    /*
    int num_edges = edgevec.size();
    int num_checkins = checkinvec.size();
    
    int degree[N];
    
    for (int v = 0; v < N; v++){
        degree[v] = 0;
    }
    
    for (int i = 0; i < num_edges; i++){
        degree[edgevec[i][1]]++; degree[edgevec[i][0]]++;
    }
    
    int max_deg = 0;
    int min_deg = N;
    double avg_deg;
    
    for (int v = 0; v < N; v++){
        if (min_deg > degree[v]) min_deg = degree[v];
        if (max_deg < degree[v]) max_deg = degree[v];
        avg_deg += degree[v];
    }
    avg_deg = avg_deg / (2.0 * N);
    
    std::cout <<"\nREADING datasets";
    std::cout << ", Statistic:";
    std::cout <<"\nNum of users: N = " <<N;
    std::cout <<"\nNum_checkins = " <<num_checkins;
    std::cout <<"\nNum_edges = " <<num_edges;
    std::cout <<"\nMax_degree = " <<max_deg;
    std::cout <<"\nMin_degree = " <<min_deg;
    std::cout <<"\navg_degree = " <<avg_deg;
*/

}

void run_euclidean_gseo(){
    
    /*
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
    */
    //------------------------    
    
    //computing a[], b[] for user's locations
    //u[M], v[M] for event's locations
    
    int count[N];
    double sum_userx[N], sum_usery[N];
    double a[N], b[N], r[N];
    double u[M], v[M];
    
    for (int i = 0; i < N; i++) count [i] = 0;
    
    for (int i = 0; i < num_checkins; i++){
        int uid = int(checkinvec[i][0] + 0.1); 
        
        count[uid] += 1;
        sum_userx[uid] += checkinvec[i][2];
        sum_usery[uid] += checkinvec[i][3];
    }
    
    for (int i = 0; i < N; i++){
        a[i] = sum_userx[i]/(double)count[i];
        b[i] = sum_usery[i]/(double)count[i];
    }
    
    //computing r
    
    double sd[N];
    for (int i = 0; i < N; i++) sd[i] = 0;
    
    for (int i = 0; i < num_checkins; i++){
        int uid = int (checkinvec[i][0]);
        double dist = pow(checkinvec[i][2] - a[uid], 2) + pow (checkinvec[i][3] - b[uid], 2);
        sd[uid] += dist;
        
    }
    double min_r = r_user[0];
    double max_r = r_user[0];
    
        
    for (int i = 0; i < N; i++){
        if (count[i] > 1) {
            r[i] = sqrt(sd[i]/(count[i] - 1));
        }else
            r[i] = r_user[0];
        
        if (r[i] < 1e-5) r[i] = r_user[0];
        
        if (min_r > r[i]) min_r = r[i];
        if (max_r < r[i]) max_r = r[i];
        
    }
    
    std::cout <<"\nmin_r " <<min_r;
    std::cout <<"\nmax_r " <<max_r;
    
    for (int i = 0; i < M; i++){
                                      //eventlocvec[i].pop_back();
        v[i] = eventlocvec[i][1];//.back(); eventlocvec[i].pop_back();
        u[i] = eventlocvec[i][0];//.back(); eventlocvec[i].pop_back();
        
    }
    
    // for calling gseo
    
    //computing the friendship weight and normalize it
    std::vector<tPair> *E = new std::vector<tPair>[N];
    
    
    for (int i = 0; i < num_edges; i++){
        E[edgevec[i][0]].push_back(tPair(edgevec[i][1], 1)); //weight = 1
        E[edgevec[i][1]].push_back(tPair(edgevec[i][0], 1));
    }
    
    //scaling friendship weight such that for each i, sum(wij) in [0, 1], j = 1..N, i and j are friend
    
    double min_w = DBL_MAX;
    double max_w = DBL_MIN;
    
    for (int v = 0; v < N; v++){
        double w = 0;

        for (tPair e: E[v]){
            w += e.val;
        }
 
        if (min_w > w) min_w = w;
        if (max_w < w) max_w = w;

    }
    
    double range_w = max_w - min_w;
    
    for (int v = 0; v< N; v++)
        E[v].clear();
    
    for (int i = 0; i < num_edges; i++){
        E[edgevec[i][0]].push_back(tPair(edgevec[i][1], 1/max_w)); //weight = 1
        E[edgevec[i][1]].push_back(tPair(edgevec[i][0], 1/max_w));
    }
    //end of scaling friendship weight;
    
    
    
    //min capacity of each p and max capacity of each p
    
    int *minp = new int[M];
    int *maxp = new int[M];
    
    for (int i = 0; i < M; i++){
        minp[i] = 1;
        maxp[i] = N / M + num_padding_users;
    }
    

    GSEO gseo;
    gseo.setN(N);
    gseo.setM(M);
    gseo.setnum_iteration(num_iteration);
    gseo.setlambda(lambda);
//    gseo.setc(c);
    gseo.setE(E);
    gseo.setminp(minp);
    gseo.setmaxp(maxp);

    
    
    //computing c[N][M]
    
    double **c = new double* [N];
    double min_x = 1;
    double max_x = -0.5 * (pow(a[1] - u[1], 2) + pow(b[1] - v[1], 2)) / pow(r[1], 2);;
    
    for (int i = 0; i < N; i++){
        c[i] = new double [M];
        for (int m = 0; m < M; m++){
            double x = -0.5 * (pow(a[i] - u[m], 2) + pow(b[i] - v[m], 2)) / pow(r[i], 2);
            if (min_x > x) min_x = x;
            if (max_x < x) max_x = x;
        }
    }
    
    double range_x = abs(max_x - min_x);
    
    std::cout <<"\nmin_x, max_x, range_x " <<min_x <<", " <<max_x <<", " <<range_x;
    
    for (int k : k_euclidean){
        std::cout <<"\n------------------------------------------------------";
        std::cout <<"\neuclidean, k = " <<k; 
        
        std::ostringstream streamk; streamk << k;           
        std::string sk = streamk.str();
    
        if (k == 1){
            for (int i = 0; i < N; i++){
                for (int m = 0; m < M; m++){
                    double x = -2 * 0.5 * ((pow(a[i] - u[m], 2) + pow(b[i] - v[m], 2)) / pow(r[i], 2)) / range_x;
                    c[i][m] = r[i] * sqrt(0.5 * M_PI) * tFunc::fL(x);
                }
            }
        }else if (k == 2){
            for (int i = 0; i < N; i++){
                c[i] = new double [M];
                for (int m = 0; m < M; m++){
                    c[i][m] = pow(a[i] - u[m], 2) + pow(b[i] - v[m], 2) + 2.0 * pow(r[i], 2);
                }
            }
        }else{

        }

        //scaling c[i][m] in [0, 1]
        double max_c = c[0][0];
        double min_c = c[0][0];

        for (int i = 0; i < N; i++)
            for (int m = 0; m < M; m++){
                if (max_c < c[i][m]) max_c = c[i][m];
                if (min_c > c[i][m]) min_c = c[i][m]; 
            }
        double nmin_c = DBL_MAX;
        double nmax_c = DBL_MIN;

        for (int i = 0; i < N; i++)
            for (int m = 0; m < M; m++){
                c[i][m] = (c[i][m] - min_c)/max_c;
                if (nmax_c < c[i][m]) nmax_c = c[i][m];
                if (nmin_c > c[i][m]) nmin_c = c[i][m]; 

            }

        std::cout <<"\nMin_c = " <<nmin_c;
        std::cout <<"\nMax_c = " <<nmax_c;
        gseo.setc(c);
        gseo.execute();
        
        std:string fout = out_folder + "z_BR_k" + sk + "_N" + sN + "_M" + sM + ".csv";
        std::cout <<"\nsaving to " <<fout;
        tFile::saveIntFile(gseo.getSolution(), N, fout, '\n');
    }
}

void run_RndPartition(){

    std::cout <<"\nrunning RandomPartition of " << N << " users in ";
    std::cout <<M <<" events";
    
    int* z = new int [N];
    int inCluster[N];
    for (int i = 0; i < N; i++)
        inCluster[i] = 0;
    
    srand(time(NULL));
    int count = 0;
    
    while (true){
        
        int s = rand() % M;
        int i = rand() % N;
        
        if (inCluster[i] == 1)
            continue;
        
        inCluster[i] = 1;
        z[i] = s;
        
        count++;
        if (count == N) break;
    }
    
    
    
    std::string z_file = out_folder + "z_RndPartition_N" + sN + "_M" + sM + ".csv";
    std::cout <<"\nSaving z into " << z_file << "\n";
    tFile::saveIntFile(z, N, z_file, '\n');
    if (z != NULL) delete[] z;
}

int main()
{
    
    read_dataset();
    //std::cout <<"\nrunning euclidean_gseo";
    //N = 5620; sN = "5620";
    //M = 20; sM = "20";
    //run_RndPartition();
    run_euclidean_gseo();
    //std::cout <<"\ndone";
   
    
    
     //clearance
    
   
           
}



