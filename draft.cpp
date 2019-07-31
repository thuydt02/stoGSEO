/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   draft.cpp
 * Author: thuydt
 *
 * Created on July 31, 2019, 4:03 PM
 */

#include <cstdlib>
#include <iostream>
#include <vector>

#include "./lib/tFile.h"

using namespace std;

/*
 * 
 */
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
double alp_mean[] = {1, 2};
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
    
    
    
    num_edges = edgevec.size();
    num_checkins = checkinvec.size();
}
void show_dataset(){
    for (int i = 0; i < num_checkins; i++){
        std::cout <<"\n";
        for (int j = 0; j < 5 ; j++)
            std::cout <<checkinvec[i][j] << " " ;
    }   
}
int main(int argc, char** argv) {
    read_dataset();
    show_dataset();
    return 0;
}

