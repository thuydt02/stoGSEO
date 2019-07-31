/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   GSEO.cpp
 * Author: thuydt
 * 
 * Created on April 12, 2019, 3:27 PM
 */
#include <iostream>
#include <cfloat>
#include <algorithm>
#include <functional>
#include <queue>
#include <unordered_set>

//#include <vector>

#include "GSEO.h"

const int num_S = 30; // randomly select distinct num_S users from Vun to assign 
using namespace std;

GSEO::GSEO() {
    tCost = NULL;
    z = NULL;
    //H = NULL;
}

GSEO::GSEO(const GSEO& orig) {
}

GSEO::~GSEO() {
    
    if (c != NULL){
        for (int i = 0; i < N; i++){
            delete[] c[i];
        }
        delete[] c;
    }
    
    if (tCost != NULL) {
        for (int i = 0; i < N; i++){
            delete[] tCost[i];
        }
        delete[] tCost;
    
    }

    if (minp != NULL)   delete[] minp;
    if (maxp != NULL)   delete[] maxp;
    if (E != NULL)      delete[] E; 
    //if (H != NULL)      delete[] H;
    if (z != NULL)      delete[] z;
}

void GSEO::execute(){             // = combined dynamics in GSEO pp
    
    initialAssignment();        //resulted in z
    
    int num_shifting = 0;   int num_swap = 0;
    int num_iter = 0;
    
    std::cout <<"\nin GSEO::execute()";
    
    while (true){
        
        std::cout <<"\niteration: " <<num_iter;
        
        num_shifting = shifting();                 //start with z, resulted in z
    
        std::cout <<"\nnum_shifting, the cost: " <<num_shifting <<", " <<getCost();
    
        num_swap = swap();
        
        std::cout <<"\nnum_swap, the cost: " <<num_swap <<", " <<getCost();
        
        if ((num_shifting == 0) && (num_swap == 0))
            break;
        
        num_iter++;
        
    } 
                    
}
void GSEO::initialAssignment(){
    
    ivector V_un;
    ivector P_op;
    ivector S;
    
    int count [M];      //num of users in each event
    
    //initialize
    
    if (z == NULL) z = new int[N];
    
    if (tCost == NULL) tCost = new double* [N];
    
    for (int v = 0; v < N; v++){
        V_un.push_back(v);                              //all users is not assigned
        tCost[v] = new double[M];
    }
    
    for (int p = 0; p < M; p++){     
        P_op.push_back(p);                              //all events is opened
        count[p] = 0;
    }
    
    double min_tcost = DBL_MAX;
    double max_tcost = DBL_MIN;
    double min_w = DBL_MAX;
    double max_w = DBL_MIN;
    
    std::priority_queue<tPair, std::vector<tPair>, CompareVal> Q[N];
    
    for (int v = 0; v < N; v++){
        for (int p = 0; p < M; p++){
            
            tCost[v][p] = lambda * c[v][p];
            double w = 0;
            
            for (tPair e: E[v]){
                w += e.val;
            }
            tCost[v][p] += (1 - lambda) * w;
            
            Q[v].push(tPair(p, tCost[v][p]));
            
            
            if (min_tcost > tCost[v][p]) {min_tcost = tCost[v][p];};
            
            if (max_tcost < tCost[v][p]) max_tcost = tCost[v][p];
            
            if (min_w > w) min_w = w;
            if (max_w < w) max_w = w;
        }
    }
    
    //phase 1
    std::cout <<"\nmin_tcost " <<min_tcost;
    std::cout <<"\nmax_tcost " <<max_tcost;
    std::cout <<"\nmin_w " <<min_w;
    std::cout <<"\nmax_w " <<max_w;
    
    //exit(1);
    
    srand (time(NULL));
    
    bool *in_S;
        
    while (P_op.size() != 0){
        
        int V_un_size = V_un.size();
        if (V_un_size == 0) break;
        
        S.clear();
            
        if (V_un_size > num_S){
            
            in_S = new bool[V_un_size];
            
            for (int t = 0; t < V_un_size; t++) in_S[t] = false;
            int i = 0;
            
            while (i < num_S) {
                
                int pos = rand() % V_un_size;
                
                if (in_S[pos] == false){
                    S.push_back(V_un[pos]);
                    in_S[pos] = true;
                    i++;
                }
            }
            delete[] in_S;
        }else{
            S = V_un;
        }
        
        double cm = DBL_MAX;        
        int vprime = -1;
        int pprime = -1;
    //    std::cout <<"\nS.size() " <<S.size();
        for (int v: S){
            //find out the min-cost user to assign
            
            
            tPair min_vp = Q[v].top(); 
            
            int p_star = min_vp.key;        //p_star is an event with min cost for v in P_op
            
            if (cm > tCost[v][p_star]){
                cm = tCost[v][p_star];
                vprime = v;
                pprime = p_star;
            }    
        }
        
    //    std::cout <<"\nvprime, pprime " <<vprime <<", " <<pprime;
        z[vprime] = pprime; //assign vprime to pprime
        
        count[pprime]++;
        
        std::vector<int>::iterator pos = std::find(V_un.begin(), V_un.end(), vprime);
        //if (pos != V_un.end())
            V_un.erase(pos);
        
        
        for (tPair e: E[vprime]){
            
            int f = e.key;
            
            pos = std::find(V_un.begin(), V_un.end(), f);
            if (pos == V_un.end())  continue;   //f was assigned
        
            tCost[f][pprime] -= (1 - lambda)*e.val;
            
            Q[f] = std::priority_queue<tPair, std::vector<tPair>, CompareVal>();
            
            for (int ev = 0; ev < M; ev++){
                Q[f].push(tPair(ev, tCost[f][ev]));
            }
        }
        
        if (count[pprime] == minp[pprime]){ //pprime needs to be closed
            pos = std::find(P_op.begin(), P_op.end(), pprime);
            
            if (pos != P_op.end()){
               P_op.erase(pos);
            }
            
            for (int v: V_un){ //remove pprime from Q[v]

                Q[v] = std::priority_queue<tPair, std::vector<tPair>, CompareVal>();
            
                for (int ev = 0; ev < M; ev++){
                    if (ev == pprime) continue;
                    Q[v].push(tPair(ev, tCost[v][ev]));
                }
            }
        }
    }
    
    //end of phase 1
    
    //phase 2
    
    //rebuild the queue
    
    for (int v: V_un){
        Q[v] = std::priority_queue<tPair, std::vector<tPair>, CompareVal>();        
        
        for (int p = 0; p < M; p++){
            
            Q[v].push(tPair(p, tCost[v][p]));
        }
    }
    
    std::cout <<"\nINIT end of phase1";
    
    //reopen all events
    
    P_op.clear();
    
    for (int p = 0; p < M; p++)
        P_op.push_back(p);
    
    //loop until V_un = empty
    //assume that sum of capacities of all event >= num of users. So we can assign all users
    //exit(1);
    
    while (true){
        
        int V_un_size = V_un.size();
        //std::cout <<"\nV_un.size() = " <<V_un_size;
    
        if (V_un_size == 0) break;
        
        S.clear();
            
        if (V_un_size > num_S){
            
            in_S = new bool[V_un_size];
            
            for (int t = 0; t < V_un_size; t++) in_S[t] = false;
            int i = 0;
            
            while (i < num_S) {
                
                int pos = rand() % V_un_size;
                
                if (in_S[pos] == false){
                    S.push_back(V_un[pos]);
                    in_S[pos] = true;
                    i++;
                }
            }
            delete[] in_S;
        }else{
            S = V_un;
        }
        
        double cm = DBL_MAX;        
        int vprime = -1;
        int pprime = -1;
        
        for (int v: S){
            
            //if (Q[v].empty()) {std::cout << "v = " <<v <<", Q[v] is empty"; exit(1);};
            
            tPair min_vp = Q[v].top(); //Q[v].pop(); //H[v].getMin();
          
            int p_star = min_vp.key;               //p_star is an event with min cost for v in P_op
            
            if (cm > tCost[v][p_star]){
                cm = tCost[v][p_star];
                vprime = v;
                pprime = p_star;
            }    
        }
        //std::cout <<"\n vprime, pprime : " <<vprime << ", " <<pprime;
        z[vprime] = pprime;
        
        count[pprime]++;
        
        std::vector<int>::iterator pos = std::find(V_un.begin(), V_un.end(), vprime);
        //if (pos != V_un.end())
            V_un.erase(pos);
        
        
        for (tPair e: E[vprime]){
            int f = e.key;
            
            pos = std::find(V_un.begin(), V_un.end(), f);
            if (pos == V_un.end())  continue; //f was assigned
            
            tCost[f][pprime] -=  (1 - lambda) * e.val;
            
            
            Q[f] = std::priority_queue<tPair, std::vector<tPair>, CompareVal>();
            
            for (int ev = 0; ev < M; ev++){
                Q[f].push(tPair(ev, tCost[f][ev]));
            }
        
        }
        
        if (count[pprime] == maxp[pprime]){ //pprime needs to be closed
            
            pos = std::find(P_op.begin(), P_op.end(), pprime);
     
            //if (pos != P_op.end()){
               P_op.erase(pos);
            //}
            
            for (int v: V_un){ //remove pprime from Q[v]

                Q[v] = std::priority_queue<tPair, std::vector<tPair>, CompareVal>();

                for (int ev = 0; ev < M; ev++){
                    if (ev == pprime) continue;
                    Q[v].push(tPair(ev, tCost[v][ev]));
                }
            }
       
        }
    }//of V_un == empty
    
    //end of phase 2
    //clear all
    
    //update_tcost();
    std::cout <<"\nEnd of initialAssignment!";
    std::cout <<"\nAssignment statistic:";
    
    std::cout <<"\nThe cost = " <<getCost();
    
    int s = 0;
    for (int p = 0; p < M; p++){
        std::cout <<"\ncount[" <<p <<"]= " <<count[p];
        s += count[p];
    }
    std::cout<<"\nN, sum_of assigned users = " <<N <<", " <<s;
}


int GSEO::shifting(){            // = UNI function in GSEO pp
    
    
    
    int count[M];
    
    for (int p = 0; p < M; p++)     count[p] = 0;
    
    for (int v = 0; v < N; v++){
        count[z[v]] ++;
    }
      
    std::priority_queue<tPair, std::vector<tPair>, CompareVal > Q[N];
  
    for (int v = 0; v < N; v++){
        for (int p = 0; p < M; p++){
            Q[v].push(tPair(p, tCost[v][p]));
        }
    }
    
    int num_move = 0;
    int sub_num_move;
    
    std::priority_queue<tPair, std::vector<tPair>, CompareVal> Qv;
    
    for (int t = 0 ; t < num_iteration; t++){
    
        sub_num_move = 0;
        
        for (int v = 0; v < N; v++){
            
            int p = z[v];
            
            if (count[p] > minp[p]){ // satisfy the constraint to move v from p to pprime 
                
                Qv = Q[v];
                
                while (true){
                
                    if (Qv.empty()) break;
                    
                    tPair min = Qv.top(); Qv.pop();
                    
                    int pprime = min.key;
                    
                    if (p == pprime) break;
                    
                    if (count[pprime] < maxp[pprime]){
                        
                        //std::cout <<"\nv, p, pprime: " <<v <<", " <<p << ", " <<pprime;
                        //std::cout <<"\nQv.size, top " <<Qv.size() <<", "<<Qv.top().key <<", " <<Qv.top().val;
                        
                        //std::cout <<"\nmin = " <<min.key <<", " <<min.val;
                        double c1 = getCost();
                        //std::cout <<"\nBEFORE shifting: The cost = "<< c1;
                        z[v] = pprime;
                        double c2 = getCost();
                        //std::cout <<"\nAFTER shifting: The cost = "<< c2;
                        
                        if ((c1 - c2) < 0){
                            z[v] = p;
                            //std::cout <<"=> no shifting";
                            break;
                        }
                        
                        count[p]--;
                        count[pprime]++;

                        for (tPair e: E[v]){
                        
                            int f = e.key;
                            
                            tCost[f][p] +=  (1 - lambda) * e.val;
                            tCost[f][pprime] -=  (1 - lambda) * e.val;
                            
                            Q[f] = std::priority_queue<tPair, std::vector<tPair>, CompareVal>();
                            
                            for (int ev = 0; ev < M; ev++)
                                Q[f].push(tPair(ev, tCost[f][ev]));
                        }
                        sub_num_move++;
                        break;
                    } //while true
                } //if count[p] > minp
            } //for v
        } //for t
        num_move += sub_num_move;
        if (sub_num_move == 0) break;
    }
    return num_move;

}
                                 
int GSEO::swap(){                // = BI function in GSEO pp
    
    int count[M];
    
    for (int i = 0; i < M; i++)     count[i] = 0;
    for (int i = 0; i < N; i++){
        count[z[i]]++;
    }
    
    //initialize
    
    //M priority queues. Each is equivalent to an event
    //
    
    std::priority_queue<tPair, std::vector<tPair>, CompareVal> Q[M][M];
    
    for (int i = 0; i < M; i++)
        for (int j = 0; j < M; j++){
            if (i == j) continue;
            for (int v = 0; v < N; v++)
                if (z[v] == i)
                    Q[i][j].push(tPair(v, tCost[v][j] - tCost[v][i]));
        }
    
    
    //main part
    
    int sub_num_move;
    int num_move = 0;
    
    std::cout <<"\nin swap()";
    int t = 0;
    while (true){
        
        sub_num_move = 0;
        
        for (int i = 0; i < M; i++){
            for (int j = 0; j < M; j++){
                
                if (i == j)     continue;
                
                tPair minpi = Q[i][j].top();    int v = minpi.key;
                tPair minpj = Q[j][i].top();    int u = minpj.key;
                double deltav = minpi.val;      double deltau = minpj.val;
                
                bool is_friend = false;
                double weightvu;
                
                if (E[v].size() < E[u].size()){
                    for (tPair e: E[v]){
                        if (e.key == u){
                            is_friend = true;
                            weightvu = e.val;
                            break;
                        }
                    }
                }else{
                    for (tPair e: E[u]){
                        if (e.key == v){
                            is_friend = true;
                            weightvu = e.val;
                            break;
                        }
                    }
                }
                
                if (is_friend){
                    deltav +=  (1-lambda) * weightvu;
                    deltau +=  (1-lambda) * weightvu;
                }
                
                //if (deltav + deltau < 0){
                if ((deltav < 0) && (deltau < 0)){ //strong stability
                    
                    //std::cout <<"\ni,j v,u: " <<i <<"," <<j <<", "<<v << ", " <<u;
                    double c1 = getCost();
                    //std::cout <<"\nBF swap(), the cost = " <<c1;
                    if ((i != z[v]) || (j != z[u])){
                        std::cout <<"\n*****************v, u, i, j, z[v], z[u] " <<v <<", " <<u <<", "<<i <<", " <<j;
                        std::cout <<", " <<z[v] <<", " <<z[u];
                        exit(1);
                    }
                    z[v] = j;   z[u] = i;
                    double c2 = getCost();
                    //std::cout <<"\nAT swap(), the cost = " <<c2;
                    
                    if (c1 < c2){
                        //std::cout <<" => no swap";
                        z[v] = i; z[u] = j;
                        continue;
                    }
                  //  std::cout <<"\n-------------------t, i,j v,u: " <<t <<", " <<i <<"," <<j <<", "<<v << ", " <<u;
                    
                    sub_num_move++;
    
                    std::unordered_set<int> friend_of_vu;
    
                    for (tPair e: E[v]){
                        
                        int f = e.key;
                        
                        friend_of_vu.insert(f);
                        
                        tCost[f][j] -=  (1-lambda)*e.val;
                        
                        tCost[f][i] +=  (1-lambda)*e.val;
                        
                    }
                    
                    for (tPair e: E[u]){
                        
                        int f = e.key;
                        
                        friend_of_vu.insert(f);
                        
                        tCost[f][j] +=  (1-lambda)*e.val;
                        
                        tCost[f][i] -=  (1-lambda)*e.val;
                        
                    }
                    
                    if (is_friend){
                        
                        tCost[v][i] -=  (1 - lambda) * weightvu; // because u ->i
                        tCost[u][i] +=  (1 - lambda) * weightvu; // because v -> j, earlier v in i
                        tCost[v][j] +=  (1 - lambda) * weightvu; //because u -> i, earlier u in j
                        tCost[u][j] -=  (1 - lambda) * weightvu; //because v ->j
                    }
                    
                    for (std::unordered_set<int>::iterator it = friend_of_vu.begin(); it != friend_of_vu.end(); it++){
                        Q[z[*it]][i] = getQueue(z[*it], i);
                        Q[z[*it]][j] = getQueue(z[*it], j);
                    }
                    
                    for (int p1 = 0; p1 < M; p1++)
                        Q[i][p1] = getQueue(i,p1);
                    
                    for (int p1 = 0; p1 < M; p1++)
                        Q[j][p1] = getQueue(j,p1);
                    
                    }
            }
        }
        t++;
        if (sub_num_move == 0) break;
        num_move += sub_num_move;
    }
    
    //clearance
    
    return num_move;
}

double GSEO::getCost(){
    
    double t, tmp;
    int f;
    double ret = 0;
    
    for (int v = 0; v < N; v++){
        
        tmp = lambda * c[v][z[v]];
        t = 0;
        for (tPair e: E[v]){
            f = e.key;
            if (z[v] != z[f])
                t += e.val;
        }
        tmp += (1-lambda) *  t;  
        ret += tmp;
    }
    return ret;
    
}

void GSEO::update_tcost(){
    
    if (z == NULL){
        std::cout <<"\nThe assignment variable z == NULL, cannot update tcost!";
        return;
    }
    
    //update individual tcost according to the assignment z, see eq (4) in the Game theory SEO pp
    //here we assume that user v is assigned to event p
    
    double t;
    int f;
    
    for (int v = 0; v < N; v++){
        for (int p = 0; p < M; p++){
            tCost[v][p] = lambda * c[v][p];
            t = 0;
            for (tPair e: E[v]){
                f = e.key;
                if (p != z[f])
                    t += e.val;
            }
            tCost[v][p] += (1-lambda) *  t;   
            //std::cout <<tCost[v][p] << "  ";
        }     
    }
}

std::priority_queue<tPair, std::vector<tPair>, CompareVal> GSEO::getQueue(int p1, int p2){
    
    std::priority_queue<tPair, std::vector<tPair>, CompareVal> Q;
    
    if (p1 == p2) return Q;
    for (int v = 0; v < N; v++)
        if (z[v] == p1)
            Q.push(tPair(v, tCost[v][p2] - tCost[v][p1]));
    return Q;
}
