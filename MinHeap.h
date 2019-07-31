/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   MinHeap.h
 * Author: thuydt
 *
 * Created on April 12, 2019, 2:25 PM
 */

#include <sstream>

#ifndef MINHEAP_H
#define MINHEAP_H

#include "./lib/declare.h"

class MinHeap {

public:
    MinHeap(){
        heap_size = 0;
        capacity = 0;
        harr = NULL;
    };
    
    void setcapacity(int cap){
        heap_size = 0; 
        capacity = cap; 
        harr = new tPair[cap];
    }
    
    MinHeap(int capacity); 

    MinHeap(const MinHeap& orig);
    
    MinHeap& operator = (const MinHeap &orig){ 
        
        //copy assignment
        
        heap_size = orig.heap_size; 
        capacity = orig.capacity;
        
        if (this->harr != NULL){
            delete[] this->harr;    
        }
        
        harr = new tPair[capacity];
        for (int i = 0; i < heap_size; i++){
            harr[i] = orig.harr[i];
        }
        return *this;
    }
  
    
    
    virtual ~MinHeap();
      
    // to heapify a subtree with the root at given index 
    void MinHeapify(int ); 
  
    int parent(int i) { return (i-1)/2; } 
  
    // to get index of left child of node at index i 
    int left(int i) { return (2*i + 1); } 
  
    // to get index of right child of node at index i 
    int right(int i) { return (2*i + 2); } 
  
    // to extract the root which is the minimum element 
    tPair extractMin(); 
  
    // Decreases key value of key at index i to new_val 
    void decreaseKey(int i, tPair new_val); 
  
    // increases key value of key at index i to new_val 
    void increaseKey(int i, tPair new_val); 
  
    // Returns the minimum key (key at root) from min heap 
    tPair getMin() { return harr[0]; } 
  
    // Deletes a key stored at index i 
    void deleteKey(int i); 
  
    // Inserts a new key 'k' 
    void insertKey(tPair k); 

    void swap(tPair *x, tPair *y);
    
    //----------------------------------------------
    void getIndex(tPair p, int from, bool* found, int* result); 
    //return index of element which has key = p.key, searching starts at from
    void getIndexLinear(tPair p, int from, bool* found, int* result);
    int getheap_size(){return heap_size;};
    
    void clearAll(){
        delete[] harr;
        harr = NULL;
        heap_size = 0;
        capacity = 0;
    };
    
    std::string toString() const
    {
        std::ostringstream oss;
        oss << capacity << "," << heap_size << "," << harr << "\n";
        return oss.str();
    }
//private:
 
    tPair *harr;  // pointer to array of elements in heap 
    int capacity; // maximum possible size of min heap 
    int heap_size; // Current number of elements in min heap 

    
};

#endif /* MINHEAP_H */

