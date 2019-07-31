/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   MinHeap.cpp
 * Author: thuydt
 * 
 * Created on April 12, 2019, 2:25 PM
 */
// A C++ program to demonstrate common Binary Heap Operations 

#include<iostream> 
#include<climits> 
#include <cfloat>

#include "MinHeap.h"
#include "./lib/declare.h"

using namespace std; 

MinHeap::MinHeap(const MinHeap& orig) {
    
    //constructor assignment
    
    heap_size = orig.heap_size; 
    capacity = orig.capacity; 
    harr = new tPair[capacity];
    for (int i = 0; i < heap_size; i++){
        harr[i] = orig.harr[i];
    }
}


MinHeap::~MinHeap() {
    if (harr != NULL){
        delete[] harr;
        harr = NULL;
    }
}
  
// Constructor: Builds a heap from a given array a[] of given size 
MinHeap::MinHeap(int cap) 
{ 
    heap_size = 0; 
    capacity = cap; 
    harr = new tPair[cap]; 
} 
  
// Inserts a new key 'k' 
void MinHeap::insertKey(tPair k) 
{ 
    if (heap_size == capacity) 
    { 
        cout << "\nOverflow: Could not insertKey\n"; 
        return; 
    } 
  
    // First insert the new key at the end 
    heap_size++; 
    int i = heap_size - 1; 
    harr[i] = k; 
  
    // Fix the min heap property if it is violated 
    while (i != 0 && harr[parent(i)].val > harr[i].val) 
    { 
       swap(&harr[i], &harr[parent(i)]); 
       i = parent(i); 
    } 
} 
  
// Decreases value of key at index 'i' to new_val.  It is assumed that 
// new_val is smaller than harr[i]. 
void MinHeap::decreaseKey(int i, tPair new_val) 
{ 
    harr[i] = new_val; 
    while (i != 0 && harr[parent(i)].val > harr[i].val) 
    { 
       swap(&harr[i], &harr[parent(i)]); 
       i = parent(i); 
    } 
} 
  
// increases value of key at index 'i' to new_val.  It is assumed that 
// new_val is bigger than harr[i]. 
void MinHeap::increaseKey(int i, tPair new_val) 
{ 
    harr[i] = new_val;
    int k = i;
    while (k <= heap_size)
    { 
        if ((harr[k].val >= harr[left(k)].val) && (harr[left(k)].val >= harr[right(k)].val)){
            swap(&harr[k], &harr[right(k)]);
            k = right(k);
        }else
        if ((harr[k].val >= harr[right(k)].val) && (harr[right(k)].val >= harr[left(k)].val)){
            swap(&harr[k], &harr[left(k)]);
            k = left(k);
        }else
        if ((harr[left(k)].val >= harr[k].val) && (harr[k].val >= harr[right(k)].val)){
            swap(&harr[k], &harr[right(k)]);
            k = right(k);
        }else
        if ((harr[right(i)].val >= harr[i].val) && (harr[i].val >= harr[left(i)].val)){
            swap(&harr[k], &harr[left(k)]);
            k = left(k);
        }else
            break;
    } 
} 
  

// Method to remove minimum element (or root) from min heap 
tPair MinHeap::extractMin() 
{ 
    tPair ret;
    ret.key = -1;
    ret.val = DBL_MAX;
    if (heap_size <= 0) 
        return ret; 
    if (heap_size == 1) 
    { 
        heap_size--; 
        return harr[0]; 
    } 
  
    // Store the minimum value, and remove it from heap 
    tPair root = harr[0]; 
    harr[0] = harr[heap_size-1]; 
    heap_size--; 
    MinHeapify(0); 
  
    return root; 
} 
  
  
// This function deletes key at index i. It first reduced value to minus 
// infinite, then calls extractMin() 
void MinHeap::deleteKey(int i) 
{ 
    tPair t(-1, DBL_MIN );
    decreaseKey(i, t); 
    extractMin(); 
} 
  
// A recursive method to heapify a subtree with the root at given index 
// This method assumes that the subtrees are already heapified 
void MinHeap::MinHeapify(int i) 
{ 
    int l = left(i); 
    int r = right(i); 
    int smallest = i; 
    if (l < heap_size && harr[l].val < harr[i].val) 
        smallest = l; 
    if (r < heap_size && harr[r].val < harr[smallest].val) 
        smallest = r; 
    if (smallest != i) 
    { 
        swap(&harr[i], &harr[smallest]); 
        MinHeapify(smallest); 
    } 
} 
  
// A utility function to swap two elements 
void MinHeap::swap(tPair *x, tPair *y) 
{ 
    tPair temp;
    
    temp.key = x->key;
    temp.val = x->val;
    
    x->key = y->key;
    x->val = y->val;
    
    y->key = temp.key;
    y->val = temp.val;
} 


void MinHeap::getIndex(tPair p, int from, bool *found, int *result){
    
    if (heap_size == 0) return;
    
    if (from > heap_size) return;
    
    if (harr[from].val > p.val) 
    { 
        /* Skip this node and its descendants, 
          as they are all >= p.val . */
        return; 
    } 
    
    if (*found == true) return;
    //if (harr[from].key == p.key) {*found = true;  *result = from; return;};
    
    if (harr[from].val == p.val) {*found = true;  *result = from; return;};
    
    getIndex(p, left(from), found, result); 
    
    getIndex(p, right(from), found, result); 
}

void MinHeap::getIndexLinear(tPair p, int from, bool *found, int *result){
    
    if (heap_size == 0) return;
    
    if (from > heap_size) return;
    
    for (int i = 0; i < heap_size; i++){
        if (harr[i].key == p.key){
            *found = true;
            *result = i;
            return;
        }
    }
}