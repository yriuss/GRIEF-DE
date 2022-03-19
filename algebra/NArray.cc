#include "NArray.h"
#include <iostream>
#include <string>


NArray::NArray(std::vector<int> N){
    shape = N;
    
    for(int i = 0; i < N.size(); i++)
        size *= N[i];
    data = new float[size];
}