#include "NArray.h"
#include <iostream>

NArray::NArray(std::vector<int> N){
    shape = N;
}


int main(){
    std::vector<int> dimensions;

    dimensions.push_back(512);
    dimensions.push_back(4);
    
    NArray array(dimensions);
    std::cout << array.shape;
    array.shape = dimensions;
    std::cout << array.shape;
    //std::cout << array.shape[0] << std::endl;


    return 0;
}