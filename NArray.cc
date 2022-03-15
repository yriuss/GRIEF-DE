#include "NArray.h"
#include <iostream>
#include <string>
NArray::NArray(std::vector<int> N){
    shape = N;
    
    for(int i = 0; i < N.size(); i++)
        size *= N[i];
    std::cout << "size is " << size << std::endl;
    std::cout << "dimension is " << shape.length() << std::endl;
    data = new int[size];
}


int main(){
    std::vector<int> dimensions;

    dimensions.push_back(4);
    dimensions.push_back(4);
    

    NArray array(dimensions);
    //std::cout << array.shape;
    
    
    array[3][0] = 8;
    array[0][0] = 5;
    std::cout << array[3][0] << std::endl;
    std::cout << array[0][1] << std::endl;

    //array.shape = dimensions;
    return 0;
}