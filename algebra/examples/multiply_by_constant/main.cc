#include "../../NArray.h"
#include <iostream>
#include <string>


int main(){
    //Defining shape of the N dimensional Array
    std::vector<int> dimensions;
    dimensions.push_back(2);
    dimensions.push_back(3);

    NArray array(dimensions);
    
    //Inserting values in the N dimensional Array
    array[0][0] = 1;
    array[0][1] = 2;
    array[0][2] = 3;
    array[1][0] = 4;
    array[1][1] = 5;
    array[1][2] = 6;

    //Defining a constant F
    int F = 5;

    std::cout << "The array before multiplication is: " << array << std::endl;

    array = array*F;

    std::cout << "The array after multiplication is: " << array << std::endl;
    
    return 0;
}