#include "../../NArray.h"
#include <iostream>
#include <string>


int main(){
    //Defining shape of the N dimensional Array
    std::vector<int> dimensions;
    dimensions.push_back(2);
    dimensions.push_back(3);

    NArray array(dimensions), array1(dimensions), array2(dimensions);
    
    //Inserting values in the N dimensional Arrays
    array[0][0] = 1;
    array[0][1] = 2;
    array[0][2] = 3;
    array[1][0] = 4;
    array[1][1] = 5;
    array[1][2] = 6;

    array1[0][0] = 0;
    array1[0][1] = 1;
    array1[0][2] = 2;
    array1[1][0] = 3;
    array1[1][1] = 4;
    array1[1][2] = 5;

    array2[0][0] = 1;
    array2[0][1] = 2;
    array2[0][2] = 3;
    array2[1][0] = 4;
    array2[1][1] = 5;
    array2[1][2] = 6;

    //Define constant to be multiplied
    int F = 5;

    //Apply sum
    array = array1 + array2 * F;

    std::cout << "The array 1: " << array1 << std::endl;

    std::cout << "The array 2: " << array2 << std::endl;

    std::cout << "The resultant array (array1 + array2 * constant): " << array << std::endl;
    
    return 0;
}