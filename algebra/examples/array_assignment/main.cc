#include "../../NArray.h"
#include <iostream>
#include <string>


int main(){
    //Defining shape of the N dimensional Arrays
    std::vector<int> dimensions = {4,3};
    std::vector<int> dimensions2 = {3};

    NArray array1(dimensions), array2(dimensions2);
    

    array1[0][0] = 0;
    array1[0][1] = 1;
    array1[0][2] = 2;
    array1[1][0] = 3;
    array1[1][1] = 4;
    array1[1][2] = 5;
    array1[2][0] = 6;
    array1[2][1] = 7;
    array1[2][2] = 8;
    array1[3][0] = 9;
    array1[3][1] = 10;
    array1[3][2] = 11;

    array2[0] = 1;
    array2[1] = 2;
    array2[2] = 3;

    std::cout << "The array 1: " << array1 << std::endl;

    std::cout << "The array 2: " << array2 << std::endl;

    array1[2] = array2;

    

    std::cout << "The array 1 after assignment: " << array1 << std::endl;
    
    return 0;
}