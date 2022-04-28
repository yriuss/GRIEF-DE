
#include <stdio.h>
#include <iostream>
#include <Eigen/Dense>
#include <vector>

int main(){
    
    Eigen::MatrixXd m1(4,3);
    Eigen::MatrixXd m2(4,4);

    m1 << 1, 3, 5, 8,
         2, 4, 6, 8,
         8, 4, 6, 8;
    
    m2 << 1, 3, 5, 8,
         2, 4, 6, 8,
         8, 4, 6, 8,
         8, 4, 6, 8;
    
    if(m1 == m2)
        std::cout << "The matrixes are equal!" << std::endl;
    else
        std::cout << "The matrixes are not equal!" << std::endl;


    return 0;
}