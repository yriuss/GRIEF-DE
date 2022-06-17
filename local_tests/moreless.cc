#include <iostream>
#include <stdio.h>
#include <Eigen/Dense>


int main()
{
    Eigen::MatrixXd x(4,3);
    x << 1,2,3,4,
         5,6,7,8,
         9,10,11,12;

    Eigen::MatrixXd a, b;
    a = x;
    b = x;
    
    for (int i = 0; i < a.rows(); i++){
        for(int j = 0; j < a.cols(); j++){
            a(i,j) = a(i,j) + 1;
            b(i,j) = b(i,j) - 1;
        }
    }


    for (int i = 0; i < a.rows(); i++){
        for(int j = 0; j < a.cols(); j++)
            std::cout << " " << a(i,j) << " ";
        std::cout << std::endl;
    }
    
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    
    for (int i = 0; i < a.rows(); i++){
        for(int j = 0; j < a.cols(); j++)
            std::cout << " " << b(i,j) << " ";
        std::cout << std::endl;
    }

    return 0;
}