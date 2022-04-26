#include <iostream>
#include <stdio.h>
#include <Eigen/Dense>
#include <vector>

Eigen::MatrixXd A(2,2),B(2,2);

std::vector<Eigen::MatrixXd> test(){
    
    for(int i = 0; i < 2; i++){
        for( int j = 0; j < 2 ; j++){
            A(i,j) = i+j;
            B(i,j) = i-j;
        }
    }
    //std::cout << "mat A: " << std::endl << A << std::endl;
    //std::cout << "mat B: " << std::endl << B << std::endl;
    std::vector<Eigen::MatrixXd> C; // <-- D(A.rows() + B.rows(), ...)
    Eigen::MatrixXd aux1(A.rows()+B.rows(), A.cols()), aux2(1,1);
    C.push_back(aux1);
    C.push_back(aux2);
    C[0] << A,B;
    C[1] << -1;
    std::cout << "mat C: " << std::endl << C[0] << std::endl;
    std::cout << "mat C: " << std::endl << C[1] << std::endl;
    return C;
}

int main(){
    test();
    return 0;
}