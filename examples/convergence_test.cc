#include "matplotlibcpp.h"
#include "DE/DE.h"
#include <stdio.h>
#include<cmath>
#include<Eigen/Dense>
#include <vector>
#include<iostream>

namespace plt = matplotlibcpp;

void plot_convergence(std::vector<int> x,std::vector<int> y){
    plt::plot(x,y);
    plt::show();
}

float evaluation(Eigen::MatrixXd individual){
    float sum = 0;
    
    for(int i = 0; i < individual.rows() ; i++){
        for(int j = 0; j < individual.cols(); j++){
            sum += pow(individual(i,j), 1);
        }
    }
    return sum;
}

int main(){
    std::vector<int> ind_shape = {2,4}, bounds = {-24, 24};
    float cr=0.5,F=0.5;

    DE::DE de(30, 
    ind_shape, 
    cr, 
    evaluation, 
    F, 
    MAXIMIZATION, 
    bounds);

    de.evolve(150);
    
    std::cout << de.get_best_fit() << std::endl;
    std::cout << de.get_best_ind() << std::endl;
    std::cout << de.get_max_elem() << std::endl;

    

    return 0;
}