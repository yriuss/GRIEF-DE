#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
// #include <functional>  // for std::greater



#define MINIMIZATION false
#define MAXIMIZATION false
#define N_PTS_CROSS 4

int main(){

    // bool problem_type = MINIMIZATION;

    // #if problem_type == MINIMIZATION && problem_type == MAXIMIZATION
    //     std::cout << "Error. Problem type MINIMIZATION and MAXIMIZATION macros was both defined as the same value." << std::endl;
    //     exit(EXIT_FAILURE);
    // #elif problem_type == MINIMIZATION
    //     std::cout << "Problema de Minimização" << std::endl;
    // #elif problem_type == MAXIMIZATION
    //     std::cout << "Problema de Maximização" << std::endl;
    // #endif

    // std::cout << "Ok." << std::endl;

    std::random_device rseed;
    std::mt19937 rng(rseed());
    std::uniform_real_distribution<float> r_dist(0,1);
    std::uniform_int_distribution<int> dist(0, 512 - 2);
    
    // std::vector<int> point_cross;
    int point_cross[N_PTS_CROSS];

    for(int i = 0; i < N_PTS_CROSS; i++)
        point_cross[i] = -1;


    for(int i = 0; i < N_PTS_CROSS; i++){

        point_cross[i] = dist(rng);		
        for(int j = 0; j < i; j++){
            if(point_cross[i] == point_cross[j]){
                point_cross[i] = dist(rng);
                j = 0;
            }
        }			
    }
	
    // std::sort(point_cross.begin(), point_cross.end());
    std::sort(point_cross, point_cross + sizeof point_cross/sizeof point_cross[0]);
        
    // for(int i=0; i<N_PTS_CROSS; i++)
    //     std::cout << point_cross[i] << std::endl;

    // std::cout << "asd" << teste<< std::endl;
    int a = 0;

    if(a) std::cout << "TRUE" << std::endl;
    else std::cout << "FALSE" << std::endl;

    return EXIT_SUCCESS;
}