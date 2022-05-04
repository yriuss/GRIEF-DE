#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#define MINIMIZATION false
#define MAXIMIZATION false

int main(){

    bool problem_type = MINIMIZATION;

    #if problem_type == MINIMIZATION && problem_type == MAXIMIZATION
        std::cout << "Error. Problem type MINIMIZATION and MAXIMIZATION macros was both defined as the same value." << std::endl;
        exit(EXIT_FAILURE);
    #elif problem_type == MINIMIZATION
        std::cout << "Problema de Minimização" << std::endl;
    #elif problem_type == MAXIMIZATION
        std::cout << "Problema de Maximização" << std::endl;
    #endif

    std::cout << "Ok." << std::endl;


    return EXIT_SUCCESS;
}