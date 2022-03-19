#ifndef DE_H_INCLUDED
#define DE_INCLUDED
#include <random>
#include <stdio.h>
#include <iostream>
#include<Eigen/Dense>

typedef float(*EvalFunction)(Eigen::MatrixXd);


std::random_device rseed;
std::mt19937 rng(rseed());

class DE{
public:
    DE(int N_pop, std::vector<int> ind_shape, float cr, EvalFunction evaluation, float F);

    void create_population();
    Eigen::MatrixXd generate_individual(std::vector<int> ind_shape);
    void crossover(int ind_idx);
    void mutate(int ind_idx);
    void repair(int ind_idx);
    void evolve(uint ng);
    void evaluate(int ind_idx);
    void selection(int ind_idx);
    float get_best_fit();
private:
    float cr;
    float F;
    Eigen::MatrixXd mutated_ind;
    std::vector<Eigen::MatrixXd> population;
    std::vector<float> fitness;
    EvalFunction eval;
    
};

#endif