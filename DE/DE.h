#ifndef DE_H_INCLUDED
#define DE_H_INCLUDED
#include <random>
#include <stdio.h>
#include <iostream>
#include<Eigen/Dense>
#include <vector>
typedef float(*EvalFunction)(Eigen::MatrixXd);
namespace DE{

#define MAXIMIZATION true
#define MINIMIZATION false

#define DEBUG 1


class DE{
public:
    DE(int N_pop, std::vector<int> ind_shape, float cr, 
    EvalFunction evaluation, float F, bool problem_type, std::vector<int> bounds);

    void create_population();
    Eigen::MatrixXd generate_individual(std::vector<int> ind_shape);
    void crossover(int ind_idx);
    void mutate(int ind_idx);
    void repair(int ind_idx);
    void evolve(uint ng);
    void evaluate(int ind_idx);
    void selection(int ind_idx);
    float get_best_fit();
    Eigen::MatrixXd get_best_ind();
    int get_max_elem();
    void get_fitness();
    bool is_infeasible(int element);
    void weibull_repair(int ind_idx);
    int get_best_idx();
    bool is_infeasible();
private:
    float cr;
    float F;
    Eigen::MatrixXd mutated_ind;
    std::vector<Eigen::MatrixXd> population;
    std::vector<float> fitness;
    EvalFunction eval;
    bool problem_type;
    int U;
    int L;
    bool infeasible = false;
};

}

#endif