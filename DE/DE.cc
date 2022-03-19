#include "matplotlibcpp.h"
#include "DE.h"

#include<cmath>

#include <vector>


namespace plt = matplotlibcpp;



DE::DE(int N_pop, std::vector<int> ind_shape, float cr, EvalFunction evaluation, float F): mutated_ind(ind_shape[0], ind_shape[1]){
    //Initialize population
    population.reserve(N_pop);
    fitness.reserve(N_pop);

    eval = evaluation;
    for(int i = 0; i < N_pop; i++){
        population.emplace_back(generate_individual(ind_shape));
        fitness[i] = eval(population[i]);
    }
    this->cr = cr;
    this->F = F;
}

void DE::evaluate(int ind_idx){
    fitness[ind_idx] = eval(population[ind_idx]);
}

Eigen::MatrixXd DE::generate_individual(std::vector<int> ind_shape){
    
    std::uniform_int_distribution<int> dist(-24,24);
    std::uniform_real_distribution<float> distr(0,1);
    Eigen::MatrixXd individual(ind_shape[0], ind_shape[1]);

    //std::cout << distr(rng) << std::endl;
    for(int i = 0; i < ind_shape[0]; i++){
        for(int j = 0; j < ind_shape[1]; j++){
            individual(i,j) = dist(rng);
        }
    }

    //std::cout << individual << std::endl;
    //std::cout << (float)(individual[0][0]) << std::endl;
    return individual;
}

void DE::mutate(int ind_idx){
    std::uniform_int_distribution<int> dist(0,population.size() - 1);
    mutated_ind = population[ind_idx] + F*(population[dist(rng)] - population[ind_idx]);
}

void DE::crossover(int ind_idx){
    std::uniform_real_distribution<float> r_dist(0,1);
    for(int i =0; i < 2; i++){
        for(int j = 0; j < 4; j++){
            if(r_dist(rng) < cr)
                mutated_ind(i,j) = (int)mutated_ind(i,j);
            else
                mutated_ind(i,j) = population[ind_idx](i,j);
        }
    }
}

void DE::repair(int ind_idx){

}

void DE::selection(int ind_idx){
    float mutated_fit = eval(mutated_ind);
    
    if(mutated_fit < fitness[ind_idx]){
        population[ind_idx] = mutated_ind;
        fitness[ind_idx] = mutated_fit;
        
    }

}


void DE::evolve(uint ng){
    for(int g = 0; g < ng; g++){
        for(int i = 0; i < population.size(); i++){
            mutate(i);
            repair(i);
            crossover(i);
            selection(i);
        }
    }
}

float DE::get_best_fit(){
    auto result = std::min_element(fitness.begin(), fitness.end());
    return *result;
}

void plot_convergence(std::vector<int> x,std::vector<int> y){
    plt::plot(x,y);
    plt::show();
}



float evaluation(Eigen::MatrixXd individual){
    float sum = 0;
    
    for(int i = 0; i < 2 ; i++){
        for(int j = 0; j < 4; j++){
            sum += pow(individual(i,j), 2);
        }
    }
    return sum;
}

int main(){
    std::vector<int> ind_shape = {2,4};

    float cr=0.5,F=0.5;

    DE de(300, ind_shape, cr, evaluation, F);

    de.evolve(100);
    
    std::cout << de.get_best_fit() << std::endl;

    //std::random_device rseed;
    //std::mt19937 rng(rseed());
    //std::uniform_int_distribution<int> dist(-24,24);
    //
    //Eigen::MatrixXd individual(ind_shape);
    //
//
    //for(int i = 0; i < ind_shape[0]; i++){
    //    for(int j = 0; j < ind_shape[1]; j++){
    //        //std::cout << i << " " << j << std::endl;
    //        individual(i,j) = dist(rng);
    //    }
    //}
    
    //std::cout << evaluation(individual);


    //plot_convergence({1,2,3,4},{3,2,4,3});

    return 0;
}