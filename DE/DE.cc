#include "DE.h"
#include <stdio.h>
#include<cmath>
#include<Eigen/Dense>
#include <vector>
#include<iostream>
namespace DE{
DE::DE(int N_pop, std::vector<int> ind_shape, float cr, EvalFunction evaluation, float F, bool problem_type, std::vector<int> bounds): mutated_ind(ind_shape[0], ind_shape[1]){
	//Initialize population
	if(N_pop > 0){
		population.reserve(N_pop);
		fitness.reserve(N_pop);

		eval = evaluation;
		for(int i = 0; i < N_pop; i++){
			population.emplace_back(generate_individual(ind_shape));
			fitness.emplace_back(eval(population[i]));
		}
		this->cr = cr;
		this->F = F;
		this->problem_type = problem_type;
		U = bounds[1];
		L = bounds[0];
	}else{

	}
}

void DE::evaluate(int ind_idx){
	fitness[ind_idx] = eval(population[ind_idx]);
}

Eigen::MatrixXd DE::generate_individual(std::vector<int> ind_shape){
	std::random_device rseed;
	std::mt19937 rng(rseed());
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
	std::random_device rseed;
	std::mt19937 rng(rseed());
	std::uniform_int_distribution<int> dist(0,population.size() - 1);
	mutated_ind = population[dist(rng)] + F*(population[dist(rng)] - population[dist(rng)]);
}

bool DE::is_infeasible(int element){
	if(element > U)
		return true;
	else
		if(element < L)
			return true;
	return false;
}

void DE::crossover(int ind_idx){
	std::random_device rseed;
	std::mt19937 rng(rseed());
	std::uniform_real_distribution<float> r_dist(0,1);
	for(int i =0; i < population[ind_idx].rows(); i++){
		for(int j = 0; j < population[ind_idx].cols(); j++){
			if(r_dist(rng) < cr){
				mutated_ind(i,j) = (int)mutated_ind(i,j);
				if(!infeasible)
					infeasible = is_infeasible(mutated_ind(i,j));
			}else
				mutated_ind(i,j) = population[ind_idx](i,j);
		}
	}
}

void DE::repair(int ind_idx){
	
	weibull_repair(ind_idx);
	infeasible = false;
}

void DE::weibull_repair(int ind_idx){
	std::random_device rseed;
	std::mt19937 rng(rseed());
	std::weibull_distribution<double> dist(2.0,23.0);
	for(int i =0; i < mutated_ind.rows(); i++){
		for(int j = 0; j < mutated_ind.cols(); j++){
			if(mutated_ind(i,j) > 0){
				while(mutated_ind(i,j) > U){
					mutated_ind(i,j) = (int)dist(rng);
				}
			}else{
				while(mutated_ind(i,j) < L){
					mutated_ind(i,j) = -(int)dist(rng);
				}
			}
		}
	}
}

void DE::selection(int ind_idx){
	float mutated_fit = eval(mutated_ind);
	
	if(!problem_type){
		if(mutated_fit < fitness[ind_idx]){
			population[ind_idx] = mutated_ind;
			fitness[ind_idx] = mutated_fit;
		}
	}else{
		if(mutated_fit > fitness[ind_idx]){
			population[ind_idx] = mutated_ind;
			fitness[ind_idx] = mutated_fit;
		}
	}

}



void DE::evolve(uint ng){
	for(int g = 0; g < ng; g++){
		for(int i = 0; i < population.size(); i++){
			mutate(i);
			crossover(i);
			if(infeasible)
				repair(i);
			selection(i);
		}
	}
}
bool DE::is_infeasible(){
	return infeasible;
}

float DE::get_best_fit(){
	if(!problem_type)
		return *std::min_element(this->fitness.begin(), this->fitness.end());
	else
		return *std::max_element(this->fitness.begin(), this->fitness.end());
}

Eigen::MatrixXd DE::get_best_ind(){
	if(!problem_type)
		return population[std::min_element(this->fitness.begin(), this->fitness.end()) - fitness.begin()];
	else
		return population[std::max_element(this->fitness.begin(), this->fitness.end()) - fitness.begin()];
}

int DE::get_best_idx(){
	if(!problem_type)
		return std::min_element(this->fitness.begin(), this->fitness.end()) - fitness.begin();
	else
		return std::max_element(this->fitness.begin(), this->fitness.end()) - fitness.begin();
}


void DE::get_fitness(){
	for(int i =0; i < 30; i++)
		std::cout << fitness[i] << std::endl;
}

int DE::get_max_elem(){
	return population[0].maxCoeff();
}
}