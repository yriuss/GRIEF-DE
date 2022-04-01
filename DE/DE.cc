#include "DE.h"
#include <stdio.h>
#include <cmath>
#include <Eigen/Dense>
#include <vector>
#include <iostream>

namespace DE{

	DE::DE( int N_pop, std::vector<int> ind_shape, float cr, EvalFunction evaluation, float F, 
			bool problem_type, std::vector<int> bounds): mutated_ind(ind_shape[0], ind_shape[1], int mutation_algorithm ){
		
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
			this->mutation_algorithm = mutation_algorithm;
			U = bounds[1];
			L = bounds[0];
			N_pop = N_pop;
		}
		else {}
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

	void DE::rand_1(int ind_idx){

		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_int_distribution<int> dist(0, N_pop - 1);

		int idx1 = NULL;
		int idx2 = NULL;
		int idx3 = NULL;
		
		do {
			idx1 = dist(rng);
		}
		while(idx1 != ind_idx);

		do {
			idx2 = dist(rng);
		}
		while(idx2 != ind_idx && idx2 != idx1 );

		do {
			idx3 = dist(rng);
		}
		while(idx3 != ind_idx && idx3 != idx2 && idx3 != idx1);

		Eigen::MatrixXd ind_idx1 = population[idx1]; 
		Eigen::MatrixXd ind_idx2 = population[idx2]; 
		Eigen::MatrixXd ind_idx3 = population[idx3]; 

		Eigen::MatrixXd ind_u = ind_idx1 + F * (ind_idx2 - ind_idx3);

		mutated_ind = ind_u;

	}

	void DE::rand_2(int ind_idx){

		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_int_distribution<int> dist(0, N_pop - 1);

		int idx1 = NULL;
		int idx2 = NULL;
		int idx3 = NULL;
		int idx4 = NULL;
		int idx5 = NULL;

		do {
			idx1 = dist(rng);
		}
		while(idx1 != ind_idx);

		do {
			idx2 = dist(rng);
		}
		while(idx2 != ind_idx && idx2 != idx1 );

		do {
			idx3 = dist(rng);
		}
		while(idx3 != ind_idx && idx3 != idx2 && idx3 != idx1);

		do {
			idx4 = dist(rng);
		}
		while(idx4 != ind_idx && idx4 != idx3 && idx4 != idx2 && idx4 != idx1);

		do {
			idx5 = dist(rng);
		}
		while(idx5 != ind_idx && idx5 != idx4 && idx5 != idx3 && idx5 != idx2 && idx5 != idx1);


		Eigen::MatrixXd ind_idx1 = population[idx1]; 
		Eigen::MatrixXd ind_idx2 = population[idx2]; 
		Eigen::MatrixXd ind_idx3 = population[idx3]; 
		Eigen::MatrixXd ind_idx4 = population[idx4]; 
		Eigen::MatrixXd ind_idx5 = population[idx5]; 

		Eigen::MatrixXd	ind_u = ind_idx1 + F * ( (idx2 - idx3) + (idx4 - idx5) );

		mutated_ind = ind_u;
	}

	void DE::randtobest_1(int ind_idx){

		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_int_distribution<int> dist(0, N_pop - 1);

		int idxb = NULL;
		int idx1 = NULL;
		int idx2 = NULL;
		int idx3 = NULL;
		
		idxb = get_best_idx();

		do {
			idx1 = dist(rng);
		}
		while(idx1 != ind_idx && idx1 != idxb);

		do {
			idx2 = dist(rng);
		}
		while(idx2 != ind_idx && idx2 != idxb && idx2 != idx1 );

		do {
			idx3 = dist(rng);
		}
		while(idx3 != ind_idx && idx3 != idxb && idx3 != idx2 && idx3 != idx1);

		Eigen::MatrixXd ind_idxb = population[idxb];
		Eigen::MatrixXd ind_idx1 = population[idx1]; 
		Eigen::MatrixXd ind_idx2 = population[idx2]; 
		Eigen::MatrixXd ind_idx3 = population[idx3]; 

		Eigen::MatrixXd ind_u = ind_idx1 + F * ( (ind_idxb - ind_idx1) + (ind_idx2 - ind_idx3) );

		mutated_ind = ind_u;
	}

	void DE::best_1(int ind_idx){

		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_int_distribution<int> dist(0, N_pop - 1);

		int idxb = NULL;
		int idx2 = NULL;
		int idx3 = NULL;
		
		idxb = get_best_idx();

		do {
			idx2 = dist(rng);
		}
		while(idx2 != ind_idx && idx2 != idxb);

		do {
			idx3 = dist(rng);
		}
		while(idx3 != ind_idx && idx3 != idxb && idx3 != idx2);

		Eigen::MatrixXd ind_idxb = population[idxb];
		Eigen::MatrixXd ind_idx2 = population[idx2]; 
		Eigen::MatrixXd ind_idx3 = population[idx3]; 

		Eigen::MatrixXd ind_u = ind_idxb + F * (ind_idx2 - ind_idx3);

		mutated_ind = ind_u;

	}

	void DE::best_2(int ind_idx){

		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_int_distribution<int> dist(0, N_pop - 1);

		int idxb = NULL;
		int idx2 = NULL;
		int idx3 = NULL;
		int idx4 = NULL;
		int idx5 = NULL;
		
		idxb = get_best_idx();

		do {
			idx2 = dist(rng);
		}
		while(idx2 != ind_idx && idx2 != idxb);

		do {
			idx3 = dist(rng);
		}
		while(idx3 != ind_idx && idx3 != idxb && idx3 != idx2);

		do {
			idx4 = dist(rng);
		}
		while(idx4 != ind_idx && idx4 != idxb && idx4 != idx3 && idx4 != idx2);

		do {
			idx5 = dist(rng);
		}
		while(idx5 != ind_idx && idx5 != idxb && idx5 != idx4 && idx5 != idx3 && idx5 != idx2);
		
		Eigen::MatrixXd ind_idxb = population[idxb];
		Eigen::MatrixXd ind_idx2 = population[idx2]; 
		Eigen::MatrixXd ind_idx3 = population[idx3]; 
		Eigen::MatrixXd ind_idx4 = population[idx4]; 
		Eigen::MatrixXd ind_idx5 = population[idx5]; 

		Eigen::MatrixXd ind_u = ind_idxb + F * (ind_idx2 - ind_idx3);

		mutated_ind = ind_u;
	}

	void DE::currenttobest_1(int ind_idx){

		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_int_distribution<int> dist(0, N_pop - 1);

		int idxb = NULL;
		int idx2 = NULL;
		int idx3 = NULL;
		
		idxb = get_best_idx();

		do {
			idx2 = dist(rng);
		}
		while(idx2 != ind_idx && idx2 != idxb);

		do {
			idx3 = dist(rng);
		}
		while(idx3 != ind_idx && idx3 != idxb && idx3 != idx2);

		Eigen::MatrixXd ind_idxi = population[ind_idx];
		Eigen::MatrixXd ind_idxb = population[idxb];
		Eigen::MatrixXd ind_idx2 = population[idx2]; 
		Eigen::MatrixXd ind_idx3 = population[idx3]; 

		Eigen::MatrixXd ind_u = ind_idxi + F * ((ind_idxb - ind_idxi) + (ind_idx2 - ind_idx3));

		mutated_ind = ind_u;
	}

	void DE::currenttorand_1(int ind_idx){

		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_int_distribution<int> dist(0, N_pop - 1);

		int idx1 = NULL;
		int idx2 = NULL;
		int idx3 = NULL;
		
		do {
			idx1 = dist(rng);
		}
		while(idx1 != ind_idx);

		do {
			idx2 = dist(rng);
		}
		while(idx2 != ind_idx && idx2 != idx1);

		do {
			idx3 = dist(rng);
		}
		while(idx3 != ind_idx && idx3 != idx1 && idx3 != idx2);

		Eigen::MatrixXd ind_idxi = population[ind_idx];
		Eigen::MatrixXd ind_idx1 = population[idx1];
		Eigen::MatrixXd ind_idx2 = population[idx2]; 
		Eigen::MatrixXd ind_idx3 = population[idx3]; 

		Eigen::MatrixXd ind_u = ind_idxi + F * ((ind_idx1 - ind_idxi) + (ind_idx2 - ind_idx3));

		mutated_ind = ind_u;
	}

	void DE::bincross(int ind_idx){

		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_real_distribution<float> r_dist(0,1);
		std::uniform_int_distribution<int> dist(0, ind_shape[1] - 1);

		for(int i = 0; i < population[ind_idx].rows(); i++){
			
			float J = dist(rng);
			for(int j = 0; j < population[ind_idx].cols(); j++){			
				if(r_dist(rng) <= cr || j == J){
					mutated_ind(i,j) = (int)mutated_ind(i,j);
					if(!infeasible)
						infeasible = is_infeasible(mutated_ind(i,j));
				}else
					mutated_ind(i,j) = population[ind_idx](i,j);
			}
		}
	}

	void DE::expcross(int ind_idx){
		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_real_distribution<float> r_dist(0,1);
		std::uniform_int_distribution<int> dist(0, ind_shape[1] - 1);

		Eigen::MatrixXd ind_cross = population[ind_idx];

		for(int i = 0; i < population[ind_idx].rows(); i++){
			int j = dist(rng);				
			int e = 0;
			
			while(r_dist(rng) <= cr && e < population[ind_idx].cols()){
				ind_cross(i,j) = mutated_ind(i,j);
				j = (j + 1) % (ind_shape[1] - 1);
				e++;
			}			
		}	
		mutated_ind = ind_cross;			
	}

	void DE::mutate(int ind_idx){
		
		switch(mutation_algorithm){

			case 0:
				rand_1(ind_idx); break;
			case 1:
				rand_2(ind_idx); break;
			case 2:
				randtobest_1(ind_idx); break;	
			case 3:
				best_1(ind_idx); break;
			case 4:
				best_2(ind_idx); break;
			case 5:
				currenttobest_1(ind_idx); break;
			case 6:
				currenttorand_1(ind_idx); break;
		}		

	}

	void DE::crossover(int ind_idx){

		switch(crossover_algorithm){
			case 0:
				bincross(ind_idx); break;
			case 1:
				expcross(ind_idx); break;
		}

	}

	// void DE::mutate(int ind_idx){
	// 	std::random_device rseed;
	// 	std::mt19937 rng(rseed());
	// 	std::uniform_int_distribution<int> dist(0,population.size() - 1);
	// 	mutated_ind = population[dist(rng)] + F*(population[dist(rng)] - population[dist(rng)]);
	// }

	bool DE::is_infeasible(int element){
		// if(element > U)
		// 	return true;
		// else
		// 	if(element < L)
		// 		return true;

		if (element > U || element < L)
			return true;
		return false;
	}

	// void DE::crossover(int ind_idx){
	// 	std::random_device rseed;
	// 	std::mt19937 rng(rseed());
	// 	std::uniform_real_distribution<float> r_dist(0,1);
	// 	for(int i = 0; i < population[ind_idx].rows(); i++){
	// 		for(int j = 0; j < population[ind_idx].cols(); j++){
	// 			if(r_dist(rng) < cr){
	// 				mutated_ind(i,j) = (int)mutated_ind(i,j);
	// 				if(!infeasible)
	// 					infeasible = is_infeasible(mutated_ind(i,j));
	// 			}else
	// 				mutated_ind(i,j) = population[ind_idx](i,j);
	// 		}
	// 	}
	// }

	void DE::repair(int ind_idx){
		
		weibull_repair(ind_idx);
		infeasible = false;
	}

	void DE::weibull_repair(int ind_idx){
		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::weibull_distribution<double> dist(2.0,23.0);
		for(int i = 0; i < mutated_ind.rows(); i++){
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
		for(int i = 0; i < 30; i++)
			std::cout << fitness[i] << std::endl;
	}

	int DE::get_max_elem(){
		return population[0].maxCoeff();
	}
}