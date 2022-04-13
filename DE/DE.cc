#include "DE.h"
#include "QS/quicksort.h"
#include <stdio.h>
#include <cmath>
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <random>


namespace DE{

	DE::DE( int N_pop, std::vector<int> ind_shape, float cr, float jr, EvalFunction evaluation, float F, 
			bool problem_type, std::vector<int> bounds, int mutation_algorithm, int crossover_algorithm): mutated_ind(ind_shape[0], ind_shape[1] ){
		
		//Initialize population
		if(N_pop > 0){
			population.reserve(N_pop);
			fitness.reserve(N_pop);
			//best_fitness.reserve(N_pop);
			eval = evaluation;
			for(int i = 0; i < N_pop; i++){
				population.emplace_back(generate_individual(ind_shape));
				fitness.emplace_back(eval(truncate_individual(ind_shape, population[i])));
			}
			this->cr = cr;
			this->F = F;
			this->jr = jr;
			this->problem_type = problem_type;
			this->mutation_algorithm = mutation_algorithm;
			this->crossover_algorithm = crossover_algorithm;
			U = bounds[1];
			L = bounds[0];
			this->N_pop = N_pop;
			this->ind_shape = ind_shape;
		}
		else {}
		
	}

	void DE::evaluate(int ind_idx){
		fitness[ind_idx] = eval(truncate_individual(ind_shape, population[ind_idx]));
	}

	Eigen::MatrixXd DE::truncate_individual(std::vector<int> ind_shape, Eigen::MatrixXd ind){
		
		Eigen::MatrixXd truncated_individual = ind;
		for (int i = 0; i < ind_shape[0]; i++){
			for (int j = 0; j < ind_shape[1]; j++)
				truncated_individual(i,j) = (int)truncated_individual(i,j);
		}

		return truncated_individual;
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

	Eigen::MatrixXd DE::generate_oppsite_individual(std::vector<int> ind_shape, int ind_idx){

		// std::cout << "Individual " << population[ind_idx] << std::endl;

		Eigen::MatrixXd opposite_individual(ind_shape[0], ind_shape[1]);

		for (int i = 0; i < ind_shape[0]; i++){
			for (int j = 0; j < ind_shape[1]; j++){
				opposite_individual(i,j) = L + U - population[ind_idx](i,j);
			}
		}

		// std::cout << "Opposite Individual " << opposite_individual << std::endl;
		// std::cout << "[ Opposite Individual Returned ]" << std::endl;

		return opposite_individual;
	}

	void DE::generate_oppsite_population(){

		// std::cout << "Generate Opposition Population Called" << std::endl;

		opposite_population.reserve(N_pop);
		opposite_fitness.reserve(N_pop);
		
		for (int i = 0; i < N_pop; i++){
			opposite_population.emplace_back(generate_oppsite_individual(ind_shape, i));
			opposite_fitness.emplace_back(eval(truncate_individual(ind_shape, opposite_population[i])));
		}

	}

	void DE::apply_opposition(){
		
		std::cout << "Apply Opposition Called" << std::endl;

		generate_oppsite_population();

		std::vector<Eigen::MatrixXd> aux_population;
		std::vector<float> aux_fitness;

		aux_population.reserve(N_pop*2);
		aux_fitness.reserve(N_pop*2);

		// std::cout << "[ Aux Population Reserved ]" << std::endl;
		
		for (int i = 0; i < N_pop; i++){
			aux_population.emplace_back(population[i]);
			aux_fitness.emplace_back(fitness[i]);
			aux_population.emplace_back(opposite_population[i]);
			aux_fitness.emplace_back(opposite_fitness[i]);
		}

		// std::cout << "[ Aux Population Created ]" << std::endl;

		std::vector<int> index_vector;
		index_vector.reserve(N_pop*2);

		for(int i = 0; i < N_pop*2; i++)
			index_vector[i] = i;

		// std::cout << "[ Index Vector Created ]" << std::endl;


		// float fitness_vector[N_pop*2];
		// for(int i = 0; i < N_pop*2; i++)
		// 	fitness[i] = aux_fitness[i];

		// std::cout << "[ Fitness Vector Created ]" << std::endl;

		
		// std::cout << "Quicksort Called" << std::endl;
		
		QS::quicksort qs;
		qs.sort( aux_fitness, index_vector, 0, (N_pop * 2) - 1 );

		// std::cout << "[ Quicksort ok ]" << std::endl;

	
		// std::cout << "[ Getting np best fitted individuals ]" << std::endl;

		#if problem_type == MINIMIZATION
			// std::cout << "MINIMIZATION" << std::endl;  
			for(int i = 0; i < N_pop; i++){
				population[i] = aux_population[index_vector[i]];
				fitness[i] = aux_fitness[i];
			}

		#elif problem_type == MAXIMIZATION
			// std::cout << "MAXIMIZATION" << std::endl;  
			for(int i = N_pop - 1; i >= 0; i--){
				population[i] = aux_population[index_vector[i]];
				fitness[i] = aux_fitness[i];
			}

		#else
			std::cout << "ERROR: Problem type was not specified. \n" << std::endl;
			exit(EXIT_FAILURE);

		#endif
		// std::cout << "[ OK ]" << std::endl;

	}

	void DE::rand_1(int ind_idx){
		
		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_int_distribution<int> dist(0, N_pop - 1);

		int idx1 = -1;
		int idx2 = -1;
		int idx3 = -1;
		
		do {
			idx1 = dist(rng);
		}
		while(idx1 == ind_idx);

		do {
			idx2 = dist(rng);
		}
		while(idx2 == ind_idx || idx2 == idx1 );

		do {
			idx3 = dist(rng);
		}
		while(idx3 == ind_idx || idx3 == idx2 || idx3 == idx1);

		mutated_ind = population[idx1] + F * (population[idx2] - population[idx3]);
		
	}

	void DE::rand_2(int ind_idx){
		
		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_int_distribution<int> dist(0, N_pop - 1);

		int idx1 = -1;
		int idx2 = -1;
		int idx3 = -1;
		int idx4 = -1;
		int idx5 = -1;

		do {
			idx1 = dist(rng);
		}
		while(idx1 == ind_idx);

		do {
			idx2 = dist(rng);
		}
		while(idx2 == ind_idx || idx2 == idx1 );

		do {
			idx3 = dist(rng);
		}
		while(idx3 == ind_idx || idx3 == idx2 || idx3 == idx1);

		do {
			idx4 = dist(rng);
		}
		while(idx4 == ind_idx || idx4 == idx3 || idx4 == idx2 || idx4 == idx1);

		do {
			idx5 = dist(rng);
		}
		while(idx5 == ind_idx || idx5 == idx4 || idx5 == idx3 || idx5 == idx2 || idx5 == idx1);
		//std::cout << "funfou";
		mutated_ind = population[idx1] + F * ( (population[idx2] - population[idx3]) + (population[idx4] - population[idx5]) );
	}

	void DE::randtobest_1(int ind_idx){

		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_int_distribution<int> dist(0, population.size() - 1);

		int idxb = -1;
		int idx1 = -1;
		int idx2 = -1;
		int idx3 = -1;
		
		idxb = get_best_idx();

		do {
			idx1 = dist(rng);
		}
		while(idx1 == ind_idx || idx1 == idxb);

		do {
			idx2 = dist(rng);
		}
		while(idx2 == ind_idx || idx2 == idxb || idx2 == idx1 );

		do {
			idx3 = dist(rng);
		}
		while(idx3 == ind_idx || idx3 == idxb || idx3 == idx2 || idx3 == idx1);

		mutated_ind = population[idx1] + F * ( (population[idxb] - population[idx1]) + (population[idx2] - population[idx3]) );
	}

	void DE::best_1(int ind_idx){

		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_int_distribution<int> dist(0, population.size() - 1);

		int idxb = -1;
		int idx2 = -1;
		int idx3 = -1;
		
		idxb = get_best_idx();

		do {
			idx2 = dist(rng);
		}
		while(idx2 == ind_idx || idx2 == idxb);

		do {
			idx3 = dist(rng);
		}
		while(idx3 == ind_idx || idx3 == idxb || idx3 == idx2);

		mutated_ind = population[idxb] + F * (population[idx2] - population[idx3]);
	}

	void DE::best_2(int ind_idx){

		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_int_distribution<int> dist(0, population.size() - 1);

		int idxb = -1;
		int idx2 = -1;
		int idx3 = -1;
		int idx4 = -1;
		int idx5 = -1;
		
		idxb = get_best_idx();

		do {
			idx2 = dist(rng);
		}
		while(idx2 == ind_idx || idx2 == idxb);

		do {
			idx3 = dist(rng);
		}
		while(idx3 == ind_idx || idx3 == idxb || idx3 == idx2);

		do {
			idx4 = dist(rng);
		}
		while(idx4 == ind_idx || idx4 == idxb || idx4 == idx3 || idx4 == idx2);

		do {
			idx5 = dist(rng);
		}
		while(idx5 == ind_idx || idx5 == idxb || idx5 == idx4 || idx5 == idx3 || idx5 == idx2);

		mutated_ind = population[idxb] + F * (population[idx2] - population[idx3]);
	}

	void DE::currenttobest_1(int ind_idx){

		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_int_distribution<int> dist(0, population.size() - 1);

		int idxb = -1;
		int idx2 = -1;
		int idx3 = -1;
		
		idxb = get_best_idx();

		do {
			idx2 = dist(rng);
		}
		while(idx2 == ind_idx || idx2 == idxb);

		do {
			idx3 = dist(rng);
		}
		while(idx3 == ind_idx || idx3 == idxb || idx3 == idx2);

		mutated_ind = population[ind_idx] + F * ((population[idxb] - population[ind_idx]) + (population[idx2] - population[idx3]));
	}

	void DE::currenttorand_1(int ind_idx){

		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_int_distribution<int> dist(0, population.size() - 1);

		int idx1 = -1;
		int idx2 = -1;
		int idx3 = -1;
		
		do {
			idx1 = dist(rng);
		}
		while(idx1 == ind_idx);

		do {
			idx2 = dist(rng);
		}
		while(idx2 == ind_idx || idx2 == idx1);

		do {
			idx3 = dist(rng);
		}
		while(idx3 == ind_idx || idx3 == idx1 || idx3 == idx2);

		mutated_ind = population[ind_idx] + F * ((population[idx1] - population[ind_idx]) + (population[idx2] - population[idx3]));
	}

	void DE::bincross(int ind_idx){
		
		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_real_distribution<float> r_dist(0,1);
		std::uniform_int_distribution<int> dist(0, ind_shape[1] - 2);
		//std::cout << "passei aqui" << std::endl;

		for(int i = 0; i < population[ind_idx].rows(); i++){
			
			float J = dist(rng);
			for(int j = 0; j < population[ind_idx].cols(); j++)
			{
				mutated_ind(i,j) = mutated_ind(i,j);	
				if(r_dist(rng) <= cr || j == J)
				{
					if(!infeasible)
						infeasible = is_infeasible(mutated_ind(i,j));
				}
				else
					mutated_ind(i,j) = population[ind_idx](i,j);
			}
		}
		
	}

	void DE::expcross(int ind_idx){
		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_real_distribution<float> r_dist(0,1);
		std::uniform_int_distribution<int> dist(0, ind_shape[1] - 2);

		Eigen::MatrixXd ind_cross = population[ind_idx];

		for(int i = 0; i < population[ind_idx].rows(); i++){
			int j = dist(rng);				
			int e = 0;
			
			while(r_dist(rng) <= cr && e < population[ind_idx].cols()){
				mutated_ind(i,j) = mutated_ind(i,j);
				if(!infeasible)
					infeasible = is_infeasible(mutated_ind(i,j));

				ind_cross(i,j) = mutated_ind(i,j);
				j = (j + 1) % (ind_shape[1]);
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
		//std::cout << "alsdjoasikl "<<crossover_algorithm << std::endl;
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

	//void DE::plot_convergence(){
	//	plt::plot(best_fitness);
	//	plt::show;
	//}

	void DE::weibull_repair(int ind_idx){
		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::weibull_distribution<double> dist(2.0,23.0);

		for(int i = 0; i < mutated_ind.rows(); i++){
			for(int j = 0; j < mutated_ind.cols(); j++){
				if(mutated_ind(i,j) > 0){
					while(mutated_ind(i,j) > U){
						mutated_ind(i,j) = dist(rng);
					}
				}else{
					while(mutated_ind(i,j) < L){
						mutated_ind(i,j) = -dist(rng);
					}
				}
			}
		}
		
	}

	uint DE::get_change_counter(){
		return change_counter;
	}

	void DE::selection(int ind_idx){
		//std::cout << mutated_ind.rows() << " "  << mutated_ind.cols() << std::endl << std::endl;
		//std::cout << mutated_ind << std::endl << std::endl;
		float mutated_fit = eval(truncate_individual(ind_shape, mutated_ind));	
		//std::cout << ind_idx << std::endl;	
		if(!problem_type){
			if(mutated_fit < fitness[ind_idx]){
				change_counter++;
				population[ind_idx] = mutated_ind;
				fitness[ind_idx] = mutated_fit;
			}
		}else{
			if(mutated_fit > fitness[ind_idx]){
				change_counter++;
				population[ind_idx] = mutated_ind;
				fitness[ind_idx] = mutated_fit;
			}
		}
	}
	//void DE::set_best_fit(){
	//	best_fitness.emplace_back(get_best_fit());
	//}
	void DE::evolve(uint ng){
		
		for(int g = 0; g < ng; g++){
			change_counter  = 0;
			for(int i = 0; i < population.size(); i++){
				mutate(i);
				crossover(i);
				if(infeasible)
					repair(i);
				selection(i);					
			}
			
			//best_fitness.emplace_back(get_best_fit());
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