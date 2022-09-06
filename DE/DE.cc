#include "DE.h"



namespace DE{
	void DE::reduce_mut(){
		//if(F_mut > 0)
		//	F_mut -= 0.008;
		//std::cout << F_mut << std::endl;
	}
	std::vector<std::vector<double>> DE::get_F(){
		return F;
	}
	std::vector<int> DE::sort_idxs(std::vector<double> v){
		std::vector<int> index(v.size(), 0);
		for (int i = 0 ; i != index.size() ; i++) {
		    index[i] = i;
		}
		std::sort(index.begin(), index.end(),
		    [&](const int& a, const int& b) {
		        return (v[a] < v[b]);
		    }
		);
		return index;
	}
	
	std::vector<Eigen::MatrixXd> DE::pop(){
		return population;
	}

	void DE::reset(){
		repair_counter = 0;
		
		mut_counter = 0;
		crossed_counter = 0;
		this->Measurements::reset();
		count_cross1 = 0;
		count_mut1 = 0;
#if SECOND_MUTATED_FIT
		count_cross2 = 0;
		count_mut2 = 0;
#endif
#if CURRENT_TO_RAND||RAND_TO_BEST_MOD
		//Initialize population
		if(N_pop > 0){
			population.clear();
			fitness.clear();
			fitness_aux.clear();
			F.clear();

			population.reserve(N_pop);
			F.reserve(N_pop);
			fitness.reserve(N_pop);
			fitness_aux.reserve(N_pop);
			//best_fitness.reserve(N_pop);
			for(int i = 0; i < N_pop; i++){
				population.emplace_back(generate_individual(ind_shape));
#if MEAN_WORST
				std::vector<double> F1;
#endif
				#if ROUND_ON_MUTATION
					
					std::vector<double> vaux = eval(population[i]), vaux2;
					for(int i = 0; i < vaux.size() - 1; i++)
						vaux2.push_back(vaux[i]);
					this->F.emplace_back(vaux2);
					actual_error = vaux[vaux.size()-1];
				#else
					this->F.emplace_back(eval(truncate_individual(ind_shape, population[i])));
				#endif
#if MEAN_WORST
				F1 = this->F[i];
				std::sort(F1.begin(), F1.end(), std::greater<int>());
#endif
				//std::cout << this->F[i](i,511);
				
				float fitness = 0;
				int min = 0;
				for(int j = 0; j < ind_shape[0]; j++){
					if(min > this->F[i][j])
						min = this->F[i][j];
				}

				this->fitness_aux.emplace_back(this->F[i][0]);

#if MEAN_WORST
				for(int j = ind_shape[0] - K; j < ind_shape[0]; j++){
					fitness+= F1[j];
				}
				fitness /= K;
				for(int j = 0; j < ind_shape[0]; j++){
					this->F[i][j] = this->F[i][j] / (min);
				}
				fitness /= K;
#else
				for(int j = 0; j < ind_shape[0]; j++){
					fitness+= this->F[i][j];
					this->F[i][j] = this->F[i][j] / (min);
					//std::cout << this->F[i](j,j);
				}
				fitness /= 512;
#endif
				
				
				std::vector<int> idxs = sort_idxs(this->F[i]);		
				
				
				this->fitness.emplace_back(fitness);
				
			}
			
			
		}
#if PLUS_BEST
			population.push_back(generate_individual(ind_shape));
			

			std::vector<double> vaux = eval(population[population.size() - 1]), vaux2;
			for(int i = 0; i < vaux.size() - 1; i++)
				vaux2.push_back(vaux[i]);
			this->F.push_back(vaux2);
			actual_error = vaux[vaux.size()-1];

			float fitness = 0;
			int min = 0;
			for(int j = 0; j < ind_shape[0]; j++){
				if(min > this->F[population.size() - 1][j])
					min = this->F[population.size() - 1][j];
			}

			this->fitness_aux.push_back(this->F[population.size() - 1][0]);

			for(int j = 0; j < ind_shape[0]; j++){
				fitness+= this->F[population.size() - 1][j];
				this->F[population.size() - 1][j] = this->F[population.size() - 1][j] / (min);
				//std::cout << this->F[i](j,j);
			}
			fitness /= 512;

			this->fitness.push_back(fitness);

#endif
#else
		//Initialize population
		if(N_pop > 0){
			population.reserve(N_pop);
			fitness.reserve(N_pop);
			//best_fitness.reserve(N_pop);
			eval = evaluation;
#if READ_BEST_IND
		read_individuals(30);
		//std::cout << population[0];
		for(int i = 0; i < N_pop; i++){
			fitness.emplace_back(eval(truncate_individual(ind_shape, population[i])));
		}
#else
			for(int i = 0; i < N_pop; i++){
				population.emplace_back(generate_individual(ind_shape));

				#if ROUND_ON_MUTATION
					fitness.emplace_back(eval(population[i]));
				#else
					fitness.emplace_back(eval(truncate_individual(ind_shape, population[i])));
				#endif

			}
#endif
		}
#endif
	}
	void DE::reduce_cr(){
		cr -= 0.68/2000;
		std::cout << cr << std::endl;
	}
	DE::DE( int N_pop, std::vector<int> ind_shape, float cr, float jr, EvalFunction evaluation, float F, 
			bool problem_type, std::vector<int> bounds, int mutation_algorithm, int crossover_algorithm, int K, int sel_type, int worsts, bool cr_reduction): mutated_ind(ind_shape[0], ind_shape[1] ){
		this->cr_reduction = cr_reduction;
		for(int i = 0; i < K; i ++)
			worst_idxs.push_back(0);
		for(int i = 0; i < 512-K; i ++)
			not_worst_idxs.push_back(0);
		count_cross1 = 0;
		count_mut1 = 0;
		this->worsts = worsts;
#if SECOND_MUTATED_FIT
		count_cross2 = 0;
		count_mut2 = 0;
#endif
#if CURRENT_TO_RAND||RAND_TO_BEST_MOD
		//Initialize population
		if(N_pop > 0){
			this->K = K;
			selection_type = sel_type;
			//population.reserve(N_pop);
			this->F.reserve(N_pop);
			fitness.reserve(N_pop);
			fitness_aux.reserve(N_pop);
			//best_fitness.reserve(N_pop);
			eval = evaluation;
			
			for(int i = 0; i < N_pop; i++){
				population.push_back(generate_individual(ind_shape));
			
#if MEAN_WORST
				std::vector<double> F1;
#endif
				#if ROUND_ON_MUTATION
					std::vector<double> vaux = eval(population[i]), vaux2;
					
					for(int i = 0; i < vaux.size() - 1; i++)
						vaux2.push_back(vaux[i]);
					
					this->F.emplace_back(vaux2);
					actual_error = vaux[vaux.size()-1];
					
				#else
					this->F.emplace_back(eval(truncate_individual(ind_shape, population[i])));
				#endif
#if MEAN_WORST
				F1 = this->F[i];
				std::sort(F1.begin(), F1.end(), std::greater<int>());
#endif
				//std::cout << this->F[i](i,511);
				
				float fitness = 0;
				int min = 0;
				for(int j = 0; j < ind_shape[0]; j++){
					if(min > this->F[i][j])
						min = this->F[i][j];
				}
				
				this->fitness_aux.emplace_back(this->F[i][0]);
				

#if MEAN_WORST
				for(int j = ind_shape[0] - K; j < ind_shape[0]; j++){
					fitness+= F1[j];
				}
				fitness /= K;
				for(int j = 0; j < ind_shape[0]; j++){
					this->F[i][j] = this->F[i][j] / (min);
				}
				fitness /= K;
#else
				for(int j = 0; j < ind_shape[0]; j++){
					fitness+= this->F[i][j];
					this->F[i][j] = this->F[i][j] / (min);
					//std::cout << this->F[i](j,j);
				}
				fitness /= 512;
#endif
				
				std::vector<int> idxs = sort_idxs(this->F[i]);		
				
				
				this->fitness.emplace_back(fitness);
				
			}

#if PLUS_BEST
			population.push_back(generate_individual(ind_shape));
						

			std::vector<double> vaux = eval(population[population.size() - 1]), vaux2;
			for(int i = 0; i < vaux.size() - 1; i++)
				vaux2.push_back(vaux[i]);
			this->F.push_back(vaux2);
			actual_error = vaux[vaux.size()-1];


			float fitness = 0;
			int min = 0;
			for(int j = 0; j < ind_shape[0]; j++){
				if(min > this->F[population.size() - 1][j])
					min = this->F[population.size() - 1][j];
			}

			this->fitness_aux.push_back(this->F[population.size() - 1][0]);

			for(int j = 0; j < ind_shape[0]; j++){
				fitness+= this->F[population.size() - 1][j];
				this->F[population.size() - 1][j] = this->F[population.size() - 1][j] / (min);
				//std::cout << this->F[i](j,j);
			}
			fitness /= 512;

			this->fitness.push_back(fitness);

#endif

			this->cr = cr;
			if(cr_reduction)
				this->cr=0.98;
			this->jr = jr;
			this->F_mut = F;
			this->problem_type = problem_type;
			this->mutation_algorithm = mutation_algorithm;
			this->crossover_algorithm = crossover_algorithm;
			U = bounds[1];
			L = bounds[0];
			this->N_pop = N_pop;
			this->ind_shape = ind_shape;
			
		}
#else
		//Initialize population
		if(N_pop > 0){
			population.reserve(N_pop);
			fitness.reserve(N_pop);
			//best_fitness.reserve(N_pop);
			eval = evaluation;
#if READ_BEST_IND
		read_individuals(30);
		//std::cout << population[0];
		for(int i = 0; i < N_pop; i++){
			fitness.emplace_back(eval(truncate_individual(ind_shape, population[i])));
		}
#else
			for(int i = 0; i < N_pop; i++){
				population.emplace_back(generate_individual(ind_shape));

				#if ROUND_ON_MUTATION
					fitness.emplace_back(eval(population[i]));
				#else
					fitness.emplace_back(eval(truncate_individual(ind_shape, population[i])));
				#endif

			}
#endif
			
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
			
			check_duplicates();
		}
#endif
		
	}


	void DE::evaluate(int ind_idx){
#if CURRENT_TO_RAND||RAND_TO_BEST_MOD
#else
		fitness[ind_idx] = eval(truncate_individual(ind_shape, population[ind_idx]));
#endif
	}

	void DE::read_individuals(int n_of_individuals){
		std::string CURRENT_DIR = get_current_dir_name();
		using namespace std;
		for(int i = 0; i < n_of_individuals; i++){
			//std::cout << CURRENT_DIR +"/../best_inds/" + std::to_string(i) + ".txt";exit(-1);
			ifstream file(CURRENT_DIR +"/../best_inds/" + std::to_string(i) + ".txt");
			Eigen::MatrixXd mat(512,4);
			std::string line;
			uint16_t k = 0, j = 0;
			bool successful=false;
			std::string cell;
			//std::cout << CURRENT_DIR +"../GRIEF_CUDA/" + fileName;
			while (std::getline(file, line)) {
				//std::cout << line;
				std::vector<int> v;
				istringstream is(line);
				while (std::getline(is, cell, ' ')) {
					//std::cout << std::stoi(cell);
					if(j < 4){
						///std::cout << k << std::endl;
						mat(k,j) = std::stoi(cell);
						
					}
					j++;
				}
				successful=true;
				k++;
				j = 0;
			}
			population.emplace_back(mat);
			file.close();
		}
	}

	Eigen::MatrixXd DE::truncate_individual(std::vector<int> ind_shape, Eigen::MatrixXd ind){
		Eigen::MatrixXd truncated_individual = ind;
		for (int i = 0; i < ind_shape[0]; i++){
			for (int j = 0; j < ind_shape[1]; j++){
				truncated_individual(i,j) = round(truncated_individual(i,j));
			}
		}
		
		return truncated_individual;
	}



	Eigen::MatrixXd DE::generate_individual(std::vector<int> ind_shape){
		
		std::random_device rseed;
		std::mt19937 rng(rseed());
		//std::uniform_int_distribution<int> dist(-24,24);
		std::uniform_real_distribution<float> distr(-24,24);
		Eigen::MatrixXd individual(ind_shape[0], ind_shape[1]);

		//std::cout << distr(rng) << std::endl;
		for(int i = 0; i < ind_shape[0]; i++){
			for(int j = 0; j < ind_shape[1]; j++){
				individual(i,j) = distr(rng);
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
#if CURRENT_TO_RAND||RAND_TO_BEST_MOD
#else
			opposite_fitness.emplace_back(eval(truncate_individual(ind_shape, opposite_population[i])));
#endif
		}

	}

#if CURRENT_TO_RAND
	
	int DE::distance(Eigen::MatrixXd m1, Eigen::MatrixXd m2){
		int dist = 0;
		for(int i = 0; i < 512; i++){
			for(int j=0; j < 4; j++){
				dist += pow(m1(i,j),2) - pow(m2(i,j),2);
			}
		}
		return sqrt(dist);
	}

	void DE::currenttorand_modified(int ind_idx){
		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_int_distribution<int> dist(0, population.size() - 1);
		std::uniform_real_distribution<float> r_dist(0,1);
		
		Eigen::MatrixXd ind1 = generate_individual(ind_shape), ind2 = generate_individual(ind_shape),ind3 = generate_individual(ind_shape), ind4 = generate_individual(ind_shape);
		
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
		

		Eigen::MatrixXd F(512,512);
		F.setZero(512,512);
		float th = 1./3;
		Eigen::MatrixXd aux = population[ind_idx];
		std::vector<int> idxs = sort_idxs(this->F[ind_idx]);
		for(int i = 0; i < K; i++){
			worst_idxs[i] = idxs[511 - i];
		}

		for(int i = 0; i < 512 - K; i++){
			not_worst_idxs[i] = idxs[i];
		}

		if(!worsts){
			for(int j = 0; j < ind_shape[0]; j++){
				if(j < 512-K){
					for(int i = 0; i < ind_shape[1]; i++){
						ind3(idxs[j],i) = 0;
					}
				}
				else{
					for(int i = 0; i < ind_shape[1]; i++){
						aux(idxs[j],i) = 0;
						
					}
	
					F(idxs[j],idxs[j]) = this->F[ind_idx][idxs[j]];
				}
			}
		}else{
			for(int j = 0; j < ind_shape[0]; j++){
				if(j < 512-K){
					if(j < 512-K - 10 && j >= 512-K - 20){
						if(r_dist(rng) < (float)th){
							for(int i = 0; i < ind_shape[1]; i++){
								aux(idxs[j],i) = 0;
						
							}
							F(idxs[j],idxs[j]) = this->F[ind_idx][idxs[j]];
						}
						else{
							for(int i = 0; i < ind_shape[1]; i++){
								ind3(idxs[j],i) = 0;
							}
							F(idxs[j],idxs[j]) = 0;
						}
					}else{
						if(j < 512-K - 20 && j >= 512-K - 30){
							if(r_dist(rng) < (float)th/3){
								for(int i = 0; i < ind_shape[1]; i++){
									aux(idxs[j],i) = 0;
						
								}
								F(idxs[j],idxs[j]) = this->F[ind_idx][idxs[j]];
							}else{
								for(int i = 0; i < ind_shape[1]; i++){
									ind3(idxs[j],i) = 0;
								}		
								F(idxs[j],idxs[j]) = 0;
							}
						}else{
							for(int i = 0; i < ind_shape[1]; i++){
								ind3(idxs[j],i) = 0;
							}
							F(idxs[j],idxs[j]) = 0;
						}
					}

				}
				else{
					for(int i = 0; i < ind_shape[1]; i++){
						aux(idxs[j],i) = 0;
						
					}
					F(idxs[j],idxs[j]) = this->F[ind_idx][idxs[j]];
				}
			}
		}
		
		//exit(-1);
		switch (mutation_algorithm)
		{
		case 0:
			mutated_ind = aux + F * (ind1 - ind2);
			break;
		case 1:
			mutated_ind = aux + ind3 + F * (ind1 - ind2);
			break;
		default:
			break;
		}
		for(int i = 0; i < mutated_ind.rows(); i++){
			for(int j = 0; j < mutated_ind.cols(); j++){
				if(!infeasible)
					infeasible = is_infeasible(mutated_ind(i,j));
			}
		}
		
		//mutated_ind = population[ind_idx] + F * ((population[idx1] - population[ind_idx]) + (population[idx2] - population[idx3]));
	}
#if SECOND_MUTATED_FIT

	void DE::currenttorand_modified2(int ind_idx){

		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_int_distribution<int> dist(0, population.size() - 1);
		std::uniform_real_distribution<float> r_dist(0,1);
		Eigen::MatrixXd ind1 = generate_individual(ind_shape), ind2 = generate_individual(ind_shape),ind3 = generate_individual(ind_shape), ind4 = generate_individual(ind_shape);

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
		Eigen::MatrixXd F(512,512);
		F.setZero(512,512);
		float th = 1./3;

		std::vector<int> idxs = sort_idxs(this->F[ind_idx]);
		for(int j = 0; j < ind_shape[0]; j++){
			if(j < 512-K){
				if(j < 512-K - 10 && j >= 512-K - 20){
					if(r_dist(rng) < th){
						F(idxs[j],idxs[j]) = this->F[ind_idx][idxs[j]];
					}
					else{
						F(idxs[j],idxs[j]) = 0;
					}
				}else{
					if(j < 512-K - 20 && j >= 512-K - 30){
						if(r_dist(rng) < (float)th/3){
							F(idxs[j],idxs[j]) = this->F[ind_idx][idxs[j]];
						}else{
							F(idxs[j],idxs[j]) = 0;
						}
					}else{
						F(idxs[j],idxs[j]) = 0;
					}
				}
				
			}
			else
				F(idxs[j],idxs[j]) = this->F[ind_idx][idxs[j]];
		}
		
		//exit(-1);

		mutated_ind2 = population[ind_idx] + F * ((ind1 - population[ind_idx]));
		//mutated_ind2 = population[ind_idx] + F * ((population[idx1] - population[ind_idx]) + (population[idx2] - population[idx3]));
		//aux = idx1;
		//std::cout << mutated_ind2.cols() << " " << mutated_ind2.rows() << std::endl;
	}
#endif
#else
#if RAND_TO_BEST_MOD
	void DE::randtobest_modified(int ind_idx){

		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_int_distribution<int> dist(0, population.size() - 1);
		std::uniform_real_distribution<float> r_dist(0,1);

		int idxb = -1;
		int idx1 = -1;
		int idx2 = -1;
		int idx3 = -1;
		
		idxb = get_best_idx();
		//std::cout << idxb;exit(-1);
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
		Eigen::MatrixXd F(512,512);
		F.setZero(512,512);
		float th = 1./3;

		std::vector<int> idxs = sort_idxs(this->F[ind_idx]);
		for(int j = 0; j < ind_shape[0]; j++){
			if(j < 512-K){
				if(j < 512-K - 10 && j >= 512-K - 20){
					if(r_dist(rng) < th){
						F(idxs[j],idxs[j]) = this->F[ind_idx][idxs[j]];
					}
					else{
						F(idxs[j],idxs[j]) = 0;
					}
				}else{
					if(j < 512-K - 20 && j >= 512-K - 30){
						if(r_dist(rng) < (float)th/3){
							F(idxs[j],idxs[j]) = this->F[ind_idx][idxs[j]];
						}else{
							F(idxs[j],idxs[j]) = 0;
						}
					}else{
						F(idxs[j],idxs[j]) = 0;
					}
				}
				
			}
			else
				F(idxs[j],idxs[j]) = this->F[ind_idx][idxs[j]];
		}
		
		//exit(-1);
		
		mutated_ind = population[idx1] + F * (population[idxb] - population[idx1]);
	}

#else

	void DE::apply_opposition(){
		
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
		#if problem_type == MINIMIZATION && problem_type == MAXIMIZATION
			std::cout << "Error. Problem type MINIMIZATION and MAXIMIZATION macros was both defined as the same value." << std::endl;
			exit(EXIT_FAILURE);

		#elif problem_type == MINIMIZATION
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
		
		#if ROUND_ON_MUTATION
			mutated_ind = truncate_individual(ind_shape, mutated_ind);
		#endif

	}


	//void DE::select_and_change(EvalRankFunction eval_and_rank){
	//	std::vector<Eigen::Matrix2Xd> C;
	//	
	//	C = eval_and_rank(mutated_ind);
	//}

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


		#if ROUND_ON_MUTATION
			mutated_ind = truncate_individual(ind_shape, mutated_ind);
		#endif
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


		#if ROUND_ON_MUTATION
			mutated_ind = truncate_individual(ind_shape, mutated_ind);
		#endif
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


		#if ROUND_ON_MUTATION
			mutated_ind = truncate_individual(ind_shape, mutated_ind);
		#endif
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


		#if ROUND_ON_MUTATION
			mutated_ind = truncate_individual(ind_shape, mutated_ind);
		#endif
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


		#if ROUND_ON_MUTATION
			mutated_ind = truncate_individual(ind_shape, mutated_ind);
		#endif
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


		#if ROUND_ON_MUTATION
			mutated_ind = truncate_individual(ind_shape, mutated_ind);
		#endif
	}

#endif
#endif

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
				//if(r_dist(rng) <= cr || j == J)
				//{
					//crossed_ind(i,j) = mutated_ind(i,j);
					if(!infeasible)
						infeasible = is_infeasible(mutated_ind(i,j));
				//}
				//else
				//	mutated_ind(i,j) = population[ind_idx](i,j);
			}
		}
	}
#if SECOND_MUTATED_FIT
	void DE::bincross_modified(int ind_idx){
		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_real_distribution<float> r_dist(0,1);
		std::uniform_int_distribution<int> dist(0, ind_shape[1] - 2);
		//std::cout << "passei aqui" << std::endl;
		crossed_ind.setZero(512,4);

		for(int i = 0; i < population[ind_idx].rows(); i++){
			float J = dist(rng);
			if(r_dist(rng) <= cr || jr == J)
			{
				for(int j = 0; j < population[ind_idx].cols(); j++)
				{
					crossed_ind(i,j) = mutated_ind(i,j);
					if(!infeasible)
						infeasible = is_infeasible(crossed_ind(i,j));
				}
			}
			else{
				for(int j = 0; j < population[ind_idx].cols(); j++)
				{
					crossed_ind(i,j) = population[ind_idx](i,j);
				}
			}
		}
		//std::cout << crossed_ind.cols() << " " << crossed_ind.rows() << std::endl;
		
	}
#else
	void DE::bincross_modified(int ind_idx){
		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_real_distribution<float> r_dist(0,1);
		std::uniform_int_distribution<int> dist(0, ind_shape[1] - 2);
		crossed_ind.setZero(512,4);
		//std::cout << "passei aqui" << std::endl;

		for(int i = 0; i < population[ind_idx].rows(); i++){
			float J = dist(rng);
			if(r_dist(rng) <= cr || jr == J)
			{
				for(int j = 0; j < population[ind_idx].cols(); j++)
				{
					crossed_ind(i,j) = mutated_ind(i,j);
				}
			}
			else{
				for(int j = 0; j < population[ind_idx].cols(); j++)
				{
					crossed_ind(i,j) = population[ind_idx](i,j);
				}
			}
		}
		//std::cout << crossed_ind.cols() << " " << crossed_ind.rows() << std::endl;
		
	}
#endif
#if SECOND_MUTATED_FIT
	void DE::bincross_modified2(int ind_idx){
		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_real_distribution<float> r_dist(0,1);
		std::uniform_int_distribution<int> dist(0, ind_shape[1] - 2);
		//std::cout << "passei aqui" << std::endl;
		crossed_ind2.setZero(512,4);

		for(int i = 0; i < population[ind_idx].rows(); i++){
			float J = dist(rng);
			if(r_dist(rng) <= cr || jr == J)
			{
				for(int j = 0; j < population[ind_idx].cols(); j++)
				{
					crossed_ind2(i,j) = mutated_ind2(i,j);
					if(!infeasible)
						infeasible = is_infeasible(crossed_ind(i,j));
				}
			}
			else{
				for(int j = 0; j < population[ind_idx].cols(); j++)
				{
					crossed_ind2(i,j) = population[ind_idx](i,j);
				}
			}
		}
		//std::cout << crossed_ind.cols() << " " << crossed_ind.rows() << std::endl;
		
	}

#endif

	void DE::aritcross(int ind_idx){
		
		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_real_distribution<float> r_dist(0,1);
		std::uniform_int_distribution<int> dist(0, ind_shape[1] - 2);
		//std::cout << "passei aqui" << std::endl;

		for(int i = 0; i < population[ind_idx].rows(); i++){
			
			float J = dist(rng);
			for(int j = 0; j < population[ind_idx].cols(); j++)
			{
				if(r_dist(rng) <= cr || j == J)
				{
					mutated_ind(i,j) = 0.5*mutated_ind(i,j) + 0.5 * population[ind_idx](i,j);
					if(!infeasible)
						infeasible = is_infeasible(mutated_ind(i,j));
				}
				else{
					mutated_ind(i,j) = population[ind_idx](i,j);
				}
			}
		}
		
	}

	void DE::aritcross_modified(int ind_idx){
		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_real_distribution<float> r_dist(0,1);
		std::uniform_int_distribution<int> dist(0, ind_shape[1] - 2);
		//std::cout << "passei aqui" << std::endl;
		crossed_ind.setZero(512,4);

		
		for(auto & i : worst_idxs){
			float J = dist(rng);
			for(int j = 0; j < population[ind_idx].cols(); j++)
			{
				if(r_dist(rng) <= cr || j == J)
				{
					crossed_ind(i,j) = (0.2+cross_changer)*mutated_ind(i,j) + (0.8 - cross_changer) * population[N_pop](i,j);
					if(!infeasible)
						infeasible = is_infeasible(mutated_ind(i,j));
				}
				else{
					crossed_ind(i,j) = population[N_pop](i,j);
				}
			}
		}

		for(auto & i : not_worst_idxs){
			for(int j = 0; j < population[ind_idx].cols(); j++)
			{
				crossed_ind(i,j) = population[N_pop](i,j);
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

#if CURRENT_TO_RAND
		currenttorand_modified(ind_idx);
#else
#if RAND_TO_BEST_MOD
		randtobest_modified(ind_idx);
#else
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
#endif
#endif
	}

	void DE::crossover(int ind_idx){
		//std::cout << "alsdjoasikl "<<crossover_algorithm << std::endl;
		switch(crossover_algorithm){
			case 0:
				bincross(ind_idx); break;
			case 1:
				expcross(ind_idx); break;
			case 2:
				aritcross(ind_idx); break;
			case 3:
				bincross_modified(ind_idx); break;
			case 4:
				aritcross_modified(ind_idx); break;
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

	void DE::dir_repair(int ind_idx){
		// int n=0;
		repair_counter++;
		float greater = 24;
		for(int i = 0; i < mutated_ind.rows(); i++){
			for(int j = 0; j < mutated_ind.cols(); j++){
				if(abs(mutated_ind(i,j)) > greater){
					greater = abs(mutated_ind(i,j));
				}
			}
			for(int j = 0; j < mutated_ind.cols(); j++){
				//std::cout << U << " " << L << std::endl;
				if(mutated_ind(i,j) > U || mutated_ind(i,j) < L){
					for(int j = 0; j < mutated_ind.cols(); j++)
						mutated_ind(i,j) /= (greater/24);
					break;
				}
			}
			greater = 24;
		}
	}

	void DE::repair(int ind_idx){

		switch (3)
		{
			case 0:
				uniform_repair_mutated(ind_idx);
#if SECOND_MUTATED_FIT
				uniform_repair_mutated2(ind_idx);
				uniform_repair_crossed(ind_idx);
#endif
				break;
			case 1:
				weibull_repair(ind_idx);
				break;
			case 2:
				F_repair(ind_idx);
				break;	
			case 3:
				dir_repair(ind_idx);
				break;	
			default:
				break;
		}
		
		infeasible = false;
	}

	//void DE::plot_convergence(){
	//	plt::plot(best_fitness);
	//	plt::show;
	//}

	void DE::check_duplicates(){
				
		std::cout << "Checking duplicates..." << std::endl;

		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_int_distribution<int> dist(-24,24);
		std::uniform_real_distribution<float> distr(0,1);

		for(int i = 0; i < N_pop - 1; i++)
		{
			for(int j = i + 1; j <= N_pop - 1; j++)
			{
				if(population[i] == population[j])
				{	
					// float p = distr(rng);
					// for(int k = 0; k < population[j].rows(); k++)
					// {
					// 	for(int l = 0; l < population[j].cols(); l++)
					// 	{
					// 		if(distr(rng) >= p)
					// 			population[j](k,l) = dist(rng);
					// 	}
					// }
					std::cout << "Duplicated idendified... generating new." << std::endl;
					population[j] = generate_individual(ind_shape);

					/* 
					 * in case of change an individual
					 * the check duplicate process
					 * is restarted 
					 */
					i = 0; 
					break;
				}
			
			}
		}
	}

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

	void DE::uniform_repair_mutated(int ind_idx){
		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_real_distribution<float> dist(0,24);
		// int n=0;
		repair_counter++;
		for(int i = 0; i < mutated_ind.rows(); i++){
			for(int j = 0; j < mutated_ind.cols(); j++){
				//std::cout << U << " " << L << std::endl;
				if(mutated_ind(i,j) > 0){
					while(mutated_ind(i,j) > U){
						mutated_ind(i,j) = rand()%24;
					}
				}else{
					while(mutated_ind(i,j) < L){
						mutated_ind(i,j) = -rand()%24;
					}
				}
			}
		}
		
	}
#if SECOND_MUTATED_FIT
	void DE::uniform_repair_mutated2(int ind_idx){
		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_real_distribution<float> dist(0,24);
		// int n=0;
		// /std::cout << mutated_ind2 << std::endl;
		for(int i = 0; i < mutated_ind2.rows(); i++){
			for(int j = 0; j < mutated_ind2.cols(); j++){
				//std::cout << U << " " << L << std::endl;
				if(mutated_ind2(i,j) > 0){
					while(mutated_ind2(i,j) > U){
						mutated_ind2(i,j) = dist(rng);
					}
				}else{
					while(mutated_ind2(i,j) < L){
						mutated_ind2(i,j) = -dist(rng);
					}
				}
			}
		}
	}

	void DE::uniform_repair_crossed(int ind_idx){
		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_real_distribution<float> dist(0,24);
		// int n=0;
		for(int i = 0; i < crossed_ind.rows(); i++){
			for(int j = 0; j < crossed_ind.cols(); j++){
				//std::cout << U << " " << L << std::endl;
				if(crossed_ind(i,j) > 0){
					while(crossed_ind(i,j) > U){
						crossed_ind(i,j) = dist(rng);
					}
				}else{
					while(crossed_ind(i,j) < L){
						crossed_ind(i,j) = -dist(rng);
					}
				}
			}
		}
	}
#endif
	void DE::F_repair(int ind_idx){
		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_real_distribution<float> dist(0,24);

		// int n=0;
		for(int i = 0; i < mutated_ind.rows(); i++){
			for(int j = 0; j < mutated_ind.cols(); j++){
				//std::cout << U << " " << L << std::endl;
				if(mutated_ind(i,j) > 0){
					while(mutated_ind(i,j) > U){
						mutated_ind(i,j) = population[ind_idx](i,j);
					}
				}else{
					while(mutated_ind(i,j) < L){
						mutated_ind(i,j) = population[ind_idx](i,j);
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
#if CURRENT_TO_RAND
#if SECOND_MUTATED_FIT
		std::vector<double> F = eval(truncate_individual(ind_shape, mutated_ind));
#else
		//std::vector<double> F = eval(truncate_individual(ind_shape, crossed_ind));
		

		std::vector<double> vaux = eval(truncate_individual(ind_shape, mutated_ind)), vaux2;
		for(int i = 0; i < vaux.size() - 1; i++)
			vaux2.push_back(vaux[i]);
		std::vector<double> F_mut = vaux2;
		
#endif
#if MEAN_WORST
		std::vector<double> F1, F2;
		F1 = F;
		std::sort(F1.begin(), F1.end(), std::greater<int>());
#endif
		int crossed_fit = 0, mutated_fit = 0;
		int min = 0, min_mut = 0;
		for(int j = 0; j < ind_shape[0]; j++){
			//if(min > F[j])
			//	min = F[j];
			if(min_mut > F_mut[j])
				min_mut = F_mut[j];
		}
		
#if MEAN_WORST
		for(int j = ind_shape[0] - K; j < ind_shape[0]; j++){
			mutated_fit+= F1[j];
		}
		for(int j = 0; j < ind_shape[0]; j++){
			F[j] = F[j] / min;
		}
		mutated_fit /= K;
#else

		for(int j = 0; j < ind_shape[0]; j++){
			//crossed_fit+= F[j];
			mutated_fit+= F_mut[j];
		}
		for(int j = 0; j < ind_shape[0]; j++){
			//F[j] = F[j] / min;
			F_mut[j] = F_mut[j] / min_mut;
		}

#endif
		crossed_fit /= 512;
		mutated_fit /= 512;
#if SECOND_MUTATED_FIT
		int second_mutated_fit = 0;
		int cross_fit2 = 0;
		int cross_fit = 0;
		std::vector<int> all_fit;
		all_fit.push_back(fitness[ind_idx]);
		all_fit.push_back(mutated_fit);
		all_fit.push_back(second_mutated_fit);
		all_fit.push_back(cross_fit);
		all_fit.push_back(cross_fit2);
		this->all_fit.push_back(all_fit);
		if(this->all_fit.size() == N_pop){
			append_fit(this->all_fit);
			this->all_fit.clear();
		}
#else
		std::vector<int> all_fit;
		all_fit.push_back(fitness[ind_idx]);
		all_fit.push_back(mutated_fit);
		//all_fit.push_back(crossed_fit);
		this->all_fit.push_back(all_fit);
#if PLUS_BEST
#else
		if(this->all_fit.size() == N_pop){
			append_fit(this->all_fit);
			this->all_fit.clear();
		}
#endif
#endif

#if SECOND_MUTATED_FIT
		currenttorand_modified2(ind_idx);
		
		uniform_repair_mutated2(ind_idx);
		
		std::vector<double> Fcross = eval(truncate_individual(ind_shape, crossed_ind));
		
		std::vector<double> F2 = eval(truncate_individual(ind_shape, mutated_ind2));
		
		bincross_modified2(ind_idx);
		std::vector<double> Fcross2 = eval(truncate_individual(ind_shape, crossed_ind2));

		
		
		min = 0;
		for(int j = 0; j < ind_shape[0]; j++){
			if(min > F2[j])
				min = F2[j];
		}
		for(int j = 0; j < ind_shape[0]; j++){
			second_mutated_fit+= F2[j];
			F2[j] = F2[j] / min;
		}
		second_mutated_fit /= 512;

		
		min = 0;
		for(int j = 0; j < ind_shape[0]; j++){
			if(min > Fcross[j])
				min = Fcross[j];
		}
		for(int j = 0; j < ind_shape[0]; j++){
			cross_fit+= Fcross[j];
			Fcross[j] = Fcross[j] / min;
		}
		cross_fit /= 512;
//////////////////
		
		min = 0;
		for(int j = 0; j < ind_shape[0]; j++){
			if(min > Fcross2[j])
				min = Fcross2[j];
		}
		for(int j = 0; j < ind_shape[0]; j++){
			cross_fit2+= Fcross2[j];
			Fcross2[j] = Fcross2[j] / min;
		}
		cross_fit2 /= 512;
/////////////////////

		//std::cout << mutated_ind << std::endl;
		
		

		if(problem_type == MINIMIZATION){
			if(mutated_fit < fitness[ind_idx]){
				this->F[ind_idx] = F;
				change_counter++;
				population[ind_idx] = mutated_ind;
				fitness[ind_idx] = mutated_fit;
			}
		}else{
			//mutated_fit << std::endl;
			if(mutated_fit > cross_fit && mutated_fit > second_mutated_fit && mutated_fit > cross_fit2 && mutated_fit > fitness[ind_idx]){
				this->F[ind_idx] = F;
				count_mut1++;
				population[ind_idx] = mutated_ind;
				fitness[ind_idx] = mutated_fit;
			}

			if(cross_fit > mutated_fit && cross_fit > second_mutated_fit && cross_fit > cross_fit2 && cross_fit > fitness[ind_idx]){
				this->F[ind_idx] = Fcross;
				count_cross1++;
				population[ind_idx] = crossed_ind;
				fitness[ind_idx] = cross_fit;
			}

			if(second_mutated_fit > cross_fit && second_mutated_fit > mutated_fit && second_mutated_fit > cross_fit2 && second_mutated_fit > fitness[ind_idx]){
				this->F[ind_idx] = F2;
				count_mut2++;
				population[ind_idx] = mutated_ind2;
				fitness[ind_idx] = second_mutated_fit;
			}

			if(cross_fit2 > cross_fit && cross_fit2 > mutated_fit && cross_fit2 > second_mutated_fit && cross_fit2 > fitness[ind_idx]){
				this->F[ind_idx] = Fcross2;
				count_cross2++;
				population[ind_idx] = crossed_ind2;
				fitness[ind_idx] = cross_fit2;
			}
		}


#else
		//std::cout << mutated_ind << std::endl;
		if(fitness[get_best_idx()] < mutated_fit)
			actual_error = vaux[vaux.size()-1];
		if(!selection_type){
			if(1){
				this->F[ind_idx] = F_mut;
				mut_counter++;
				population[ind_idx] = mutated_ind;
				this->fitness[ind_idx] = mutated_fit;
			}else{
				//this->F[ind_idx] = F;
				crossed_counter++;
				population[ind_idx] = crossed_ind;
				this->fitness[ind_idx] = crossed_fit;
			}
		}else{
			if(1){
				this->F[ind_idx] = F_mut;
				mut_counter++;
				population[ind_idx] = mutated_ind;
				this->fitness[ind_idx] = mutated_fit;
			}else{
				//this->F[ind_idx] = F;
				crossed_counter++;
				population[ind_idx] = crossed_ind;
				this->fitness[ind_idx] = crossed_fit;
			}
		}
		
#endif

		

		

#else
#if RAND_TO_BEST_MOD
		Eigen::MatrixXd F = eval(truncate_individual(ind_shape, mutated_ind));


		int mutated_fit = 0;

		for(int j = 0; j < ind_shape[0]; j++){
			mutated_fit += F(j,j);
		}
		mutated_fit = mutated_fit/ind_shape[0];

		if(problem_type == MINIMIZATION){
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
#else
		#if ROUND_ON_MUTATION
			float mutated_fit = eval(mutated_ind);	
		#else
			float mutated_fit = eval(truncate_individual(ind_shape, mutated_ind));	
		#endif


		//std::cout << ind_idx << std::endl;	
		if(problem_type == MINIMIZATION){
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
#endif
#endif
	}
	int DE::get_counter(int counter){
		if(!counter)
			return mut_counter;
		else
			return crossed_counter;
	}
	void DE::bincross_best(){
		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_real_distribution<float> r_dist(0,1);
		std::uniform_int_distribution<int> dist(0, ind_shape[1] - 2);
		crossed_best.setZero(512,4);
		//std::cout << "passei aqui" << std::endl;

		for(auto & i : worst_idxs){
			float J = dist(rng);
			if(r_dist(rng) <= cr || jr == J)
			{
				for(int j = 0; j < best_ind.cols(); j++)
				{
					crossed_best(i,j) = (0.8 - cross_changer)*best_ind(i,j) + (0.2+cross_changer) * population[get_best_idx()](i,j);;
					if(!infeasible)
						infeasible = is_infeasible(crossed_best(i,j));
				}
			}
			else{
				for(int j = 0; j < best_ind.cols(); j++)
				{
					crossed_best(i,j) = best_ind(i,j);
				}
			}
		}
		for(auto & i : not_worst_idxs){
			for(int j = 0; j < best_ind.cols(); j++)
			{
				crossed_best(i,j) = best_ind(i,j);
			}
		}
		//std::cout << crossed_ind.cols() << " " << crossed_ind.rows() << std::endl;
		
	}
	void DE::currenttorand_modified(Eigen::MatrixXd ind){
		
		std::random_device rseed;
		std::mt19937 rng(rseed());
		std::uniform_int_distribution<int> dist(0, population.size() - 1);
		std::uniform_real_distribution<float> r_dist(0,1);
		
		Eigen::MatrixXd ind1 = generate_individual(ind_shape), ind2 = generate_individual(ind_shape),ind3 = generate_individual(ind_shape), ind4 = generate_individual(ind_shape);
		
		int idx1 = -1;
		int idx2 = -1;
		int idx3 = -1;
		

		Eigen::MatrixXd F(512,512);
		F.setZero(512,512);
		float th = 1./3;

		std::vector<int> idxs = sort_idxs(this->F[get_best_idx()]);
		Eigen::MatrixXd aux = ind;
		for(int i = 0; i < K; i++){
			worst_idxs[i] = idxs[511 - i];
		}
		for(int i = 0; i < 512-K; i++){
			not_worst_idxs[i] = idxs[i];
		}
		for(int j = 0; j < ind_shape[0]; j++){
			if(j < 512-K){
			}
			else{
				F(idxs[j],idxs[j]) = this->F[get_best_idx()][idxs[j]];
				for(int i = 0; i < 4; i++){
					aux(idxs[j], i) = 0;
				}
			}
		}
		
		best_ind = aux + F/N_pop * (ind1 - ind2);
		
		//mutated_ind = population[ind_idx] + F * ((population[idx1] - population[ind_idx]) + (population[idx2] - population[idx3]));
	}

	Eigen::MatrixXd DE::extra_dir_repair(Eigen::MatrixXd mutated_ind){
		// int n=0;
		repair_counter++;
		float greater = 24;
		for(int i = 0; i < mutated_ind.rows(); i++){
			for(int j = 0; j < mutated_ind.cols(); j++){
				if(abs(mutated_ind(i,j)) > greater){
					greater = abs(mutated_ind(i,j));
				}
			}
			for(int j = 0; j < mutated_ind.cols(); j++){
				//std::cout << U << " " << L << std::endl;
				if(mutated_ind(i,j) > U || mutated_ind(i,j) < L){
					for(int j = 0; j < mutated_ind.cols(); j++)
						mutated_ind(i,j) /= (greater/24);
					break;
				}
			}
			
			greater = 24;
		}
		return mutated_ind;
	}

	void DE::extra_selection(int ind_idx){
		
		//if(fitness[N_pop] < fitness[get_best_idx()]){
		//	population[N_pop] = population[get_best_idx()];
		//	fitness[N_pop] = fitness[get_best_idx()];
		//	best_ind = population[get_best_idx()];
		//}
		//std::cout << get_best_idx() << std::endl;
		//std::cout << "best ind is" << population[N_pop] << std::endl;
		currenttorand_modified(truncate_individual(ind_shape, population[N_pop]));
		best_ind = population[N_pop];
		//best_ind = extra_dir_repair(best_ind);
		bincross_best();
		//crossed_best = extra_dir_repair(crossed_best);
		
		std::vector<double> vaux = eval(truncate_individual(ind_shape, crossed_best)), vaux2;
		for(int i = 0; i < vaux.size() - 1; i++)
			vaux2.push_back(vaux[i]);
		std::vector<double> F = vaux2;
		


		int crossed_fit = 0;
		int min = 0;
		for(int j = 0; j < ind_shape[0]; j++){
			if(min > F[j])
				min = F[j];
		}

		for(int j = 0; j < ind_shape[0]; j++){
			crossed_fit+= F[j];
			F[j] = F[j] / min;
		}

		crossed_fit /= 512;
		
		

		std::vector<int> all_fit;
		all_fit.push_back(fitness[N_pop]);
		all_fit.push_back(crossed_fit);

		
		this->all_fit.push_back(all_fit);

		append_fit(this->all_fit);

		this->all_fit.clear();
		
		if(problem_type == MINIMIZATION){
			if(crossed_fit < fitness[N_pop]){
				this->F[N_pop] = F;
				change_counter++;
				population[N_pop] = best_ind;
				fitness[N_pop] = crossed_fit;
				actual_error = vaux[vaux.size()-1];
				//best_ind = crossed_best;
			}else{
				this->F[N_pop] = this->F[N_pop];
				change_counter++;
				population[N_pop] = population[N_pop];
				fitness[N_pop] = fitness[N_pop];
				best_ind = population[N_pop];
			}
		}else{
			//mutated_fit << std::endl;
			if(crossed_fit > fitness[N_pop]){
				this->F[N_pop] = F;
				change_counter++;
				population[N_pop] = best_ind;
				fitness[N_pop] = crossed_fit;
				best_ind = crossed_best;
				actual_error = vaux[vaux.size()-1];
			}else{
				this->F[N_pop] = this->F[N_pop];
				change_counter++;
				fitness[N_pop] = fitness[get_best_idx()];
				population[N_pop] = population[get_best_idx()];
				best_ind = population[N_pop];
			}
		}

	}

	void DE::set_change_counter(uint value){
		change_counter = value;
	}

	void DE::change_cross(int g, int ng){
		cross_changer = (0.8-0.2)*(g+1)/ng;
		if(cross_changer < 0)
			cross_changer = 0 ;
		std::cout << "g is " << g << " and ng " << ng << std::endl;
		std::cout << "cross is " << cross_changer << std::endl;
	}
	//void DE::set_best_fit(){
	//	best_fitness.emplace_back(get_best_fit());
	//}

	void DE::evolve(uint ng){
		
		for(int g = 0; g < ng; g++){
			change_counter  = 0;
			check_duplicates();

			for(int i = 0; i < population.size(); i++){
				mutate(i);
				if(infeasible)
					repair(i);
				crossover(i);
				selection(i);					
			}
			
			//best_fitness.emplace_back(get_best_fit());
		}
	}

	bool DE::is_infeasible(){
		return infeasible;
	}

#if CURRENT_TO_RAND
	int DE::get_best_idx(){
		if(problem_type ==MINIMIZATION)
			return std::min_element(this->fitness.begin(), this->fitness.end()) - this->fitness.begin();
		else
			return std::max_element(this->fitness.begin(), this->fitness.end()) - this->fitness.begin();
	}
	Eigen::MatrixXd DE::get_best_ind(){
		if(problem_type == MINIMIZATION)
			return population[std::min_element(this->fitness.begin(), this->fitness.end()) - fitness.begin()];
		else
			return population[std::max_element(this->fitness.begin(), this->fitness.end()) - fitness.begin()];
	}

	float DE::get_best_fit(){
		if(problem_type == MINIMIZATION)
			return *std::min_element(this->fitness.begin(), this->fitness.end());
		else
			return *std::max_element(this->fitness.begin(), this->fitness.end());
	}

#else

#if RAND_TO_BEST_MOD
	float DE::get_best_fit(){
		if(problem_type == MINIMIZATION)
			return *std::min_element(this->fitness.begin(), this->fitness.end());
		else
			return *std::max_element(this->fitness.begin(), this->fitness.end());
	}

	Eigen::MatrixXd DE::get_best_ind(){
		if(problem_type == MINIMIZATION)
			return population[std::min_element(this->fitness.begin(), this->fitness.end()) - fitness.begin()];
		else
			return population[std::max_element(this->fitness.begin(), this->fitness.end()) - fitness.begin()];
	}
#else
	float DE::get_best_fit(){
		if(problem_type == MINIMIZATION)
			return *std::min_element(this->fitness.begin(), this->fitness.end());
		else
			return *std::max_element(this->fitness.begin(), this->fitness.end());
	}

	Eigen::MatrixXd DE::get_best_ind(){
		if(problem_type == MINIMIZATION)
			return population[std::min_element(this->fitness.begin(), this->fitness.end()) - fitness.begin()];
		else
			return population[std::max_element(this->fitness.begin(), this->fitness.end()) - fitness.begin()];
	}
#endif
#if RAND_TO_BEST_MOD
	int DE::get_best_idx(){
		if(problem_type ==MINIMIZATION)
			return std::min_element(this->fitness_aux.begin(), this->fitness_aux.end()) - fitness_aux.begin();
		else
			return std::max_element(this->fitness_aux.begin(), this->fitness_aux.end()) - fitness_aux.begin();
	}
#else
	int DE::get_best_idx(){
		if(problem_type == MINIMIZATION)
			return std::min_element(this->fitness.begin(), this->fitness.end()) - fitness.begin();
		else
			return std::max_element(this->fitness.begin(), this->fitness.end()) - fitness.begin();
	}
#endif

	void DE::get_fitness(){
		for(int i = 0; i < 30; i++)
			std::cout << fitness[i] << std::endl;
	}

	int DE::get_max_elem(){
		return population[0].maxCoeff();
	}

#endif

}