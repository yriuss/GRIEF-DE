#ifndef DE_H_INCLUDED
#define DE_H_INCLUDED

#include <random>
#include <stdio.h>
#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <unistd.h>
#include <fstream>
#include <cmath>
#include "QS/quicksort.h"
#include "../measurements/measurements.h"

#define CURRENT_TO_RAND true
#define READ_BEST_IND false
#define RAND_TO_BEST_MOD false
#define MEAN_WORST false
#define SECOND_MUTATED_FIT false
 
#define BIN_CROSS_GENE true

#if ( CURRENT_TO_RAND || RAND_TO_BEST_MOD ) && BIN_CROSS_GENE
	typedef void (*EvalFunction)(Eigen::MatrixXd, std::vector<double> &, std::vector<float> &);
#elif CURRENT_TO_RAND || RAND_TO_BEST_MOD
	typedef std::vector<double>(*EvalFunction)(Eigen::MatrixXd);
#elif BIN_CROSS_GENE
	typedef void (*EvalFunction)(Eigen::MatrixXd, float &, std::vector<float> &);
#else
	typedef float(*EvalFunction)(Eigen::MatrixXd);
#endif


namespace DE {


/* DEFINITION OF PROBLEM TYPE */
#define MAXIMIZATION true
#define MINIMIZATION false
#define DEBUG 1

/* DEFINITION OF MUTATION ALGORITHM */
#define            RAND_1 0 
#define            RAND_2 1
#define      RAND_TO_BEST 2
#define            BEST_1 3
#define            BEST_2 4
#define CURRENT_TO_BEST_1 5
#define CURRENT_TO_RAND_1 6

/* DEFINITION OF CROSSOVER ALGORITHM */
#define EXP   1
#define ARIT  2
#define BIN   3
#define BIN_G 4

#define OPPOSITION_LEARNING false
#define ROUND_ON_MUTATION   true
	
	class DE: public Measurements{

		public:
			DE(int N_pop, std::vector<int> ind_shape, float cr, float jr,
			EvalFunction evaluation, float F, bool problem_type, std::vector<int> bounds, 
			int mutation_algorithm, int crossover_algorithm, int K);
			void reset();
			
			Eigen::MatrixXd generate_individual(std::vector<int> ind_shape);
			Eigen::MatrixXd generate_oppsite_individual(std::vector<int> ind_shape, int ind_idx);
			Eigen::MatrixXd truncate_individual(std::vector<int> ind_shape, Eigen::MatrixXd ind);
			void F_repair(int ind_idx);
			void generate_oppsite_population();
			void apply_opposition();
			void crossover(int ind_idx);
			void mutate(int ind_idx);
			void repair(int ind_idx);
			void evolve(uint ng);
			void evaluate(int ind_idx);
			void selection(int ind_idx);
			void create_population();
			void get_fitness();
			void weibull_repair(int ind_idx);
			void uniform_repair_mutated(int ind_idx);
			void bincross_modified(int ind_idx);
			void bincross(int ind_idx);
			void aritcross_modified(int ind_idx);
			void expcross(int ind_idx);
			void aritcross(int ind_idx);
			float get_best_fit();
			int get_best_idx();
			void bincross_modified2(int ind_idx);
			int get_max_elem();
			std::vector<Eigen::MatrixXd> pop();
			bool is_infeasible(int element);
			bool is_infeasible();
			// void select_and_change(EvalRankFunction evaluation);
			//void plot_convergence();
			//void set_best_fit();
			float jr;
			uint get_change_counter();
			void set_change_counter(uint value);
			void check_duplicates();
			
			Eigen::MatrixXd get_best_ind();

			void read_individuals(int n_of_individuals);


			#if BIN_CROSS_GENE
				
				float normalize(float x, float min, float max);
				void  bincrossnorm(int ind_idx);

			#endif


			#if CURRENT_TO_RAND

				#if SECOND_MUTATED_FIT

					void uniform_repair_mutated2(int ind_idx);
					void uniform_repair_crossed(int ind_idx);
					void currenttorand_modified2(int ind_idx);

				#endif

					void currenttorand_modified(int ind_idx);

			#else
				#if RAND_TO_BEST_MOD

					void randtobest_modified(int ind_idx);

				#else
					
					void rand_1(int ind_idx);
					void rand_2(int ind_idx);
					void randtobest_1(int ind_idx);
					void best_1(int ind_idx);
					void best_2(int ind_idx);
					void currenttobest_1(int ind_idx);
					void currenttorand_1(int ind_idx);

				#endif
			#endif

			int count_mut1 = 0;
			int count_cross1 = 0;
			int count_mut2 = 0;
			int count_cross2 = 0;

		private:

			EvalFunction eval;
			Eigen::MatrixXd mutated_ind;

			

			#if SECOND_MUTATED_FIT
				Eigen::MatrixXd mutated_ind2;
			#endif

			Eigen::MatrixXd crossed_ind;
			Eigen::MatrixXd crossed_ind2;
			std::vector<Eigen::MatrixXd> population;

			#if OPPOSITION_LEARNING
				std::vector<Eigen::MatrixXd> opposite_population;
				std::vector<float> opposite_fitness;
			#endif

			#if OPPOSITION_LEARNING && BIN_CROSS_GENE
				std::vector<vector<float>> opposite_gene_fitness;
			#endif

			// std::vector<float> ;
			//std::vector<float> best_fitness;
			std::vector<int> ind_shape;
			float cr;
			std::vector<int> sort_idxs(std::vector<double> v);

			#if CURRENT_TO_RAND || RAND_TO_BEST_MOD
			
				//std::vector<Eigen::MatrixXd> F;
				float K = 10;
				std::vector<std::vector<double>> F;
				float F_mut;
				std::vector<float> fitness;
				std::vector<float> fitness_aux;

			#else

				float F;
				std::vector<float> fitness;

			#endif


			#if BIN_CROSS_GENE				

				std::vector<std::vector<float>> gene_fitness;

			#endif

			uint change_counter = 0;
			// float jr;
			bool problem_type;
			bool infeasible = false;
			
			#if SECOND_MUTATED_FIT
				bool infeasible2 = false;
				bool infeasible3 = false;
			#endif

			int U;
			int aux;
			int L;
			int mutation_algorithm;
			int crossover_algorithm;
			int N_pop;
			
			
	};
}

#endif