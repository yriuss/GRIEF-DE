#ifndef DE_H_INCLUDED
#define DE_H_INCLUDED

#include <random>
#include <stdio.h>
#include <iostream>
#include <Eigen/Dense>
#include <vector>

#define CURRENT_TO_RAND false

#if CURRENT_TO_RAND
typedef Eigen::MatrixXd(*EvalFunction)(Eigen::MatrixXd);
#else
typedef float(*EvalFunction)(Eigen::MatrixXd);
#endif


namespace DE {



/* DEFINITION OF PROBLEM TYPE */
#define MAXIMIZATION false
#define MINIMIZATION true
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
#define BIN 0
#define EXP 1
#define ARIT 2

#define OPPOSITION_LEARNING false
#define ROUND_ON_MUTATION  true
	
	class DE{

		public:
			DE(int N_pop, std::vector<int> ind_shape, float cr, float jr,
			EvalFunction evaluation, float F, bool problem_type, std::vector<int> bounds, 
			int mutation_algorithm, int crossover_algorithm);

			
			Eigen::MatrixXd generate_individual(std::vector<int> ind_shape);
			Eigen::MatrixXd generate_oppsite_individual(std::vector<int> ind_shape, int ind_idx);
			Eigen::MatrixXd truncate_individual(std::vector<int> ind_shape, Eigen::MatrixXd ind);

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
			void uniform_repair(int ind_idx);
			
			void bincross(int ind_idx);
			void expcross(int ind_idx);
			void aritcross(int ind_idx);
			float get_best_fit();
			int get_best_idx();
			int get_max_elem();
			bool is_infeasible(int element);
			bool is_infeasible();
			// void select_and_change(EvalRankFunction evaluation);
			//void plot_convergence();
			//void set_best_fit();
			float jr;
			uint get_change_counter();
			void set_change_counter(uint value);
			void check_duplicates();
			float indv_mean();
			float indv_variance();
			float indv_standard_deviation();
			Eigen::MatrixXd get_best_ind();

#if CURRENT_TO_RAND

			void currenttorand_modified(int ind_idx);
#else
			
			void rand_1(int ind_idx);
			void rand_2(int ind_idx);
			void randtobest_1(int ind_idx);
			void best_1(int ind_idx);
			void best_2(int ind_idx);
			void currenttobest_1(int ind_idx);
			void currenttorand_1(int ind_idx);
#endif

		private:
			EvalFunction eval;
			Eigen::MatrixXd mutated_ind;
			std::vector<Eigen::MatrixXd> population;
			std::vector<Eigen::MatrixXd> opposite_population;
			
			std::vector<float> opposite_fitness;
			// std::vector<float> ;
			//std::vector<float> best_fitness;
			std::vector<int> ind_shape;
			float cr;

#if CURRENT_TO_RAND

			std::vector<Eigen::MatrixXd> F;
			std::vector<float> fitness;
#else

			float F;
			std::vector<float> fitness;

#endif
			uint change_counter = 0;
			// float jr;
			bool problem_type;
			bool infeasible = false;
			int U;
			int L;
			int mutation_algorithm;
			int crossover_algorithm;
			int N_pop;
	};
}

#endif