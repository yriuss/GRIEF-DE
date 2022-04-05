#ifndef DE_H_INCLUDED
#define DE_H_INCLUDED

#include <random>
#include <stdio.h>
#include <iostream>
#include <Eigen/Dense>
#include <vector>

typedef float(*EvalFunction)(Eigen::MatrixXd);

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
#define BIN 0
#define EXP 1

#define OPPOSITION_LEARNING true
#define LOCAL_SEARCH true

	class DE{

		public:
			DE(int N_pop, std::vector<int> ind_shape, float cr, 
			EvalFunction evaluation, float F, bool problem_type, std::vector<int> bounds, 
			int mutation_algorithm, int crossover_algorithm);

			Eigen::MatrixXd get_best_ind();
			Eigen::MatrixXd generate_individual(std::vector<int> ind_shape);
			Eigen::MatrixXd generate_oppsite_individual(std::vector<int> ind_shape, int ind_idx);

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
			void rand_1(int ind_idx);
			void rand_2(int ind_idx);
			void randtobest_1(int ind_idx);
			void best_1(int ind_idx);
			void best_2(int ind_idx);
			void currenttobest_1(int ind_idx);
			void currenttorand_1(int ind_idx);
			void bincross(int ind_idx);
			void expcross(int ind_idx);
			float get_best_fit();
			int get_best_idx();
			int get_max_elem();
			bool is_infeasible(int element);
			bool is_infeasible();
			//void plot_convergence();
			//void set_best_fit();

		private:
			EvalFunction eval;
			Eigen::MatrixXd mutated_ind;
			std::vector<Eigen::MatrixXd> population;
			std::vector<Eigen::MatrixXd> opposite_population;
			std::vector<float> fitness;
			std::vector<float> opposite_fitness;
			//std::vector<float> best_fitness;
			std::vector<int> ind_shape;
			float cr;
			float F;
			bool problem_type;
			int U;
			int L;
			bool infeasible = false;
			int mutation_algorithm;
			int crossover_algorithm;
			int N_pop;

	};

}

#endif