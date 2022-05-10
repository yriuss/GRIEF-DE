#include <iostream>
#include <stdio.h>
#include <Eigen/Dense>
#include <vector>
#include <unistd.h>
#include <fstream>

Eigen::MatrixXd mutated_ind;

std::vector<Eigen::MatrixXd> population;


void read_individuals(int n_of_individuals){
	std::string CURRENT_DIR = get_current_dir_name();
	using namespace std;
	for(int i = 0; i < n_of_individuals; i++){
		//std::cout << CURRENT_DIR +"/../best_inds/" + std::to_string(i) + ".txt";exit(-1);
		ifstream file(CURRENT_DIR +"/best_inds/" + std::to_string(i*7) + ".txt");
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

void rand_1(){
    float F = 0.8;
    
	mutated_ind = population[0] + F * (population[1] - population[0]);
    std::cout << "mutated individual: " << std::endl <<  mutated_ind;
}


int main(){
    population.reserve(3);
    read_individuals(3);
    std::cout << "current individual: " << std::endl << population[0] << std::endl;
    rand_1();
    return 0;
}