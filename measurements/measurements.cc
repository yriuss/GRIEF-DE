#include "measurements.h"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

bool dir_exist(const std::string &s)
{
  struct stat buffer;
  return (stat (s.c_str(), &buffer) == 0);
}

void _mkdir(const std::string &s){
	mkdir((s).c_str(), S_IRUSR | S_IWUSR | S_IXUSR);
}

Eigen::MatrixXd Measurements::mean(std::vector<Eigen::MatrixXd> x, bool save){
	Eigen::MatrixXd mean(512,4);
	mean << Eigen::MatrixXd::Zero(512, 4);

	for(int k = 0; k < x.size(); k++)
		for(int i = 0; i < x[k].rows(); i++)
			for(int j = 0; j < x[k].cols(); j++)
				mean(i,j) += x[k](i,j);
	mean /= x.size();

	if(save)
		means.push_back(mean);
	return mean;
}

double Measurements::std_dev(std::vector<Eigen::MatrixXd> x, bool save){
	Eigen::MatrixXd mean = this->mean(x);
	Eigen::MatrixXd std_dev(512, 4);
	std_dev << Eigen::MatrixXd::Zero(512, 4);
	long double std = 0;
	//std::cout << x[0].rows() << x[0].cols() << std::endl;
	for(int k = 0; k < x.size(); k++)
		for(int i = 0; i < x[k].rows(); i++)
			for(int j = 0; j < x[k].cols(); j++)
				std_dev(i,j) += (x[k](i,j)-mean(i,j))*(x[k](i,j)-mean(i,j));
	std_dev /= (x.size()-1);
	//std::cout << std_dev << std::endl;
	for(int i = 0; i < std_dev.rows(); i++){
		for(int j = 0; j < std_dev.cols(); j++){
			std_dev(i,j) = sqrt(std_dev(i,j));
		}
	}
	

	for(int i = 0; i < std_dev .rows(); i++)
		for(int j = 0; j < std_dev .cols(); j++)
			std += std_dev(i,j);
	//std::cout << std << std::endl;
	std /= 512*4;
	
	if(save)
		std_devs.push_back(std);
	
	return std;
}

void Measurements::make_boxes(){
	//for(int i = 0; i < means.size(); i++){
	//	boxes.push_back(std::vector<double>{means[i] - std_devs[i], means[i], means[i] + std_devs[i]});
	//}
	
}

void Measurements::append_fit1(std::vector<int> all_fit){
	all_fits_ind1.push_back(all_fit);
}
void Measurements::append_fit2(std::vector<int> all_fit){
	all_fits_ind2.push_back(all_fit);
}
void Measurements::append_fit3(std::vector<int> all_fit){
	all_fits_ind3.push_back(all_fit);
}
void Measurements::append_fit4(std::vector<int> all_fit){
	all_fits_ind4.push_back(all_fit);
}

void Measurements::reset(){
	boxes.clear();
	means.clear();
	std_devs.clear();
}

void Measurements::save_data(std::vector<float> y, const std::string &dataset, const std::string &exp, Eigen::MatrixXd best_individual, int count1, int count2, int count3, int count4){
	

	if(!dir_exist(CURRENT_DIR+"/../results/"))
		_mkdir(CURRENT_DIR+"/../results/");
	
	if(!dir_exist(CURRENT_DIR+"/../results/"  + dataset))
		_mkdir(CURRENT_DIR+"/../results/" + dataset);
	
	if(!dir_exist(CURRENT_DIR+"/../results/"  + dataset + "/" + exp))
		_mkdir(CURRENT_DIR+"/../results/" + dataset+ "/" + exp);
	
	plt::plot(y);
	plt::title("Convergence " + dataset);
	plt::xlabel("Gerações");
	plt::ylabel("Fitness");
	plt::save(CURRENT_DIR+"/../results/" + dataset + "/" + exp + "/" + "convergence.png");
	plt::cla();

	plt::plot(std_devs);
	plt::title("Convergence " + dataset);
	plt::xlabel("Gerações");
	plt::ylabel("std dev");
	plt::save(CURRENT_DIR+"/../results/" + dataset + "/" + exp + "/" + "behavior.png");
	plt::cla();
	boxes.clear();
	std::ofstream f1(CURRENT_DIR +"/../results/" + dataset+ "/" + exp + "/" + "convergence.txt"), f2(CURRENT_DIR +"/../results/" + dataset+ "/" + exp + "/" + "best_individual.txt");
	std::ofstream f3(CURRENT_DIR +"/../results/" + dataset+ "/" + exp + "/" + "changes.txt"), f4(CURRENT_DIR +"/../results/" + dataset+ "/" + exp + "/" + "std_devs.txt");
	std::ofstream f5(CURRENT_DIR +"/../results/" + dataset+ "/" + exp + "/" + "all_fits1.txt");
	std::ofstream f6(CURRENT_DIR +"/../results/" + dataset+ "/" + exp + "/" + "all_fits2.txt");
	std::ofstream f7(CURRENT_DIR +"/../results/" + dataset+ "/" + exp + "/" + "all_fits3.txt");
	std::ofstream f8(CURRENT_DIR +"/../results/" + dataset+ "/" + exp + "/" + "all_fits4.txt");
	for(std::vector<float>::const_iterator i = y.begin(); i != y.end(); ++i) {
    	f1 << *i << '\n';
	}

	

	if (f2.is_open())
  	{
    	f2 << best_individual;
  	}

	if (f3.is_open())
  	{
		  f3 << "mut 1" << " " << "mut 2" << " " << "cross 1" << " " << "cross 2" << "\n";
    	f3 << count1 << " " << count2 << " " << count3 << " " << count4;
  	}


	for(std::vector<double>::const_iterator i = std_devs.begin(); i != std_devs.end(); ++i) {
    	f4 << *i << '\n';
	}

	for(int i = 0; i < all_fits_ind1.size(); i++){
		for(int j = 0; j < all_fits_ind1[i].size(); j++){
			f5 << all_fits_ind1[i][j] << " ";
		}
		f5 << '\n';
	}

	for(int i = 0; i < all_fits_ind2.size(); i++){
		for(int j = 0; j < all_fits_ind2[i].size(); j++){
			f6 << all_fits_ind2[i][j] << " ";
		}
		f6 << '\n';
	}

	for(int i = 0; i < all_fits_ind3.size(); i++){
		for(int j = 0; j < all_fits_ind3[i].size(); j++){
			f7 << all_fits_ind3[i][j] << " ";
		}
		f7 << '\n';
	}

	for(int i = 0; i < all_fits_ind4.size(); i++){
		for(int j = 0; j < all_fits_ind4[i].size(); j++){
			f8 << all_fits_ind4[i][j] << " ";
		}
		f8 << '\n';
	}
	//plt::show();
}
