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

void Measurements::append_fit(std::vector<std::vector<int>> all_fit){
	all_fits_ind.push_back(all_fit);
}

void Measurements::reset(){
	boxes.clear();
	means.clear();
	std_devs.clear();
	all_fits_ind.clear();
	gen = 1;
}

void Measurements::save_data(std::vector<float> y, const std::string &dataset, const std::string &exp, Eigen::MatrixXd best_individual, int count1, int count2, int count3, int count4, std::vector<std::vector<double>> F, int repair_counter){
	
	
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
	std::ofstream f5(CURRENT_DIR +"/../results/" + dataset+ "/" + exp + "/" + "all_fits.txt");
	std::ofstream f6;
	f6.open(CURRENT_DIR +"/../results/" + dataset+ "/" + exp + "/" + "all_F.txt", std::ios_base::app);

	for(std::vector<float>::const_iterator i = y.begin(); i != y.end(); ++i) {
    	f1 << *i << '\n';
	}

	

	if (f2.is_open())
  	{
    	f2 << best_individual;
  	}

	if (f3.is_open())
  	{
		  f3 << "mut 1" << " " << "mut 2" << " " << "cross 1" << " " << "cross 2" << " repair_counter" << "\n";
    	f3 << count1 << " " << count2 << " " << count3 << " " << count4 << " " << repair_counter;
  	}


	for(std::vector<double>::const_iterator i = std_devs.begin(); i != std_devs.end(); ++i) {
    	f4 << *i << '\n';
	}

	

	gen = 1;
	for(int i = 0; i < all_fits_ind.size(); i++){
		f5 << "Gen " << gen << std::endl;
		for(int j = 0; j < all_fits_ind[i].size();j++){
			f5 << "ind " << j + 1 << ": ";
			for(int k = 0; k < all_fits_ind[i][j].size(); k++){
				f5 << all_fits_ind[i][j][k] << " ";
			}
			f5 << '\n';
		}
		gen++;
	}
	gen--;

	f6 << "Gen " << gen << std::endl;
	for(int i = 0; i < F.size(); i++){
		f6 << "ind " << i << ": ";
		for(int j = 0; j < F[i].size(); j++){
			f6 << F[i][j] << " ";
		}
		f6 << "\n";

	}
	gen++;

	
	//plt::show();
}
