// weibull_distribution
#include <iostream>
#include <random>


int main()
{
	const int nrolls=100000;  // number of experiments
	const int nstars=1000;    // maximum number of stars to distribute

	std::random_device rseed;
	std::mt19937 rng(rseed());
	std::weibull_distribution<double> distribution(2.0,23.0);

	int p[24]={};

	for (int i=0; i<nrolls; ++i) {
		double number = distribution(rng);
		if (number<24) ++p[int(number)];
	}

	std::cout << "weibull_distribution (2.0,4.0):" << std::endl;

	for (int i=0; i<24; ++i) {
		std::cout << i << "-" << (i+1) << ": ";
		std::cout << std::string(p[i]*nstars/nrolls,'*') << std::endl;
	}

	return 0;
}