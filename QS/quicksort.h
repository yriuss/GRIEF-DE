#include <stdio.h>
#include <iostream>
#include <vector>

#ifndef QUICKSORT_H_INCLUDED
#define QUICKSORT_H_INCLUDED

namespace QS{
	class quicksort
	{
		public:
			void sort(std::vector<float> &fitness_vector, std::vector<int> &population_indexes, int left, int right);
		
		private:
			void partition(std::vector<float> &fitness_vector, std::vector<int> &population_indexes, int left, int right, int *i, int *j);	
	};
}

#endif