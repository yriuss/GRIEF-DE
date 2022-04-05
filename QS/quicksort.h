#include <stdio.h>
#include <iostream>

#ifndef QUICKSORT_H_INCLUDED
#define QUICKSORT_H_INCLUDED

namespace QS{
	class quicksort
	{
		public:
			int sort(float fitness_vector[], int population_indexes[], int left, int right);
		
		private:
			void partition(float fitness_vector[], int population_indexes[], int left, int right, int *i, int *j);	
	};
}

#endif