#include <stdlib.h>
#include <iostream>
#include <vector>
#include "quicksort.h"

int main(){

	int i = 0, n = 10, op = 0;
	
	std::vector<int> population_indexes;
	population_indexes.reserve(n);

	for(i = 0; i < n; i++)
		population_indexes[i] = i;

	std::vector<float> fitness_vector = {100, 8, 90, -1, 3, 0, 5, 1, -5, -10};

	for (i = 0; i < n; i++)
		printf("%f  %d \n", fitness_vector[i], population_indexes[i]);
	
	printf("\n\n");

	QS::quicksort qs;
	qs.sort(fitness_vector, population_indexes, 0, n-1);

	for (int i = 0; i < n; i++)
		printf("%f  %d \n", fitness_vector[i], population_indexes[i]);
	
	return 0;
}
