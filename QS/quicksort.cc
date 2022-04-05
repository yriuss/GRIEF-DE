#include "quicksort.h"

namespace QS{

    int quicksort::sort(float fitness_vector[], int population_indexes[], int left, int right){
        int i = 0, j = 0;

        partition(fitness_vector, population_indexes, left, right, &i, &j);

        if (left < j)
            sort(fitness_vector, population_indexes, left, j);
        if (i < right)
            sort(fitness_vector, population_indexes, i, right);

        return (fitness_vector[0], population_indexes[0]);
    }

    void quicksort::partition(float fitness_vector[], int population_indexes[], int left, int right, int *i, int *j){
        
        int x = 0, aux = 0;
        *i = left;
        *j = right;
        x = fitness_vector[(*i + *j)/2]; 

        do{
            while (x > fitness_vector[*i])
                (*i)++;

            while (x < fitness_vector[*j])
                (*j)--;

            if(*i <= *j){
                aux = fitness_vector[*i];
                fitness_vector[*i] = fitness_vector[*j];
                fitness_vector[*j] = aux;

                aux = population_indexes[*i];
                population_indexes[*i] = population_indexes[*j];
                population_indexes[*j] = aux;

                (*i)++;
                (*j)--;
            }
        }while(*i <= *j);
    }
}