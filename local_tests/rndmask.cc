#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>

using namespace std;

int main(){

    std::random_device rseed;
    std::mt19937 rng(rseed());
    std::uniform_real_distribution<float> r_dist(0,1);
    std::uniform_int_distribution<int> dist(0, 1);
    
    for(int i=0; i<10; i++)
        std::cout << dist(rng) << std::endl;

    return 0;
}