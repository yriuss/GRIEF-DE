#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <array>
#include <vector>
using namespace std;

float normalize(float x, float min, float max){
    return (x - min) / (max - min);
}

int main(){

    // std::array<int,5> v {-25000, -65000, -30000, -10000, -36000};
    std::vector<float> v = {-25000, -65000, -30000, -10000, -36000, 1000};
    std::vector<float> nv;
    nv.reserve(v.size());

    auto minmax = std::minmax_element(v.begin(), v.end());

    for(int i=0; i<v.size(); i++){
        nv[i] = normalize(v[i], *minmax.first, *minmax.second);
    }

    std::cout << "Min: " << *minmax.first << std::endl;
    std::cout << "Max: " << *minmax.second << std::endl;

    for(int i=0; i<v.size(); i++){
        std::cout << " Value: " << v[i] << " Norm: " << nv[i] << std::endl;
    }

    return 0;
}