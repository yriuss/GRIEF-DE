#include "DE.h"
#include <random>
#include "matplotlibcpp.h"
#include <vector>
namespace plt = matplotlibcpp;

void DE::DE(){
    
}




void plot_convergence(std::vector<int> x,std::vector<int> y){
    plt::plot(x,y);
    plt::show();
}

static void evaluation(){
    
}

int main(){
    std::vector<int> x, y;

    std::random_device rseed;
    std::mt19937 rng(rseed());
    std::uniform_int_distribution<int> dist(-24,24);

    std::cout << dist(rng) << std::endl;
    std::cout << dist(rng) << std::endl;

    if(x.size() > 0 && y.size() > 0)
        plot_convergence(x,y);
    
    
    return 0;
}