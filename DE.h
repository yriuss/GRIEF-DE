#include <stdio.h>
#include <iostream>

typedef std::vector<int> Vector;
typedef std::vector<Vector> IndMat;
typedef std::vector<IndMat> PopMat;

struct IndShape
{
    uint x;
    uint y;
};

struct PopShape
{
    uint x;
    uint y;
    uint z;
};



class Pop: public Individual{
public:
private:
    std::vector<int> shape(3);
}



class DE{
public:
    DE(int N_pop, IndShape ind_shape);

    void create_population();
    void generate
    int crossover();
    int mutation();
    int selection();
private:
    float cr;
    Pop population;

};