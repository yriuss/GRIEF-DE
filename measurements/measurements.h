#ifndef MEASUREMENTS_H_INCLUDED
#define MEASUREMENTS_H_INCLUDED

#include <random>
#include <stdio.h>
#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <unistd.h>
#include <fstream>
#include <cmath>
#include <fstream>
#include <sys/stat.h>
#include <fstream>
#include<math.h>
#define SAVE true

//namespace plt = matplotlibcpp;
class Measurements{
public:
    Eigen::MatrixXd mean(std::vector<Eigen::MatrixXd> x, bool save = false);
    double std_dev(std::vector<Eigen::MatrixXd> x, bool save = false);
    void reset();
    void append_fit(std::vector<std::vector<int>> all_fit);
    void make_boxes();
    void save_data(std::vector<float> y, const std::string &dataset, const std::string &exp, Eigen::MatrixXd best_individual, int count1, int count2, int count3, int count4, std::vector<std::vector<double>> F, int repair_counter);
private:
    std::string CURRENT_DIR = get_current_dir_name();
    std::vector<Eigen::MatrixXd> means;
    std::vector<double> std_devs;
    std::vector<std::vector<double>> boxes;
    std::vector<std::vector<std::vector<int>>> all_fits_ind;
    
    int gen = 1;
};
#endif