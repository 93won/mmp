#ifndef UTILS_H
#define UTILS_H

#include <cmath>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <memory>
#include <algorithm>
#include <numeric>
#include <igraph.h>
#include <istream>
#include <iterator> 
#include <boost/tokenizer.hpp>
#include <fstream>
#include "dbscan.h"
#include <sstream>
#include <chrono>
#include <stdlib.h>
#include <omp.h>
#include <mutex>
#include <random>

#include "dbscan.h"

using namespace std;




//auto rng2 = default_random_engine {};

class Gaussian{
public:

    Gaussian(bool isNull){
        //cout<<"Null Gaussian"<<endl;
        this->isNull = true;
    };

    Gaussian(int _dim, double _weight, std::vector<double>& _mean, std::vector<double>& _cov){
        this->dim = _dim;
        this->weight = _weight;
        this->mean = _mean;
        this->cov = _cov;
        this->isNull = false;
    };   

    Gaussian(){};

    void showInfo(){

        if(isNull){
            cout<<"Null Gaussian"<<endl;
        }
        else{
            cout<<"x = "<<this->mean[0]<<" y = "<<this->mean[1]<<" h = "<<this->mean[2]<<endl;
        }

    }

    int dim;
    std::vector<double> mean;
    std::vector<double> cov;
    double weight;
    bool isNull = false;
};



const Eigen::Matrix3d v2t(std::vector<double>& vec);
const std::vector<double> t2v(const Eigen::Ref<Eigen::Matrix3d>& mtx);

// multiplication of two multinominal Gaussian distribution (assume diagonal covariance) 
Gaussian multGaussians(const Gaussian& g1, const Gaussian& g2);

// calculate pdf of multivariate normal distribution
double calcNormalPDF(std::vector<double> x, std::vector<double> mean, std::vector<double> cov, int dim);

void showMatrix(Eigen::Matrix3d _T);

std::vector<std::vector<int>> cartProduct(const std::vector<std::vector<int>>& v);

std::vector<Gaussian> exactSampling(std::vector<std::vector<Gaussian>>& mixtures, int dim, bool reparam, std::vector<string>& types, bool showMode);

void intersection(std::vector<int> &v1, std::vector<int> &v2, std::vector<int> &v3);

void copy_vector(igraph_vector_t *v, std::vector<int>& vec);

void print_vector(igraph_vector_t *v, FILE *f);

void print_std_vector(std::vector<int>& vec);

void readCSV(string path, std::vector<std::vector<int>>& _idx_v, std::vector<std::vector<double>>& _odom);

double error_per(std::vector<double>& est, std::vector<double>& gt);

void getMode(std::vector<std::vector<double>>& means, std::vector<std::vector<double>>& covs, std::vector<double>& ws, int max_iter, int dim, int nb_bin);

void vec2pts(std::vector<std::vector<double>>& pts_in, std::vector<Point>& pts_out);

std::vector<double> calcDist(const std::vector<double>& p1, const std::vector<double>& p2);

void readCSV_MH(string path, std::vector<std::vector<int>>& _idx_v, std::vector<std::vector<std::vector<double>>>& _odom, int last_idx);

void readCSV_MH_GT(string path, std::vector<std::vector<double>>& _odom, int last_idx);

#endif
