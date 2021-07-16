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

    Gaussian(int _dim, double _weight, vector<double>& _mean, vector<double>& _cov){
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
    vector<double> mean;
    vector<double> cov;
    double weight;
    bool isNull = false;
};



const Eigen::Matrix3d v2t(vector<double>& vec);
const vector<double> t2v(const Eigen::Ref<Eigen::Matrix3d>& mtx);

// multiplication of two multinominal Gaussian distribution (assume diagonal covariance) 
Gaussian multGaussians(const Gaussian& g1, const Gaussian& g2);

// calculate pdf of multivariate normal distribution
double calcNormalPDF(vector<double> x, vector<double> mean, vector<double> cov, int dim);

void showMatrix(Eigen::Matrix3d _T);

vector<vector<int>> cartProduct(const vector<vector<int>>& v);

vector<Gaussian> exactSampling(vector<vector<Gaussian>>& mixtures, int dim, bool reparam, vector<string>& types, bool showMode);

vector<int> intersection(vector<int> &v1, vector<int> &v2);

void copy_vector(igraph_vector_t *v, vector<int>& vec);

void print_vector(igraph_vector_t *v, FILE *f);

void print_std_vector(vector<int>& vec);

void readCSV(string path, vector<vector<int>>& _idx_v, vector<vector<double>>& _odom);

double error_per(vector<double>& est, vector<double>& gt);

void getMode(vector<vector<double>>& means, vector<vector<double>>& covs, vector<double>& ws, int max_iter, int dim, int nb_bin);

void vec2pts(vector<vector<double>>& pts_in, vector<Point>& pts_out);

vector<double> calcDist(const vector<double>& p1, const vector<double>& p2);

void readCSV_MH(string path, vector<vector<int>>& _idx_v, vector<vector<vector<double>>>& _odom, int last_idx);

void readCSV_MH_GT(string path, vector<vector<double>>& _odom, int last_idx);

#endif
