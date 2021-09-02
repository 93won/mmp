#include "utils.h"
#include "factor_graph.h"
#include <time.h>
#include <omp.h>
#include <iostream>
#include <chrono>

#pragma warning( disable : 4100 )

vector<double> prior_cov = {0.0001, 0.0001, 0.0001};
vector<double> meas_cov = {0.01, 0.01, 0.002};
vector<double> odom_cov = {0.05, 0.05, 0.01};

int dim = 3;

clock_t start, finish;
double duration;

mutex m;

int nb_var = 3500;

int s = 0;


int main(int argc, char** argv){
    

    string seq = "M3500";
    // read data
    string path = "../data/"+seq+"/data_"+to_string(s)+".csv";
    //string path = "../data/"+seq+"/slam.csv";

    vector<vector<int>> idxs_v;
    vector<vector<vector<double>>> odom;
    readCSV_MH(path, idxs_v, odom, nb_var-1);
    string path_gt = "../data/"+seq+"/gt2.csv";
    //string path_gt = "../data/"+seq+"/gt.csv";
    vector<vector<double>> gt;


    readCSV_MH_GT(path_gt, gt, nb_var-1);
    
    int nb_data = odom.size();

    cout<<nb_data<<endl;

    // initialize graph
    FactorGraph graph(dim);

    // first pose
    vector<double> initial = {0, 0, 0};
    vector<Gaussian> prior = {Gaussian(dim, 1.0, initial, prior_cov)};
    graph.addVariable(0, initial, prior_cov);
    graph.addFactor(0, prior, "prior", 0, 0, m);

    vector<vector<double>> poses;
    poses.push_back(initial);
    
    cout<<"Read Data ..."<<endl;

    vector<vector<double>> result;

    int to_ref = 0;

    for(int i=0; i<nb_data; i++){

     
        int idx_from = idxs_v[i][0];
        int idx_to = idxs_v[i][1];

        vector<vector<double>> o = odom[i];

        bool isLoop = (abs(idx_to - idx_from) > 1);
        
        if(!isLoop){

            vector<Gaussian> z;

            int idx = graph.vars.size();  
            vector<double> pose = graph.vars[idx-1]->mean;

            vector<double> _cov = odom_cov;


            Eigen::Matrix3d T;

            if(o.size() == 2){
                z = {Gaussian(dim, 1/2, o[0], _cov), Gaussian(dim, 1/2, o[1], _cov)};
                T = (v2t(pose)*v2t(o[1]));

            }
            else{
                z = {Gaussian(dim, 1, o[0], _cov)};
                T = (v2t(pose)*v2t(o[0]));
            }
            
            pose = move(t2v(T));
            poses.push_back(pose);
            

            graph.addVariable(idx, pose, odom_cov);
            graph.addFactor(graph.factors.size(), z, "between", idx_from, idx_to, m);
        }
        else{


            vector<Gaussian> z;
            vector<double> pose_from = graph.vars[idx_from]->mean;
            vector<double> pose_to = graph.vars[idx_to]->mean;

            vector<double> _cov = meas_cov;

            Eigen::Matrix3d T_NULL = (v2t(pose_from).inverse()*v2t(pose_to));
            vector<double> o_NULL = move(t2v(T_NULL));



            if(o.size() == 2){
                z = {Gaussian(dim, 1/2, o[0], _cov), Gaussian(dim, 1/2, o[1], _cov)}; 
                // Gaussian(dim, 1/3, o_NULL, _cov)}; // z + null loop
                // z = {Gaussian(dim, 1/2, o[0], _cov), Gaussian(dim, 1/2, o[1], _cov)}; // z + null loop
            }
            else{
                z = {Gaussian(dim, 1, o[0], _cov)};//, Gaussian(dim, 1/2, o_NULL, _cov)}; // z + null loop
            }



            graph.addFactor(graph.factors.size(), z, "loop", idx_from, idx_to, m);
            
            if(to_ref != idx_to){
                to_ref = idx_to;
                std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

                //graph.getUpdateOrdetJT();
                // omp_set_num_threads(5);
                // #pragma omp parallel for
                for(int iter=0; iter < 3; iter++){
                    graph.propagateMsgAll(false, m);

                }

                    graph.updatePoseAll();
                
                //duration = (double)(finish - start)/CLOCKS_PER_SEC;

                std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
                double _e = graph.getError(gt);
                cout<<idx_to<<" Update : time = "<<sec.count()<<" err = "<<_e<<" delta = "<<graph.mean_delta/graph.vars.size()<<endl;
                graph.mean_delta = 0.0;
                vector<double> temp_result = {(double)i, sec.count(), _e};

            }

        }
    }


    std::ofstream outFile("../result/"+seq+"/result.txt");
    for (const auto &e : result) outFile << to_string((int)e[0]) <<" "<<to_string(e[1])<<" "<<to_string(e[2]) << "\n";

    cout<<"../result/"+seq+"/pose"+to_string(s)+".txt"<<endl;
    std::ofstream outFile2("../result/"+seq+"/pose"+to_string(s)+".txt");
    for (int i=0; i<graph.vars.size(); i++) {
        vector<double> pose = graph.getVarPose(i);
        outFile2 << to_string(pose[0]) <<" "<<to_string(pose[1])<<" "<<to_string(pose[2]) << "\n";
    }
    std::ofstream outFile3("../result/"+seq+"/pose_before.txt");
    for (int i=0; i<poses.size(); i++) {
        vector<double> pose = poses[i];
        outFile2 << to_string(pose[0]) <<" "<<to_string(pose[1])<<" "<<to_string(pose[2]) << "\n";
    }
}
