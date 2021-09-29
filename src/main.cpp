#include "utils.h"
#include "factor_graph.h"
#include <time.h>
#include <omp.h>
#include <iostream>
#include <chrono>

#pragma warning( disable : 4100 )

std::vector<double> prior_cov = {0.00001, 0.00001, 0.00001};
std::vector<double> meas_cov = {0.01, 0.01, 0.002};
std::vector<double> odom_cov = {0.05, 0.05, 0.01};
std::vector<double> odom_double_cov = {0.1, 0.1, 0.02};
std::vector<double> meas_double_cov = {0.02, 0.02, 0.004};

int dim = 3;

clock_t start, finish;
double duration;

std::mutex m;

int nb_var = 3500;

int jt_update_term = 100;

int idx_to_ref = 0;


int s = 10;

int main(int argc, char** argv){

    // /////////////////////////////////////////////////////////////////////////////////////////////////
    
    string seq = "M3500";
    // read data
    string path = "../data/"+seq+"/data_high_noise.csv";

    std::vector<std::vector<int>> idxs_v;
    std::vector<std::vector<std::vector<double>>> odom;
    readCSV_MH(path, idxs_v, odom, nb_var-1);
    string path_gt = "../data/"+seq+"/gt.csv";
    std::vector<std::vector<double>> gt;


    readCSV_MH_GT(path_gt, gt, nb_var-1);
    int nb_data = odom.size();

    cout<<nb_data<<endl;

    // initialize graph
    FactorGraph graph(dim);

    // first pose
    std::vector<double> initial = {0, 0, 0};
    std::vector<Gaussian> prior = {Gaussian(dim, 1.0, initial, prior_cov)};
    graph.addVariable(0, initial, prior_cov);
    graph.addFactor(0, prior, "prior", 0, 0, m);

    std::vector<std::vector<double>> poses;
    poses.push_back(initial);
    
    cout<<"Read Data ..."<<endl;

    std::vector<std::vector<double>> result;

    int to_ref = 0;
    int jt_update = 0;

    for(int i=0; i<nb_data; i++){
        

        int idx_from = idxs_v[i][0];
        int idx_to = idxs_v[i][1];

        
        std::vector<std::vector<double>> o = odom[i];

        bool isLoop = (abs(idx_to - idx_from) > 1);
        
        if(!isLoop){
            
            

            std::vector<Gaussian> z;

            int idx = graph.vars.size();  
            std::vector<double> pose = graph.vars[idx-1]->mean;

            std::vector<double> _cov = odom_cov;


            Eigen::Matrix3d T;

            if(o.size() == 2){
                z = {Gaussian(dim, 1/2, o[0], odom_double_cov), Gaussian(dim, 1/2, o[1], odom_double_cov)};
                T = (v2t(pose)*v2t(o[1]));

            }
            else{
                z = {Gaussian(dim, 1, o[0], odom_cov)};
                T = (v2t(pose)*v2t(o[0]));
            }
            
            pose = move(t2v(T));
            poses.push_back(pose);
            

            graph.addVariable(idx, pose, odom_cov);
            graph.addFactor(graph.factors.size(), z, "between", idx_from, idx_to, m);

            

        }
        else{

            std::vector<Gaussian> z;
            std::vector<double> pose_from = graph.vars[idx_from]->mean;
            std::vector<double> pose_to = graph.vars[idx_to]->mean;
            std::vector<double> _cov = odom_cov;

            Eigen::Matrix3d T_NULL = (v2t(pose_from).inverse()*v2t(pose_to));
            std::vector<double> o_NULL = move(t2v(T_NULL));

            if(o.size() == 2){
                z = {Gaussian(dim, 1/2, o[0], _cov), Gaussian(dim, 1/2, o[1], _cov)}; 
            }
            else{
                z = {Gaussian(dim, 1, o[0], _cov)};
            }

            graph.addFactor(graph.factors.size(), z, "loop", idx_from, idx_to, m);

            if(to_ref != idx_to){
                
                to_ref = idx_to;

                if((idx_to - jt_update) > jt_update_term)
                {
                    graph.getUpdateOrdetJT(m);
                    jt_update = idx_to;
                    std::cout<<"################################# Juntion tree updated #################################"<<std::endl;
                }
                std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
                
                for(int iter=0; iter < 3; iter++){
                    graph.propagateMsgAll(true, m);

                }

                graph.updatePoseAll();
                std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;

                double _e = graph.getError(gt);
                cout<<idx_to<<" Update : time = "<<sec.count()<<" err = "<<_e<<" delta = "<<graph.mean_delta/graph.vars.size()<<" "<<graph.vars.size()<<endl;
                graph.mean_delta = 0.0;
                std::vector<double> temp_result = {(double)idx_to, sec.count(), _e};
                result.push_back(temp_result);

            }

        }

        if(idx_to != idx_to_ref){
            idx_to_ref = idx_to;
            std::ofstream outFile_("../result_for_plot/"+to_string(idx_to)+".txt");
            for (int k=0; k<graph.vars.size(); k++) {
                std::vector<double> pose = graph.getVarPose(k);
                outFile_ << to_string(pose[0]) <<" "<<to_string(pose[1])<<" "<<to_string(pose[2]) << "\n";
            }
        }
    }

    /*std::ofstream outFile("../result/"+seq+"/result"+to_string(s)+".txt");
    for (const auto &e : result) outFile << to_string((int)e[0]) <<" "<<to_string(e[1])<<" "<<to_string(e[2]) << "\n";


    cout<<"NB vars : "<<graph.vars.size()<<endl;
    cout<<"../result/"+seq+"/pose"+to_string(s)+".txt"<<endl;
    std::ofstream outFile2("../result/"+seq+"/pose"+to_string(s)+".txt");
    for (int i=0; i<graph.vars.size(); i++) {
        std::vector<double> pose = graph.getVarPose(i);
        outFile2 << to_string(pose[0]) <<" "<<to_string(pose[1])<<" "<<to_string(pose[2]) << "\n";
    }
    std::ofstream outFile3("../result/"+seq+"/pose_before.txt");
    for (int i=0; i<poses.size(); i++) {
        std::vector<double> pose = poses[i];
        outFile2 << to_string(pose[0]) <<" "<<to_string(pose[1])<<" "<<to_string(pose[2]) << "\n";
    }*/
}
