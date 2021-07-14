#include "utils.h"
#include "factor_graph.h"
#include <time.h>
#include <omp.h>
#include <iostream>
#include <chrono>

#pragma warning( disable : 4100 )

vector<double> prior_cov = {0.0001, 0.0001, 0.0001};
vector<double> meas_cov = {0.05, 0.05, 0.01};
vector<double> odom_cov = {0.05, 0.05, 0.01};

int dim = 3;

clock_t start, finish;
double duration;

mutex m;

int nb_var = 3500;


int main(int argc, char** argv){
    

    string seq = "M3500";
    // read data
    string path = "../data/"+seq+"/data_1.csv";
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

            vector<double> _cov = odom_cov;

            Eigen::Matrix3d T_NULL = (v2t(pose_from).inverse()*v2t(pose_to));
            vector<double> o_NULL = move(t2v(T_NULL));

            if(o.size() == 2){
                //z = {Gaussian(dim, 1/3, o[0], _cov), Gaussian(dim, 1/3, o[1], _cov), Gaussian(dim, 1/3, o_NULL, _cov)}; // z + null loop
                z = {Gaussian(dim, 1/3, o[0], _cov), Gaussian(dim, 1/3, o[1], _cov), Gaussian(true)}; // z + null loop

            }
            else{
                z = {Gaussian(dim, 1/2, o[0], _cov), Gaussian(true)}; // z + null loop
            }


            graph.addFactor(graph.factors.size(), z, "loop", idx_from, idx_to, m);
          
            
            if(to_ref != idx_to){
                to_ref = idx_to;
                std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

                // graph.getUpdateOrdetJT();
                // omp_set_num_threads(5);
                // #pragma omp parallel for
                for(int iter=0; iter < 1; iter++){
                    graph.propagateMsgAll(false, m);
                }
                std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
                //duration = (double)(finish - start)/CLOCKS_PER_SEC;

                graph.updatePoseAll();
                double _e = graph.getError(gt);
                cout<<idx_to<<" Update : time = "<<sec.count()<<" err = "<<_e<<" delta = "<<graph.mean_delta/graph.vars.size()<<endl;
                graph.mean_delta = 0.0;
                vector<double> temp_result = {(double)i, sec.count(), _e};

            }

        }
    }
    
    // for(int i=0; i<500; i++){
    //     std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    //     for(int iter=0; iter < 5; iter++){
    //         graph.propagateMsgAll(false, m);
    //     }
    //     std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;

    //     graph.updatePoseAll();
    //     double _e = graph.getError(gt);
    //     cout<<3500<<" Update : time = "<<sec.count()<<" err = "<<_e<<endl;
    //     vector<double> temp_result = {(double)i, sec.count(), _e};

    // }

    //graph.getUpdateOrdetJT();
    
    // cout<<-1<<" Update : "<<" err = "<<graph.getError(gt)<<endl;


    // for(int i=0; i<10000; i++){

    //     std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
    //     graph.propagateMsgAll(false, m);
    //     std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
        
    //     if(i%10 == 0){
    //         graph.updatePoseAll();
    //         double _e = graph.getError(gt);
    //         cout<<i<<" Update : time = "<<sec.count()<<" err = "<<_e<<endl;
    //         vector<double> temp_result = {(double)i, sec.count(), _e};
    //     }
    // }

    std::ofstream outFile("../result/"+seq+"/result.txt");
    for (const auto &e : result) outFile << to_string((int)e[0]) <<" "<<to_string(e[1])<<" "<<to_string(e[2]) << "\n";

    std::ofstream outFile2("../result/"+seq+"/pose.txt");
    for (int i=0; i<graph.vars.size(); i++) {
        vector<double> pose = graph.getVarPose(i);
        outFile2 << to_string(pose[0]) <<" "<<to_string(pose[1])<<" "<<to_string(pose[2]) << "\n";
    }
    std::ofstream outFile3("../result/"+seq+"/pose_before.txt");
    for (int i=0; i<poses.size(); i++) {
        vector<double> pose = poses[i];
        outFile2 << to_string(pose[0]) <<" "<<to_string(pose[1])<<" "<<to_string(pose[2]) << "\n";
    }







    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // string seq = "M3500";
    // // read data
    // string path = "../data/"+seq+"/slam.csv";

    // vector<vector<int>> idxs_v;
    // vector<vector<double>> odom;


    // string path_gt = "../data/"+seq+"/gt.csv";
    // vector<vector<int>> idxs_v_gt;
    // vector<vector<double>> odom_gt;

    // readCSV(path, idxs_v, odom);
    // readCSV(path_gt, idxs_v_gt, odom_gt);
    // int nb_data = odom.size();

    // cout<<"Test"<<endl;

    // // initialize graph
    // FactorGraph graph(dim);

    // // first pose
    // vector<double> initial = {0, 0, 0};
    // vector<Gaussian> prior = {Gaussian(dim, 1.0, initial, prior_cov)};
    // graph.addVariable(0, initial, prior_cov);
    // graph.addFactor(0, prior, "prior", 0, 0, m);
    
    // cout<<"Read Data ..."<<endl;

    // vector<double> pose_gt = {0, 0, 0};

    // vector<vector<double>> gt;
    // gt.push_back(pose_gt);

    // /* KITTI GT */
    // for(int i=0; i<odom_gt.size(); i++){
    //     vector<double> o = odom_gt[i];
    //     Eigen::Matrix3d T = (v2t(pose_gt)*v2t(o));
    //     pose_gt = move(t2v(T));
    //     gt.push_back(pose_gt);
    // }

  
    // vector<vector<double>> result;
    // // vector<double> temp_result = {(double)i, duration, e};
    // // //         result.push_back(temp_result);
    // //  std::ofstream outFile("../result/"+seq+"/result.txt");
    // // // for (const auto &e : result) outFile << to_string((int)e[0]) <<" "<<to_string(e[1])<<" "<<to_string(e[2]) << "\n";

    // // odom
    // int to_ref = 0;
    // for(int i=0; i<nb_data; i++){

    //     int idx_from = idxs_v[i][0];
    //     int idx_to = idxs_v[i][1];

    //     vector<double> o = odom[i];

    //     bool isLoop = (abs(idx_to - idx_from) > 1);

    //     vector<vector<double>> result;
    //     if(!isLoop){


    //         vector<Gaussian> z;
    //         vector<double> o_false = {0.0, 0.0, 0.0};

    //         z = {Gaussian(dim, 1.0, o, odom_cov)};

    //         int idx = graph.vars.size();  
    //         vector<double> pose = graph.vars[idx-1]->mean;

    //         Eigen::Matrix3d T = (v2t(pose)*v2t(o));
    //         pose = move(t2v(T));

    //         graph.addVariable(idx, pose, odom_cov);
    //         graph.addFactor(graph.factors.size(), z, "between", idx_from, idx_to, m);

    //     }
    //     else{
    //         vector<double> pose_from = graph.vars[idx_from]->mean;
    //         vector<double> pose_to = graph.vars[idx_to]->mean;

    //         Eigen::Matrix3d T_NULL = (v2t(pose_from).inverse()*v2t(pose_to));
    //         vector<double> o_NULL = move(t2v(T_NULL));

    //         vector<Gaussian> z = {Gaussian(dim, 1, o, meas_cov), Gaussian(true)}; // z + null loop
    //         //vector<Gaussian> z = {Gaussian(dim, 1.0, o, meas_cov)}; // z + null loop

    //         graph.addFactor(graph.factors.size(), z, "loop", idx_from, idx_to, m);

    //         if(to_ref != idx_to){
    //             to_ref = idx_to;
    //             std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    //             omp_set_num_threads(5);
    //             #pragma omp parallel for
    //             for(int iter=0; iter < 5; iter++){
    //                 graph.propagateMsgAll(false, m);
    //             }
    //             std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
    //             //duration = (double)(finish - start)/CLOCKS_PER_SEC;

    //             graph.updatePoseAll();
    //             double _e = graph.getError(gt);
    //             cout<<idx_to<<" Update : time = "<<sec.count()<<" err = "<<_e<<endl;
    //             vector<double> temp_result = {(double)i, sec.count(), _e};

    //         }
    //     }
    // }

    // std::ofstream outFile3("../result/"+seq+"/result.txt");
    // for (const auto &e : result) outFile3 << to_string((int)e[0]) <<" "<<to_string(e[1])<<" "<<to_string(e[2]) << "\n";

    // std::ofstream outFile("../result/"+seq+"/gt.txt");
    // for (int i=0; i<gt.size(); i++) {
    //     vector<double> pose = gt[i];
    //     outFile << to_string(pose[0]) <<" "<<to_string(pose[1])<<" "<<to_string(pose[2]) << "\n";
    // }

    // std::ofstream outFile2("../result/"+seq+"/pose.txt");
    // for (int i=0; i<graph.vars.size(); i++) {
    //     vector<double> pose = graph.getVarPose(i);
    //     outFile2 << to_string(pose[0]) <<" "<<to_string(pose[1])<<" "<<to_string(pose[2]) << "\n";
    // }

}
