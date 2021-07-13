#include "utils.h"
#include "factor_graph.h"
#include <time.h>
#include <omp.h>
#include <iostream>

#pragma warning( disable : 4100 )

vector<double> prior_cov = {0.0001, 0.0001, 0.0001};
vector<double> meas_cov = {0.05, 0.05, 0.01};
vector<double> odom_cov = {0.1, 0.1, 0.02};

int dim = 3;

clock_t start, finish;
double duration;

mutex m;

int nb_var = 1000;


int main(int argc, char** argv){
    


    string seq = "city10000";
    // read data
    string path = "../data/"+seq+"/slam.csv";

    vector<vector<int>> idxs_v;
    vector<vector<vector<double>>> odom;
    readCSV_MH(path, idxs_v, odom, nb_var-1);

    string path_gt = "../data/"+seq+"/gt.csv";
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
    
    cout<<"Read Data ..."<<endl;

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

    // odom
    for(int i=0; i<nb_data; i++){

        

        int idx_from = idxs_v[i][0];
        int idx_to = idxs_v[i][1];

        vector<vector<double>> o = odom[i];

        bool isLoop = (abs(idx_to - idx_from) > 1);

        if(i%1 == 0){
            //cout<<i<<" th line loading... / isLoop : "<<isLoop<<endl;
        }

        
        if(!isLoop){

            vector<Gaussian> z;
            vector<double> o_false = {0.0, 0.0, 0.0};

            if(o.size() == 2){
                z = {Gaussian(dim, 0.5, o[0], odom_cov), Gaussian(dim, 0.5, o[1], odom_cov)};
            }
            else{
                z = {Gaussian(dim, 1.0, o[0], odom_cov)};
            }

            int idx = graph.vars.size();  
            vector<double> pose = graph.vars[idx-1]->mean;

            Eigen::Matrix3d T = (v2t(pose)*v2t(o[0]));
            pose = move(t2v(T));

            graph.addVariable(idx, pose, odom_cov);
            graph.addFactor(graph.factors.size(), z, "between", idx_from, idx_to, m);
        }
        else{
            vector<Gaussian> z;
            vector<double> pose_from = graph.vars[idx_from]->mean;
            vector<double> pose_to = graph.vars[idx_to]->mean;

            Eigen::Matrix3d T_NULL = (v2t(pose_from).inverse()*v2t(pose_to));
            vector<double> o_NULL = move(t2v(T_NULL));

            if(o.size() == 2){

                z = {Gaussian(dim, 1/3, o[0], odom_cov), Gaussian(dim, 1/3, o[1], odom_cov), Gaussian(dim, 1/3, o_NULL, odom_cov)}; // z + null loop

            }
            else{
                z = {Gaussian(dim, 1/2, o[0], odom_cov), Gaussian(dim, 1/2, o_NULL, odom_cov)}; // z + null loop
            }

            graph.addFactor(graph.factors.size(), z, "loop", idx_from, idx_to, m);
          
            start = clock();
            for(int iter=0; iter < 1; iter++){
                graph.propagateMsgAll(false, m);
            }

            finish = clock();
            duration = (double)(finish - start)/CLOCKS_PER_SEC;

            graph.updatePoseAll();
            cout<<i<<" Update : time = "<<duration<<" err = "<<graph.getError(gt)<<endl;

        }

    }

    graph.getUpdateOrdetJT();

    cout<<"Initial err = "<<graph.getError(gt)<<endl;

    for(int t=0; t<1000; t++){
        start = clock();
        for(int iter=0; iter < 1; iter++){
            graph.propagateMsgAll(true, m);
        }

        finish = clock();
        duration = (double)(finish - start)/CLOCKS_PER_SEC;

        graph.updatePoseAll();
        if(t%1 == 0)
            cout<<t<<" Update : time = "<<duration<<" err = "<<graph.getError(gt)<<endl;
    }


    cout<<" err = "<<graph.getError(gt)<<endl;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // string seq = "seq5";
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

    // // odom
    // for(int i=0; i<nb_data; i++){

    //     int idx_from = idxs_v[i][0];
    //     int idx_to = idxs_v[i][1];

    //     vector<double> o = odom[i];

    //     bool isLoop = (abs(idx_to - idx_from) > 1);

        
    //     if(!isLoop){


    //         vector<Gaussian> z;
    //         vector<double> o_false = {0.0, 0.0, 0.0};

    //         if(i%300 == 0){
    //             z = {Gaussian(dim, 0.5, o, odom_cov), Gaussian(dim, 0.5, o_false, odom_cov)};
    //         }
    //         else{
    //             z = {Gaussian(dim, 1.0, o, odom_cov)};
    //         }

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

    //         vector<Gaussian> z = {Gaussian(dim, 0.5, o, meas_cov), Gaussian(dim, 0.5, o_NULL, meas_cov)}; // z + null loop
    //         //vector<Gaussian> z = {Gaussian(dim, 1.0, o, meas_cov)}; // z + null loop

    //         graph.addFactor(graph.factors.size(), z, "loop", idx_from, idx_to, m);
          
    //         start = clock();
    //         for(int iter=0; iter < 50; iter++){
    //             graph.propagateMsgAll(false, m);
    //         }

    //         finish = clock();
    //         duration = (double)(finish - start)/CLOCKS_PER_SEC;

    //         graph.updatePoseAll();
    //         cout<<"Update : time = "<<duration<<" err = "<<graph.getError(gt)<<endl;

    //     }

    // }


    // cout<<" err = "<<graph.getError(gt)<<endl;
}
