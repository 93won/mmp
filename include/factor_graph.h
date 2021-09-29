#ifndef FACTOR_GRAPH_H
#define FACTOR_GRAPH_H

#include "utils.h"
#include "graph_eliments.h"
#include <assert.h>

#include <random>

class FactorGraph{
public:
    FactorGraph(int _dim){
        this->dim = _dim;

        igraph_vector_init(&e, 100000);
        
    };

    FactorGraph(){};

    // igraph
    igraph_t g;
    igraph_vector_t e;
    int nb_between_edge = 0;

    int dim;
    std::vector<std::shared_ptr<Variable>> vars;
    std::vector<std::shared_ptr<Factor>> factors;
    std::vector<std::shared_ptr<Edge>> edges;

    std::vector<int> jt_order; // update order of factors

    void addVariable(int _idx, std::vector<double>& initial, std::vector<double>& _cov);

    void addFactor(int _idx, const std::vector<Gaussian>& _zs, 
                   const string _type, 
                   const int _idx_from_var, 
                   const int _idx_to_var, mutex& m);
    
    void buildMsgFromVarToFactor(int _idx_factor, string _type, mutex& m);
    void buildMsgFromFactorToVar(int _idx_factor, string _type, mutex& m);
    
    void propagateMsg(int _idx_factor, mutex& m);
    void updateVar(int _idx_var);

    void getUpdateOrdetJT(std::mutex& m);

    void propagateMsgAll(bool JT, mutex& m);
    void updatePoseAll();

    // functions for debug
    std::vector<string> var_keys;
    std::vector<string> factor_keys;
    std::vector<string> edge_keys;

    double mean_delta = 0.0;

    std::vector<double> getVarPose(int _idx){
        return this->vars[_idx]->mean;
    }

    double getError(std::vector<std::vector<double>>& gt){

        int nb_v = this->vars.size();

        double error = 0.0;
        for(int i=0; i<nb_v; i++){
            double err = error_per(this->vars[i]->mean, gt[i]);
            error += err;
        }

        return error/((double)nb_v);
        
    }

    void showVarPose(int _idx){

        assert(this->vars.size() >= _idx);
        //cout<<"No var x"+to_string(_idx)<<endl;


        double x = this->vars[_idx]->mean[0];
        double y = this->vars[_idx]->mean[1];
        double a = this->vars[_idx]->mean[2];

        cout<<"Pose of "<<this->vars[_idx]->key<<" x : "+to_string(x)+" y :"+to_string(y)+" theta : "+to_string(a)<<endl;
    }

    void showVarPoseAll(){
        int nb_vars = this->vars.size();
        for(int i=0; i<nb_vars; i++){
            this->showVarPose(i);
        }
    }

    void showMessage(int _idx){

        bool valid_from = this->factors[_idx]->from_edge->valid;
        bool valid_to = this->factors[_idx]->to_edge->valid;
        if(valid_from){
            string var_key = this->factors[_idx]->from_edge->var->key;
            string factor_key = this->factors[_idx]->from_edge->factor->key;

            Gaussian v2f_g = this->factors[_idx]->from_edge->msg_var_to_factor->gs[0];
            Gaussian f2v_g = this->factors[_idx]->from_edge->msg_factor_to_var->gs[0];

            cout<<"Message from "<<var_key<<" to "<<factor_key<<" : x="<<v2f_g.mean[0]<<" y="<<v2f_g.mean[1]<<" h="<<v2f_g.mean[2]<<
            " // cov : x="<<v2f_g.cov[0]<<" y="<<v2f_g.cov[1]<<" z="<<v2f_g.cov[2]<<endl;
            cout<<"Message from "<<factor_key<<" to "<<var_key<<" : x="<<f2v_g.mean[0]<<" y="<<f2v_g.mean[1]<<" h="<<f2v_g.mean[2]<<
            " // cov : x="<<f2v_g.cov[0]<<" y="<<f2v_g.cov[1]<<" z="<<f2v_g.cov[2]<<endl;
        }

        if(valid_to){
            string var_key2 = this->factors[_idx]->to_edge->var->key;
            string factor_key2 = this->factors[_idx]->to_edge->factor->key;

            Gaussian v2f_g2 = this->factors[_idx]->to_edge->msg_var_to_factor->gs[0];
            Gaussian f2v_g2 = this->factors[_idx]->to_edge->msg_factor_to_var->gs[0];

            cout<<"Message from "<<var_key2<<" to "<<factor_key2<<" : x="<<v2f_g2.mean[0]<<" y="<<v2f_g2.mean[1]<<" h="<<v2f_g2.mean[2]<<
            " // cov : x="<<v2f_g2.cov[0]<<" y="<<v2f_g2.cov[1]<<" z="<<v2f_g2.cov[2]<<endl;
            cout<<"Message from "<<factor_key2<<" to "<<var_key2<<" : x="<<f2v_g2.mean[0]<<" y="<<f2v_g2.mean[1]<<" h="<<f2v_g2.mean[2]<<
            " // cov : x="<<f2v_g2.cov[0]<<" y="<<f2v_g2.cov[1]<<" z="<<f2v_g2.cov[2]<<endl;
        }

    }

    void showVarConnection(int _idx){
        assert(this->vars.size() >= _idx);
        //cout<<"No var x"+to_string(_idx)<<endl;


        cout<<"### Debug variable " + this->vars[_idx]->key + " ###"<<endl;
        cout<<"--- Neighbor factors  ---"<<endl;

        string factors_key = "";

        for(auto iter = this->vars[_idx]->neighbor_factors.begin(); iter != this->vars[_idx]->neighbor_factors.end(); iter++){

            factors_key += (*iter)->key + " ";
        }
        
        cout<<factors_key<<endl;

        cout<<"---  Neighbor Edges  ---"<<endl;

        string edges_key = "";

        for(auto iter = this->vars[_idx]->neighbor_edges.begin(); iter != this->vars[_idx]->neighbor_edges.end(); iter++){

            edges_key += (*iter)->key + " ";
        }

        cout<<edges_key<<endl;

        cout<<"#####   End Debug   #####"<<endl;
        cout<<endl;

    }

    

};



#endif