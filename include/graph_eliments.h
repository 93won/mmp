#ifndef GRAPH_ELIMENTS_H
#define GRAPH_ELIMENTS_H

#include "utils.h"

class Edge;
class Factor;
class Variable;
class Message;


// (from_var) --(from_edge)-- (factor) --(to_edge)-- (to_var)

class Variable{

public:

    Variable(const int _idx, const string _key, const vector<double> _mean, vector<double> _cov){
        this->idx = _idx;
        this->key = _key;
        this->mean = _mean;
        this->cov = _cov;
    };

    Variable(){};
    
    int idx;
    string key;
    vector<double> mean;
    vector<double> cov;
    vector<shared_ptr<Factor>> neighbor_factors; // neighbor factors
    vector<shared_ptr<Edge>> neighbor_edges;
    
    void addNeighbor(const shared_ptr<Edge>& edge, const shared_ptr<Factor>& factor){
        this->neighbor_edges.push_back(edge);
        this->neighbor_factors.push_back(factor);   
    }

    void update(const Gaussian &g){

        vector<double> delta = calcDist(this->mean, g.mean);

        double trans = sqrt(pow(delta[0], 2) + pow(delta[1], 2));
        double rot = delta[2];

        if(trans > 0.0 || rot > 0.0){
            this->mean = g.mean;
            this->cov = g.cov;
        }
    }
};

class Factor{
public:
    
    Factor(const int _idx, const string _key, const vector<Gaussian> _zs, const string _type){
        this->idx = _idx;
        this->key = _key;
        this->zs = _zs;
        this->type = _type;
    };

    Factor(){};

    int idx;
    string key;
    vector<Gaussian> zs; // measurements
    string type; // prior or between or loop

    shared_ptr<Edge> from_edge; // from edge
    shared_ptr<Edge> to_edge; // to edge

    shared_ptr<Variable> from_var; // from var
    shared_ptr<Variable> to_var; // to var

    void setEdge(const shared_ptr<Variable>& _var,
                 const shared_ptr<Edge>& _edge, string _type){
        
        bool isFromSet(_type == "from");

        if(isFromSet){
            this->from_var = _var;
            this->from_edge = _edge;
        }
        else{
            this->to_var = _var;
            this->to_edge = _edge;
        }
    }
};

class Edge{
public:

    Edge(const int _idx_edge, const string _key_edge,
         const shared_ptr<Variable>& _var,
         const shared_ptr<Factor>& _factor, string _type){

        this->idx = _idx_edge;
        this->key = _key_edge;

        this->var = _var;
        this->factor = _factor;

        this->type = _type;
        this->valid = true;

    };

    Edge(){};

    // (var) -- (edge) -- (factor)

    bool valid = false;

    int idx;
    string key;
    string type;


    shared_ptr<Variable> var; // variable of edge
    shared_ptr<Factor> factor; // factor of edge
    shared_ptr<Message> msg_var_to_factor; // message from variable to factor
    shared_ptr<Message> msg_factor_to_var; // message from factor to variable

};


class Message{

public:
    Message(const string _type,
            const shared_ptr<Variable>& _var,
            const shared_ptr<Factor>& _factor){
        this->type = _type;
        this->var = _var;
        this->factor = _factor;
        this->valid = true;
    };

    Message(){};

    void setGaussians(vector<Gaussian>& _gs){
        this->gs = _gs;
    }

    bool valid = false;

    shared_ptr<Variable> var;
    shared_ptr<Factor> factor;

    vector<Gaussian> gs; // obejct -> move
    vector<double> ws;
    string type;

    void showInfo(){

        cout<<"Debug Message between "<<this->var->key<<" and "<<
        this->factor->key<<" // Type : "+type<<endl;

        for(int i=0; i<gs.size(); i++){
            cout<<"X : "<<this->gs[i].mean[0]
            <<" Y : "<<this->gs[i].mean[1]<<" H : "<<this->gs[i].mean[2]<<endl;
        }

    }
  

};

#endif