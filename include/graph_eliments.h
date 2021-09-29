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

    Variable(const int _idx, const string _key, const std::vector<double> _mean, std::vector<double> _cov){
        this->idx = _idx;
        this->key = _key;
        this->mean = _mean;
        this->cov = _cov;
    };

    Variable(){};
    
    int idx;
    string key;
    std::vector<double> mean;
    std::vector<double> cov;
    std::vector<std::shared_ptr<Factor>> neighbor_factors; // neighbor factors
    std::vector<std::shared_ptr<Edge>> neighbor_edges;
    
    void addNeighbor(const std::shared_ptr<Edge>& edge, const std::shared_ptr<Factor>& factor){
        this->neighbor_edges.push_back(edge);
        this->neighbor_factors.push_back(factor);   
    }

    void update(const Gaussian &g){

        std::vector<double> delta = calcDist(this->mean, g.mean);

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
    
    Factor(const int _idx, const string _key, const std::vector<Gaussian> _zs, const string _type){
        this->idx = _idx;
        this->key = _key;
        this->zs = _zs;
        this->type = _type;
    };

    Factor(){};

    int idx;
    string key;
    std::vector<Gaussian> zs; // measurements
    string type; // prior or between or loop

    std::shared_ptr<Edge> from_edge; // from edge
    std::shared_ptr<Edge> to_edge; // to edge

    std::shared_ptr<Variable> from_var; // from var
    std::shared_ptr<Variable> to_var; // to var

    void setEdge(const std::shared_ptr<Variable>& _var,
                 const std::shared_ptr<Edge>& _edge, string _type){
        
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
         const std::shared_ptr<Variable>& _var,
         const std::shared_ptr<Factor>& _factor, string _type){

        this->idx = _idx_edge;
        this->key = _key_edge;

        this->var = _var;
        this->factor = _factor;

        this->type = _type;
        this->valid = true;

    };

    Edge(){};

    void showEdgeInfo(){
        cout<<"Edge "<<this->key<<endl;
        cout<<this->var->key <<" "<<this->factor->key<<endl;
    }

    // (var) -- (edge) -- (factor)

    bool valid = false;

    int idx;
    string key;
    string type;


    std::shared_ptr<Variable> var; // variable of edge
    std::shared_ptr<Factor> factor; // factor of edge
    std::shared_ptr<Message> msg_var_to_factor; // message from variable to factor
    std::shared_ptr<Message> msg_factor_to_var; // message from factor to variable

};


class Message{

public:
    Message(const string _type,
            const std::shared_ptr<Variable>& _var,
            const std::shared_ptr<Factor>& _factor){
        this->type = _type;
        this->var = _var;
        this->factor = _factor;
        this->valid = true;
    };

    Message(){};

    void setGaussians(std::vector<Gaussian>& _gs){
        this->gs = _gs;
    }

    bool valid = false;

    std::shared_ptr<Variable> var;
    std::shared_ptr<Factor> factor;

    std::vector<Gaussian> gs; // obejct -> move
    std::vector<double> ws;
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