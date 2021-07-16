#include "factor_graph.h"
#include "utils.h"
#include <random>

void FactorGraph::addVariable(int _idx, vector<double>& initial, vector<double>& _cov){
    string key = "x" + to_string(_idx);
    this->vars.emplace_back(new Variable(_idx, key, initial, _cov));
}


void FactorGraph::addFactor(int _idx, const vector<Gaussian>& _zs, 
                            const string _type, 
                            const int _idx_from_var, 
                            const int _idx_to_var, mutex& m){

    bool isPrior = (_type == "prior");
    bool isBetween = (_type == "between");
    bool isLoop = (_type == "loop");

    string key_factor = "f"+to_string(_idx);
    int idx_factor = _idx;

    if(!isPrior){

        
        
        // for igraph
        this->nb_between_edge += 1;
        
        // from edge

        int idx_from_edge = this->edges.size();
        string key_from_edge = "e"+to_string(idx_from_edge);

        int idx_from_var = _idx_from_var;

        // make links
        this->factors.emplace_back(new Factor(idx_factor, key_factor, _zs, _type));
        this->edges.emplace_back(new Edge(idx_from_edge, key_from_edge, this->vars[idx_from_var], this->factors[idx_factor], "from"));
        this->factors[idx_factor]->setEdge(this->vars[idx_from_var], this->edges[idx_from_edge], "from");
        this->vars[idx_from_var]->addNeighbor(this->edges[idx_from_edge], this->factors[idx_factor]);

        // to edge

        int idx_to_edge = this->edges.size();
        string key_to_edge = "e"+to_string(idx_to_edge);

        int idx_to_var = _idx_to_var;

        // make links

        this->edges.emplace_back(new Edge(idx_to_edge, key_to_edge, this->vars[idx_to_var], this->factors[idx_factor], "to"));
        this->factors[idx_factor]->setEdge(this->vars[idx_to_var], this->edges[idx_to_edge], "to");
        this->vars[idx_to_var]->addNeighbor(this->edges[idx_to_edge], this->factors[idx_factor]);


        // propagate msg

        this->propagateMsg(idx_factor, m);

        // igraph
        igraph_vector_init(&e, 2);
        VECTOR(e)[0] = idx_from_var;
        VECTOR(e)[1] = idx_to_var;

        if(this->nb_between_edge == 1){
            igraph_create(&g, &e, 0, 1);
        }

        else{
            if(isBetween){
                igraph_add_vertices(&this->g, 1, 0);
            }
            igraph_add_edges(&g, &e, 0);
        }

    }
    else{
        // to edge

        int idx_to_edge = this->edges.size();
        string key_to_edge = "e"+to_string(idx_to_edge);

        int idx_to_var = _idx_to_var;

        // make links

        this->factors.emplace_back(new Factor(idx_factor, key_factor, _zs, _type));
        this->edges.emplace_back(new Edge(idx_to_edge, key_to_edge, this->vars[idx_to_var], this->factors[idx_factor], "to"));
        this->factors[idx_factor]->setEdge(this->vars[idx_to_var], this->edges[idx_to_edge], "to");
        this->vars[idx_to_var]->addNeighbor(this->edges[idx_to_edge], this->factors[idx_factor]);

        // propagate msg

        this->propagateMsg(idx_factor, m);
    }
}


void FactorGraph::buildMsgFromVarToFactor(int _idx_factor, string _type, mutex& m){

    shared_ptr<Edge> edge(new Edge());

    if(_type == "from"){
        edge = this->factors[_idx_factor]->from_edge;
    }
    else{
        edge = this->factors[_idx_factor]->to_edge;
    }

    //cout<<"Build msg from var "<<edge->var->key<<" to factor "<<edge->factor->key<<"  //  Type : "<<_type<<endl;

    bool isValid = (edge->valid);

    if(isValid){

        //shared_ptr<Message> msg(new Message(_type, edge->var, this->factors[_idx_factor]));
        shared_ptr<Message> msg(new Message(_type, edge->var, this->factors[_idx_factor]));

        bool isDirect = ((edge->var->neighbor_edges).size() == 1);

        // No other source except target factor -> direct message
        
        if(isDirect){

            vector<Gaussian> _gs = {Gaussian(this->dim, 1.0, edge->var->mean, edge->var->cov)};
            msg->setGaussians(_gs);

            edge->msg_var_to_factor = msg;
        }
        
        
        // More than one source except target factor -> multiply messages from other factors

        else{
            
            // pull messages from neighbor edges of variable of target edge exept target edge.

            //       * : target, // ~ : other source
            //
            //     (~edge)-(~msg)-(var)-(*msg)-(*edge)
            //                      |
            //                   (~msg)
            //                      |
            //                   (~edge)

            vector<shared_ptr<Message>> msgs;
            vector<string> types;

            for(auto iter = edge->var->neighbor_edges.begin(); iter != edge->var->neighbor_edges.end(); iter++){

                bool isOtherSource =((*iter)->key != edge->key);

                types.emplace_back((*iter)->factor->type);

                if(isOtherSource){
                    msgs.push_back((*iter)->msg_factor_to_var);
                }
            }


            if(msgs.size()==1){
                edge->msg_var_to_factor = msgs[0];
            }
            else{

                vector<vector<Gaussian>> mixtures(msgs.size());

                int nb_msgs = msgs.size();

                for(int i=0; i<nb_msgs; i++){
                    mixtures[i] = msgs[i]->gs;
                }
                vector<Gaussian> exact_product;

                if(mixtures.size() > 1){
                    exact_product = move(exactSampling(mixtures, this->dim, true, types, false));
                }
                else{
                    exact_product = mixtures[0];

                }
                msg->setGaussians(exact_product);
                //
                //{
                //lock_guard<mutex> lock_guard(m);
                edge->msg_var_to_factor = msg;
                //}
                //
                //cout<<"Exact sampling is done"<<endl;
            }
        }
        
    }
}


void FactorGraph::buildMsgFromFactorToVar(int _idx_factor, string _type, mutex& m){
    
    shared_ptr<Edge> edge(new Edge());
    shared_ptr<Edge> edge_source(new Edge());

    vector<Gaussian> zs = this->factors[_idx_factor]->zs; // copy

    // for(int i=0; i<zs.size(); i++){
    //     cout<<zs[i].mean[0]<<" "<<zs[i].mean[1]<<" "<<zs[i].mean[2]<<endl;
    // }
    if(_type == "from"){
        edge = this->factors[_idx_factor]->from_edge;
        edge_source = this->factors[_idx_factor]->to_edge;

        for(int i=0; i<zs.size(); i++){
            if(!zs[i].isNull){
                Eigen::Matrix3d T = v2t(zs[i].mean);
                T = T.inverse().eval();
                zs[i].mean = t2v(T);
            }
        }
    }

    else{
        edge = this->factors[_idx_factor]->to_edge;
        edge_source = this->factors[_idx_factor]->from_edge;
    }

    //cout<<"Build msg from factor "<<edge->factor->key<<" to var "<<edge->var->key<<"  //  Type : "<<_type<<endl;

    bool isValid = (edge->valid);

    // {
    //     lock_guard<mutex> lock_guard(m);

    if(isValid){
    
            shared_ptr<Message> msg;
            msg = edge_source->msg_var_to_factor;
            shared_ptr<Message> msg_new(new Message(_type, edge->var, this->factors[_idx_factor]));

            vector<Gaussian> gs;
            vector<Gaussian> gs_source = edge_source->msg_var_to_factor->gs;

            double w_sum = 0;
            
            for(auto& z : zs){
                if(!z.isNull){
                    for(auto& g : gs_source){
                        if(!g.isNull){
                                {
                                    //lock_guard<mutex> lock_guard(m);
                                    vector<double> mean_source = {g.mean[0], g.mean[1], g.mean[2]};
                                    Eigen::Matrix3d T = (v2t((g.mean)))*(v2t(z.mean));
                                    vector<double> mean = t2v(T);
                                    vector<double> cov = {g.cov[0]+z.cov[0], g.cov[1]+z.cov[1], g.cov[2]+z.cov[2]};
                                    gs.emplace_back(Gaussian(this->dim, g.weight*z.weight, mean, cov)); // weight mul is right?
                                    w_sum += g.weight*(z.weight);
                                }
                        }
                    } 
                }
                else{
                    
                }
                
            }
            

            for(auto& g : gs){
                g.weight /= w_sum;
            }

            msg_new->setGaussians(gs);
            edge->msg_factor_to_var = msg_new;
            //cout<<"DEBUG ~ Sum of weights : "<<w_sum<<" "<<w_sum_norm<<endl;
    }
    
}


void FactorGraph::propagateMsg(int _idx_factor, mutex& m){

    bool isPrior = (this->factors[_idx_factor]->type == "prior");
    bool isLoop = (this->factors[_idx_factor]->type == "loop");

    if(isPrior){

        // build message

        shared_ptr<Message> msg(new Message("to",
                                                 this->factors[_idx_factor]->to_edge->var,
                                                 this->factors[_idx_factor]));

        msg->setGaussians(this->factors[_idx_factor]->zs);

        //cout<<"Build msg from var "<<this->factors[_idx_factor]->to_edge->var->key<<" to factor "<<this->factors[_idx_factor]->to_edge->factor->key<<"  //  Type : "<<"to"<<endl;
        this->factors[_idx_factor]->to_edge->msg_var_to_factor = msg;

        //cout<<"Build msg from factor "<<this->factors[_idx_factor]->to_edge->factor->key<<" to var "<<this->factors[_idx_factor]->to_edge->var->key<<"  //  Type : "<<"to"<<endl;
        this->factors[_idx_factor]->to_edge->msg_factor_to_var = msg;

        // msg->showInfo();

    }

    else{
        {
            lock_guard<mutex> lock_guard(m);
            this->buildMsgFromVarToFactor(_idx_factor, "from", m);
            this->buildMsgFromFactorToVar(_idx_factor, "to", m);
            this->buildMsgFromVarToFactor(_idx_factor, "to", m);
            this->buildMsgFromFactorToVar(_idx_factor, "from", m);
        }
    }     
}

void FactorGraph::updateVar(int _idx_var){
    
    vector<vector<Gaussian>> mixtures(this->vars[_idx_var]->neighbor_edges.size());
    vector<string> types(this->vars[_idx_var]->neighbor_edges.size());
    int nb_neighbor_edges = this->vars[_idx_var]->neighbor_edges.size();

    for(int i=0; i<nb_neighbor_edges; i++){

        mixtures[i] = this->vars[_idx_var]->neighbor_edges[i]->msg_factor_to_var->gs;
        types[i] = this->vars[_idx_var]->neighbor_edges[i]->msg_factor_to_var->factor->type;
    }

    vector<Gaussian> products = exactSampling(mixtures, this->dim, false, types, false);
    vector<double> ws(products.size());

    for(int i=0; i<products.size(); i++){
        ws[i] = products[i].weight;
    }

    double w_sum = accumulate(ws.begin(), ws.end(), 0.0);

    for(int i=0; i<products.size(); i++){
        products[i].weight /= w_sum;
    }

    int max_idx = max_element(ws.begin(), ws.end()) - ws.begin();

    vector<double> temp_vec =  this->vars[_idx_var]->mean;

    this->vars[_idx_var]->update(products[max_idx]);

    //this->mean_delta += calcDist(temp_vec, this->vars[_idx_var]->mean);
    temp_vec = calcDist(temp_vec, this->vars[_idx_var]->mean);

    this->mean_delta += sqrt(pow(temp_vec[0], 2) + pow(temp_vec[1],2) + pow(temp_vec[2],2));

}

void FactorGraph::getUpdateOrdetJT(){

    // Step 1 : make graph chordal
    igraph_vector_t alpha; //
    igraph_vector_t alpham1; //
    igraph_bool_t chordal;
    igraph_vector_t fill_in; //
    igraph_t new_graph; //
    
    igraph_vector_init(&alpha, 1);
    igraph_vector_init(&alpham1, 1);
    igraph_vector_init(&fill_in, 1);

    igraph_maximum_cardinality_search(&this->g, &alpha, &alpham1);
    igraph_is_chordal(&this->g, &alpha, &alpham1, &chordal, &fill_in, &new_graph);

    igraph_vector_destroy(&alpha);
    igraph_vector_destroy(&alpham1);
    igraph_destroy(&new_graph);


    // fill in
    int nb_fill_in = igraph_vector_size(&fill_in);

    cout<<"Fill in : "<<nb_fill_in<<endl;
    vector<int> fill_in_vec(nb_fill_in);
    copy_vector(&fill_in, fill_in_vec);


    /* To Do */
    // Find shortest path and add edge in the factor graph

    // igraph_vector_t vetrices;
    // igraph_vector_t edges;
    // igraph_integer_t from;
    // igraph_integer_t to;
    // igraph_integer_t to;
    

    
    // Temporary solution
    igraph_vector_init(&e, nb_fill_in);

    for(int i=0; i<nb_fill_in/2; i++){
        VECTOR(e)[2*i] = fill_in_vec[2*i];
        VECTOR(e)[2*i+1] = fill_in_vec[2*i+1];
    }

    igraph_add_edges(&this->g, &e, 0); 
    igraph_vector_destroy(&e);


    cout<<"Fill in is done"<<endl;
    // Temporary solution


    //cout<<fill_in_vec[0]<<" "<<fill_in_vec[1]<<endl;

    igraph_vector_destroy(&fill_in);

    // check fill in -- later
    // print_vector(&fill_in, stdout);

    // Step 2 : build clique graph
    igraph_vector_ptr_t res; //
    igraph_vector_ptr_init(&res, 1);
    igraph_integer_t min_size = 2;
    igraph_integer_t max_size = 10000;

    igraph_maximal_cliques(&this->g, &res, min_size, max_size);
    int n = igraph_vector_ptr_size(&res);

    vector<vector<int>> v_idxs_cliques(n);

    //cout<<"###################"<<endl;
    
    for (int i = 0; i < n; i++) {
        igraph_vector_t* v2 = (igraph_vector_t*) igraph_vector_ptr_e(&res, i); //
        int m = igraph_vector_size(v2);
        vector<int> v_idxs_clique(m);

        for (int j = 0; j < m; j++) {
            v_idxs_clique[j] = (int)VECTOR(*v2)[j];
        }

        // sorting for find intersection
        sort(v_idxs_clique.begin(), v_idxs_clique.end());
        v_idxs_cliques[i] = v_idxs_clique;

        igraph_vector_destroy(v2);
        igraph_free(v2);

        //print_std_vector(v_idxs_clique);
    }

    // for(int i=0; i<v_idxs_cliques.size(); i++){
    //     print_std_vector(v_idxs_cliques[i]);
    // }

    //cout<<"###################"<<endl;

    igraph_vector_ptr_destroy(&res); 

    vector<vector<int>> idx_pairs_cliques;
    vector<vector<int>> seperators_each_pair;

    for(int i=0; i<n; i++){
        for(int j=i+1; j<n; j++){
            if(i!=j){
                vector<int> seperator = intersection(v_idxs_cliques[i], v_idxs_cliques[j]);
                if(seperator.size()>0){
                    vector<int> pair = {i, j};
                    idx_pairs_cliques.push_back(move(pair));
                    seperators_each_pair.push_back(move(seperator));
                }
            }
        }
    }

    // clique graph
    igraph_t g_clique;
    igraph_vector_t e_clique;
    igraph_vector_t e_weight;
    igraph_vector_t idx_m_spanning_tree_edge;

    int nb_clique_edges = idx_pairs_cliques.size();


    igraph_vector_init(&e_clique, 2*nb_clique_edges);
    igraph_vector_init(&e_weight, nb_clique_edges);
    igraph_vector_init(&idx_m_spanning_tree_edge, 1);

    for(int i=0; i<nb_clique_edges; i++){
        VECTOR(e_clique)[2*i] = idx_pairs_cliques[i][0];
        VECTOR(e_clique)[2*i+1] = idx_pairs_cliques[i][1];
        VECTOR(e_weight)[i] = (1.0)/((double)seperators_each_pair[i].size());
    }
    igraph_create(&g_clique, &e_clique, 0, 1);

    igraph_minimum_spanning_tree(&g_clique, &idx_m_spanning_tree_edge, &e_weight);

    //print_vector(&idx_m_spanning_tree_edge, stdout);

    int l = igraph_vector_size(&idx_m_spanning_tree_edge);

    vector<vector<int>> junction_edge(l);

    igraph_vector_t e_tree;
    igraph_vector_init(&e_tree, 2*l);

    for (int i = 0; i < l; i++) {
        int idx_edge = (long int)VECTOR(idx_m_spanning_tree_edge)[(long int)i];
        junction_edge[i] = idx_pairs_cliques[idx_edge];
        VECTOR(e_tree)[2*i] = idx_pairs_cliques[idx_edge][0];
        VECTOR(e_tree)[2*i+1] = idx_pairs_cliques[idx_edge][1];

        //cout<<idx_pairs_cliques[idx_edge][0]<<" "<<idx_pairs_cliques[idx_edge][1]<<endl;
    }

    /*
    
    Let G any BN minimal I-map for H. Then G is necessarily chordal.
    Every undirected chordal graph H has a clique tree T
    Let H be a chordal MN(H). Then there exists a BN G such that I(H) = I(G)
    
    */

    // for topological sort
    igraph_t g_tree;
    igraph_bool_t isDAG;
    igraph_vector_t order;

    igraph_vector_init(&order, n);

    igraph_create(&g_tree, &e_tree, 0, 1);

    igraph_is_dag(&g_tree, &isDAG);

    igraph_topological_sorting(&g_tree, &order, IGRAPH_IN);

    vector<int> vec_order(n);
    copy_vector(&order, vec_order);

    igraph_vector_destroy(&e_clique);
    igraph_vector_destroy(&e_weight);
    igraph_vector_destroy(&idx_m_spanning_tree_edge);
    igraph_vector_destroy(&e_tree);
    igraph_vector_destroy(&order);
    igraph_destroy(&g_tree);
    igraph_destroy(&g_clique);

    jt_order.clear();

    vector<bool> used(this->factors.size(), false);

    for(int i=0; i<n; i++){
        vector<int> idx_var = v_idxs_cliques[vec_order[i]];
        vector<int> idx_factor;
        for(int j=0; j<idx_var.size(); j++){

            for(auto& nf : this->vars[idx_var[j]]->neighbor_factors)
            {   
                //cout<<!count(idx_factor.begin(), idx_factor.end(), nf->idx)<<endl;
                if(!count(idx_factor.begin(), idx_factor.end(), nf->idx)){
                    //idx_factor.emplace_back(nf->idx);
                    if(!used[nf->idx]){
                        idx_factor.emplace_back(nf->idx);
                        used[nf->idx] = true;
                    }
                    //cout<<"Addition : "<<nf->idx<<endl;
                }
            }
        }
        this->jt_order.push_back(move(idx_factor));
    }
}

void FactorGraph::propagateMsgAll(bool JT, mutex& m){
    if(JT){
        for(auto& clique : this->jt_order){
            for(auto& idx_factor : clique){
                this->propagateMsg(idx_factor, m);
            }
        }
    }

    else{
        int nb_factors = this->factors.size();

        auto randeng = std::default_random_engine {};
        vector<int> order(nb_factors);
        iota(order.begin(), order.end(), 0);
        shuffle(order.begin(), order.end(), randeng);

        // omp_set_num_threads(20);
        // #pragma omp parallel for
        for(int i=0; i<nb_factors; i++){
            //cout<<i<<endl;
            this->propagateMsg(order[i], m);           
            
            // {
            // lock_guard<mutex> lock_guard(m);
            // }
        }
    }
}

void FactorGraph::updatePoseAll(){
    int nb_vars = this->vars.size();
    for(int i=0; i<nb_vars; i++){
            this->updateVar(i);
    }
}