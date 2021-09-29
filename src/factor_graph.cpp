#include "factor_graph.h"
#include "utils.h"
#include <random>

void FactorGraph::addVariable(int _idx, std::vector<double>& initial, std::vector<double>& _cov){
    string key = "x" + to_string(_idx);
    this->vars.emplace_back(new Variable(_idx, key, initial, _cov));
}


void FactorGraph::addFactor(int _idx, const std::vector<Gaussian>& _zs, 
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

        //cout<<"Edge addition "<<idx_from_var<<" "<<idx_to_var<<endl;

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

    std::shared_ptr<Edge> edge(new Edge());

    if(_type == "from"){
        edge = this->factors[_idx_factor]->from_edge;
    }
    else{
        edge = this->factors[_idx_factor]->to_edge;
    }

    //cout<<"Build msg from var "<<edge->var->key<<" to factor "<<edge->factor->key<<"  //  Type : "<<_type<<endl;

    bool isValid = (edge->valid);

    if(isValid){

        //std::shared_ptr<Message> msg(new Message(_type, edge->var, this->factors[_idx_factor]));
        std::shared_ptr<Message> msg(new Message(_type, edge->var, this->factors[_idx_factor]));

        bool isDirect = ((edge->var->neighbor_edges).size() == 1);

        // No other source except target factor -> direct message
        
        if(isDirect){
            std::vector<Gaussian> _gs = {Gaussian(this->dim, 1.0, edge->var->mean, edge->var->cov)};
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

            std::vector<std::shared_ptr<Message>> msgs;
            std::vector<string> types;

            //cout<<"#################  "<<edge->var->key<<"  "<<edge->factor->key<<endl;
            for(auto iter = edge->var->neighbor_edges.begin(); iter != edge->var->neighbor_edges.end(); iter++){

                bool isOtherSource =((*iter)->key != edge->key);

                if(isOtherSource){
                    types.emplace_back((*iter)->factor->type);
                    //cout<<(*iter)->factor->key<<" "<<(*iter)->factor->type<<endl;
                    // if(edge->var->key=="x2" and edge->factor->key=="f2"){
                    //     cout<<(*iter)->factor->key<<" "<<(*iter)->factor->type<<endl;
                    //     cout<<"from edge connecting "<<(*iter)->factor->key<<" and "<<(*iter)->var->key<<endl;
                    // }
                    std::shared_ptr<Message> m_temp(new Message(_type, edge->var, edge->factor));
                    m_temp->setGaussians((*iter)->msg_factor_to_var->gs);
                    msgs.push_back(move(m_temp));



                }
            }

            if(msgs.size()==1){
                edge->msg_var_to_factor = msgs[0];
            }
            else{

                std::vector<std::vector<Gaussian>> mixtures(msgs.size());

                int nb_msgs = msgs.size();

                for(int i=0; i<nb_msgs; i++){
                    mixtures[i] = msgs[i]->gs;
                }
                std::vector<Gaussian> exact_product;

                if(mixtures.size() > 1){
                    
                    // if(edge->var->key=="x2" and edge->factor->key=="f2"){

                    //     //cout<<types[0]<<" "<<types[1]<<" "<<types[2]<<endl;

                    // exact_product = move(exactSampling(mixtures, this->dim, true, types, true));

                    // }
                    // else{

                    // exact_product = move(exactSampling(mixtures, this->dim, true, types, false));
                    // }

                    exact_product = move(exactSampling(mixtures, this->dim, true, types, false));
                }
                else{
                    exact_product = mixtures[0];

                }
                msg->setGaussians(exact_product);
                edge->msg_var_to_factor = msg;
            } 
        }
        
    }

   
}


void FactorGraph::buildMsgFromFactorToVar(int _idx_factor, string _type, mutex& m){
    
    std::shared_ptr<Edge> edge(new Edge());
    std::shared_ptr<Edge> edge_source(new Edge());

    std::vector<Gaussian> zs = this->factors[_idx_factor]->zs; // copy

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
    
            std::shared_ptr<Message> msg;
            msg = edge_source->msg_var_to_factor;
            std::shared_ptr<Message> msg_new(new Message(_type, edge->var, this->factors[_idx_factor]));

            std::vector<Gaussian> gs;
            std::vector<Gaussian> gs_source = edge_source->msg_var_to_factor->gs;


            double w_sum = 0;
            
            for(auto& z : zs){
                if(!z.isNull){
                    for(auto& g : gs_source){
                        if(!g.isNull){
                                {
                                    //lock_guard<mutex> lock_guard(m);
                                    std::vector<double> mean_source = {g.mean[0], g.mean[1], g.mean[2]};
                                    Eigen::Matrix3d T = (v2t((g.mean)))*(v2t(z.mean));
                                    std::vector<double> mean = t2v(T);
                                    //std::vector<double> cov = {g.cov[0]+z.cov[0], g.cov[1]+z.cov[1], g.cov[2]+z.cov[2]};
                                    std::vector<double> cov = {z.cov[0], z.cov[1], z.cov[2]};
                                    //gs.emplace_back(Gaussian(this->dim, g.weight*z.weight, mean, cov)); // weight mul is right?
                                    gs.emplace_back(Gaussian(this->dim, z.weight, mean, cov)); // weight mul is right?
                                    //w_sum += g.weight*(z.weight);
                                    w_sum += z.weight;

                                    // cout<<"Mean source : "<<g.mean[0]<<endl;
                                    // cout<<"Mean measurement : "<<z.mean[0]<<endl;
                                    // cout<<"Mean result : "<<mean[0]<<endl;
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

    //cout<<"Factor "<<_idx_factor<<" Propagation"<<endl;

    bool isPrior = (this->factors[_idx_factor]->type == "prior");
    bool isLoop = (this->factors[_idx_factor]->type == "loop");

    if(isPrior){

        // build message

        std::shared_ptr<Message> msg(new Message("to",
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
            //lock_guard<mutex> lock_guard(m);
            this->buildMsgFromVarToFactor(_idx_factor, "from", m);
            this->buildMsgFromFactorToVar(_idx_factor, "to", m);
            this->buildMsgFromVarToFactor(_idx_factor, "to", m);
            this->buildMsgFromFactorToVar(_idx_factor, "from", m);
        }
    }     
}

void FactorGraph::updateVar(int _idx_var){

    std::vector<std::vector<Gaussian>> mixtures(this->vars[_idx_var]->neighbor_edges.size());
    std::vector<string> types(this->vars[_idx_var]->neighbor_edges.size());
    int nb_neighbor_edges = this->vars[_idx_var]->neighbor_edges.size();

    for(int i=0; i<nb_neighbor_edges; i++){

        mixtures[i] = this->vars[_idx_var]->neighbor_edges[i]->msg_factor_to_var->gs;
        types[i] = this->vars[_idx_var]->neighbor_edges[i]->msg_factor_to_var->factor->type;
    }
    std::vector<Gaussian> products;
    // if(_idx_var == 1)
    //     products = exactSampling(mixtures, this->dim, false, types, true);

    // else
    products = exactSampling(mixtures, this->dim, false, types, false);

    std::vector<double> ws(products.size());

    for(int i=0; i<products.size(); i++){
        ws[i] = products[i].weight;
    }

    double w_sum = accumulate(ws.begin(), ws.end(), 0.0);

    for(int i=0; i<products.size(); i++){
        products[i].weight /= w_sum;
    }

    int max_idx = max_element(ws.begin(), ws.end()) - ws.begin();

    std::vector<double> temp_vec =  this->vars[_idx_var]->mean;

    this->vars[_idx_var]->update(products[max_idx]);

    //this->mean_delta += calcDist(temp_vec, this->vars[_idx_var]->mean);
    temp_vec = calcDist(temp_vec, this->vars[_idx_var]->mean);

    this->mean_delta += sqrt(pow(temp_vec[0], 2) + pow(temp_vec[1],2) + pow(temp_vec[2],2));
}

void FactorGraph::getUpdateOrdetJT(std::mutex& m){

    std::vector<std::chrono::system_clock::time_point> ts;

    ts.push_back(std::chrono::system_clock::now());
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

    ts.push_back(std::chrono::system_clock::now());


    // fill in
    // int nb_fill_in = igraph_vector_size(&fill_in);
    // std::vector<int> fill_in_vec(nb_fill_in);
    // copy_vector(&fill_in, fill_in_vec);


    /* To Do */
    // Find shortest path and add edge in the factor graph

    // igraph_vector_t vetrices;
    // igraph_vector_t edges;
    // igraph_integer_t from;
    // igraph_integer_t to;
    // igraph_integer_t to;
    

    
    // Temporary solution
    // igraph_vector_init(&e, nb_fill_in);

    // int nb_fill_in_valid = 0;
    
    // for(int i=0; i<nb_fill_in/2; i++){
    //     if(fill_in_vec[2*i] != fill_in_vec[2*i+1]){
    //     VECTOR(e)[2*i] = fill_in_vec[2*i];
    //     VECTOR(e)[2*i+1] = fill_in_vec[2*i+1];
    //     nb_fill_in_valid += 1;
    //     }
    // }

    igraph_add_edges(&this->g, &fill_in, 0); 
    //igraph_vector_destroy(&e);

    igraph_vector_destroy(&fill_in);

    //cout<<"Fill in : "<<nb_fill_in_valid<<endl;
    // Temporary solution



    //// check all edge
    // igraph_vector_t edges;
    // igraph_vector_init(&edges, 22);
    // igraph_get_edgelist(&this->g, &edges, 0);

    // int nb_edges = igraph_vector_size(&edges);
    // std::vector<int> all_edges(nb_edges);
    // copy_vector(&edges, all_edges);

    // igraph_vector_destroy(&edges);

    // for(int i=0; i<all_edges.size()/2; i++){
    //     std::cout<<all_edges[2*i]<<" -- "<<all_edges[2*i+1]<<endl;
    // }
    // std::cout<<all_edges.size()/2<<endl;
    

    // check fill in -- later
    // print_vector(&fill_in, stdout);

    ts.push_back(std::chrono::system_clock::now());
    // Step 2 : build clique graph
    igraph_vector_ptr_t res; //
    igraph_vector_ptr_init(&res, 1);
    igraph_integer_t min_size = 2;
    igraph_integer_t max_size = 10000;

    igraph_maximal_cliques(&this->g, &res, min_size, max_size);

    ts.push_back(std::chrono::system_clock::now());
    int n = igraph_vector_ptr_size(&res);

    std::vector<std::vector<int>> v_idxs_cliques(n);

    //cout<<"###################"<<endl;
    
    for (int i = 0; i < n; i++) {
        igraph_vector_t* v2 = (igraph_vector_t*) igraph_vector_ptr_e(&res, i); //
        int m = igraph_vector_size(v2);
        std::vector<int> v_idxs_clique(m);

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

    // //// check cliques -- debug step 1
    // for(int i=0; i<v_idxs_cliques.size(); i++){
    //     std::cout<<"Clique "<<i<<" : ";
    //     print_std_vector(v_idxs_cliques[i]);
    // }

    //cout<<"###################"<<endl;

    ts.push_back(std::chrono::system_clock::now());
    igraph_vector_ptr_destroy(&res); 

    std::vector<std::vector<int>> idx_pairs_cliques;
    std::vector<std::vector<int>> separators_each_pair;
    
    //std::cout<<"NB cliques : "<<n<<std::endl;

    for(int i=0; i<n; i++){
        for(int j=i+1; j<n; j++){
            if(i!=j){
                        
                //lock_guard<mutex> lock_guard(m);
                std::vector<int> separator;
                intersection(v_idxs_cliques[i], v_idxs_cliques[j], separator);
                
                if(separator.size()>0){
                    std::vector<int> pair = {i, j};

                    //// check separator -- debug step 2
                    // cout<<"## Seperator between clique C"<<pair[0]<<" and C"<<pair[1]<<" : ";
                    // print_std_vector(separator);
                    {
                    //lock_guard<mutex> lock_guard(m);
                    idx_pairs_cliques.push_back(move(pair));
                    separators_each_pair.push_back(move(separator));
                    }
                }
            }
        }
    }

    ts.push_back(std::chrono::system_clock::now());
    // build clique graph
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
        VECTOR(e_weight)[i] = (1.0)/((double)separators_each_pair[i].size());
    }
    igraph_create(&g_clique, &e_clique, 0, 1);

    ts.push_back(std::chrono::system_clock::now());
    // // check all edge
    // igraph_vector_t edges;
    // igraph_vector_init(&edges, 1);
    // igraph_get_edgelist(&g_clique, &edges, 0);

    // int nb_edges = igraph_vector_size(&edges);
    // std::vector<int> all_edges(nb_edges);
    // copy_vector(&edges, all_edges);

    // igraph_vector_destroy(&edges);

    // for(int i=0; i<all_edges.size()/2; i++){
    //     std::cout<<all_edges[2*i]<<" -- "<<all_edges[2*i+1]<<endl;
    // }
    // std::cout<<all_edges.size()/2<<endl;

    // maximum spanning tree of clique tree is junction tree. Here, weight is 1/weight -> minumum spanning tree

    igraph_minimum_spanning_tree(&g_clique, &idx_m_spanning_tree_edge, &e_weight);

    //print_vector(&idx_m_spanning_tree_edge, stdout);

    int l = igraph_vector_size(&idx_m_spanning_tree_edge);

    std::vector<std::vector<int>> junction_edge(l);

    igraph_vector_t e_tree;
    igraph_vector_init(&e_tree, 2*l);

    for (int i = 0; i < l; i++) {
        int idx_edge = (long int)VECTOR(idx_m_spanning_tree_edge)[(long int)i];
        junction_edge[i] = idx_pairs_cliques[idx_edge];
        VECTOR(e_tree)[2*i] = idx_pairs_cliques[idx_edge][0];
        VECTOR(e_tree)[2*i+1] = idx_pairs_cliques[idx_edge][1];

        //cout<<idx_pairs_cliques[idx_edge][0]<<" "<<idx_pairs_cliques[idx_edge][1]<<endl;
    }

    ts.push_back(std::chrono::system_clock::now());
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

    igraph_topological_sorting(&g_tree, &order, IGRAPH_IN); // --------------------------------------> root selecting (some algorithm) and make tree
    
    ts.push_back(std::chrono::system_clock::now());

    // // check all edge
    // igraph_vector_t edges;
    // igraph_vector_init(&edges, 1);
    // igraph_get_edgelist(&g_tree, &edges, 0);

    // int nb_edges = igraph_vector_size(&edges);
    // std::vector<int> all_edges(nb_edges);
    // copy_vector(&edges, all_edges);

    // igraph_vector_destroy(&edges);

    // for(int i=0; i<all_edges.size()/2; i++){
    //     std::cout<<all_edges[2*i]<<" -- "<<all_edges[2*i+1]<<endl;
    // }
    // std::cout<<all_edges.size()/2<<endl;

    // int nb_order = igraph_vector_size(&order);
    // std::vector<int> orders(nb_order);
    // copy_vector(&order, orders);
    // print_std_vector(orders);

    //std::vector<int> vec_order(n);
    //copy_vector(&order, vec_order);
    //VECTOR(e)[2*i]
    //VECTOR(order)[i]
    igraph_vector_destroy(&e_clique);
    igraph_vector_destroy(&e_weight);
    igraph_vector_destroy(&idx_m_spanning_tree_edge);
    igraph_vector_destroy(&e_tree);
    igraph_destroy(&g_tree);
    igraph_destroy(&g_clique);

    jt_order.clear();

    ts.push_back(std::chrono::system_clock::now());
    std::vector<bool> used(this->factors.size(), false);

    for(int i=0; i<n; i++){
        for(int j=0; j<v_idxs_cliques[VECTOR(order)[i]].size(); j++){
            for(auto& nf : this->vars[v_idxs_cliques[VECTOR(order)[i]][j]]->neighbor_factors)
            {   
                if(!count(v_idxs_cliques[VECTOR(order)[i]].begin(), v_idxs_cliques[VECTOR(order)[i]].end(), nf->idx)){
                    {
                        //lock_guard<mutex> lock_guard(m);
                        if(!used[nf->idx]){
                            this->jt_order.emplace_back(nf->idx);
                            used[nf->idx] = true;
                        }
                    }
                }
            }
        }
    }
    igraph_vector_destroy(&order);
    
    ts.push_back(std::chrono::system_clock::now());
    // double tt = 0;
    // for(int i=0; i<ts.size()-1; i++){
    //     std::chrono::duration<double> sec = ts[i+1]-ts[i];
    //     std::cout<<sec.count()<<" ";
    //     tt += sec.count();
    // }
    // std::cout<<" "<<tt<<std::endl;

    


}

void FactorGraph::propagateMsgAll(bool JT, mutex& m){

    int cnt = 0;

    if(JT){
        for(int i=0; i<this->jt_order.size(); i++){
            try{
                this->propagateMsg(jt_order[i], m);
            }
            catch(...){cout<<"Unknown error happen"<<endl;}
                cnt += 1;
        }
    }

    else{
        int nb_factors = this->factors.size();

        auto randeng = std::default_random_engine {};
        std::vector<int> order(nb_factors);
        iota(order.begin(), order.end(), 0);
        shuffle(order.begin(), order.end(), randeng);

        // omp_set_num_threads(10);
        // #pragma omp parallel for
        for(int i=0; i<nb_factors; i++){
            this->propagateMsg(order[i], m);           
            //this->propagateMsg(i, m);           
        }
    }

    //std::cout<<this->factors.size()<<" "<<cnt<<std::endl;


}

void FactorGraph::updatePoseAll(){
    int nb_vars = this->vars.size();
    for(int i=0; i<nb_vars; i++){
            //cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Update x"<<i<<endl;
            this->updateVar(i);
    }
}